"""Getting inside the feasible set before trying to optimise within it.

The measurement this answers
----------------------------
§6.9 reports that manifold-projected restarts reached feasibility in **0 of 6
attempts** on the range problem at an affordable budget, which is why §3.9's
multistart is inconclusive there: it is not that the restarts found worse
optima, it is that they never found the feasible set at all, so there was
nothing to compare.  §7.3 asks for a restoration phase before each restart.

The reason a restart lands outside is §3.4's: the relaxed equalities make the
feasible set a thin shell, and a design scattered from a good one by any
useful amount is outside it.  A range-maximising solve started there spends its
whole budget on the constraint violation and has none left for the objective --
and SLSQP in particular, which needs a feasible point to build a useful QP
around, frequently never recovers.

The restoration problem
-----------------------
Separate the two jobs.  Before optimising, solve

.. math:: \\max_{X,\\,t} \\; t \\quad \\text{s.t.} \\quad c_i(X) \\ge t, \\;
          X \\in [X_{lb}, X_{ub}]

with :math:`c` the constraints on the "positive means satisfied" convention of
:func:`exlink.synthesis._range_constraints`.  The solution maximises the
*smallest* margin, so it does not merely reach the feasible set but heads for
its interior, which is where the next solve wants to start.

Why the epigraph form rather than a penalty
-------------------------------------------
The natural objective is :math:`\\min_i c_i(X)`, which is non-smooth wherever
two constraints tie -- and on this problem they tie constantly, because the two
sides of each relaxed band are the same residual negated.  Introducing the
scalar :math:`t` and moving the minimum into the constraints makes an
eleven-variable non-smooth problem into a twelve-variable smooth one, at no
cost in the answer.  A weighted penalty would also be smooth but would need
weights, and the weights are exactly what is not known before the run.

What restoration cannot do
--------------------------
It maximises a margin, so it stops at the interior point nearest its start.  It
does not search for the *global* feasible region, and if a restart is in a
basin whose closest feasible point is poor, restoration will hand the optimizer
a poor start -- promptly, rather than after spending the whole budget failing.
That is still the improvement worth having: it turns a run that reports nothing
into one that reports something.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design
from .performance import evaluate
from .synthesis import NUMBER_OF_CONSTRAINTS, TargetMotion, _range_constraints

FloatArray = NDArray[np.float64]

DEFAULT_MARGIN = 1.0e-4
"""Smallest margin counted as being *inside* rather than merely on the boundary.

Zero would accept a point sitting exactly on a constraint, which is the state
§6.2 shows fails half the builds and §6.8 shows breaks the reliability
estimate.  A restoration that stops at the boundary has not done its job.
"""


@dataclass(frozen=True)
class Restoration:
    """One restoration solve."""

    design: Design
    margin: float
    """The smallest constraint margin reached; positive means feasible."""

    start_margin: float
    """What it was before, for comparison."""

    evaluations: int
    converged: bool

    @property
    def feasible(self) -> bool:
        """Whether the point ended inside the feasible set."""
        return self.margin >= DEFAULT_MARGIN

    @property
    def improved(self) -> bool:
        """Whether the solve moved the margin in the right direction at all."""
        return self.margin > self.start_margin


def restore(
    target: TargetMotion,
    start: Design,
    bounds: Bounds = GLOBAL_BOUNDS,
    speed_rpm: float = 1000.0,
    band: float = 0.15,
    max_iterations: int = 60,
    module: float | None = None,
    teeth: int | None = None,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
) -> Restoration:
    """Push a design into the interior of the feasible set.

    Args:
        target: Supplies the two equality targets the bands are taken about.
        start: The design to restore, feasible or not.
        bounds: The design box.
        speed_rpm: Analysis speed [rev/min].
        band: Half-width of the two relaxed equalities.
        max_iterations: SLSQP budget.
        module: Gear module to pin [mm].
        teeth: Tooth count to pin.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        The restored design and the margin it reached.
    """
    from scipy.optimize import minimize

    calls = {"n": 0}
    cache: dict[bytes, FloatArray] = {}

    def margins(vector: FloatArray) -> FloatArray:
        key = np.ascontiguousarray(vector, dtype=float).tobytes()
        hit = cache.get(key)
        if hit is not None:
            return hit
        calls["n"] += 1
        try:
            performance = evaluate(
                Design.from_array(vector),
                speed_rpm=speed_rpm,
                samples=target.samples,
                module=module,
                teeth=teeth,
                spec=spec,
            )
            rows = _range_constraints(performance, target, band, targets)
        except (ValueError, FloatingPointError):
            rows = np.full(NUMBER_OF_CONSTRAINTS, -1.0e3)
        if len(cache) >= 4 * len(VARIABLE_NAMES):
            cache.pop(next(iter(cache)))
        cache[key] = rows
        return rows

    lower = np.array(bounds.lower, dtype=float)
    upper = np.array(bounds.upper, dtype=float)
    begin = np.clip(start.to_array(), lower, upper)
    if module is not None and teeth is not None:
        from .gears import lattice_inter_axle

        index = VARIABLE_NAMES.index("I")
        pinned = float(lattice_inter_axle(module, teeth))
        lower[index] = upper[index] = pinned
        begin[index] = pinned

    start_margin = float(np.min(margins(begin)))
    # The epigraph variable rides along as the twelfth unknown, bounded by the
    # worst margin any design in the box could plausibly show and the best.
    augmented = np.concatenate([begin, [start_margin]])
    box = [*zip(lower, upper, strict=True), (-1.0e3, 1.0e3)]

    outcome = minimize(
        lambda z: -z[-1],
        augmented,
        method="SLSQP",
        bounds=box,
        jac=lambda z: np.concatenate([np.zeros(z.size - 1), [-1.0]]),
        constraints=[{"type": "ineq", "fun": lambda z: margins(z[:-1]) - z[-1]}],
        options={"maxiter": int(max_iterations), "ftol": 1.0e-10},
    )
    design = Design.from_array(np.asarray(outcome.x[:-1], dtype=float))
    return Restoration(
        design=design,
        margin=float(np.min(margins(design.to_array()))),
        start_margin=start_margin,
        evaluations=calls["n"],
        converged=bool(outcome.success),
    )


def format_restoration(results: list[Restoration]) -> str:
    """Render a set of restorations as an aligned table."""
    lines = ["feasibility restoration", "=" * 56, ""]
    lines.append(f"  {'start margin':>14}{'end margin':>14}{'evals':>8}  verdict")
    for item in results:
        verdict = "feasible" if item.feasible else ("improved" if item.improved else "stuck")
        lines.append(
            f"  {item.start_margin:>14.3e}{item.margin:>14.3e}{item.evaluations:>8}  {verdict}"
        )
    restored = sum(item.feasible for item in results)
    lines.append("")
    lines.append(f"  {restored} of {len(results)} reached the feasible set")
    return "\n".join(lines)
