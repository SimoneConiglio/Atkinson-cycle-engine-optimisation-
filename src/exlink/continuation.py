"""Reaching the reliable region by walking to it.

The problem this solves
-----------------------
Section 3.10 records a failure that section 6.2 shows is expensive: from the
*deterministic* optimum, SLSQP with a reliability constraint attached returns
its starting point unchanged -- not only for a demanding target, but for
``beta >= -0.2``, a step of 0.17 from where it stands.  Meanwhile sampling in a
neighbourhood of that same point finds fully feasible designs that halve the
failure probability at no cost in range.  Better designs are there; descent
cannot get to them.

The cause is not the reliability gradient, which is smooth.  It is the
thinness of the feasible set: a step of 0.05 mm along the normalised
``grad beta`` leaves the geometric constraint set entirely, so every direction
that improves reliability also violates something, and the QP subproblem has
nothing admissible to offer.

Why a continuation is the right instrument
------------------------------------------
The obstruction is *local*, and it is local to the target rather than to the
problem: at ``beta_target = 3`` nothing near the start is admissible, so the
line search fails at the first step; at ``beta_target = -0.3`` the current
point is admissible and there is somewhere to go.  Solving a sequence of
problems whose targets rise from where the design already stands, warm-starting
each from the last, replaces one impossible step with many possible ones.

That is a homotopy in the constraint level rather than in the objective, and
the property it relies on is that the reliable set is *connected to* the
starting point through the feasible set even though it is not reachable from it
in one step.  Whether that holds is exactly what this module measures; it is
not assumed.

What it does not do
-------------------
A continuation cannot cross a gap.  If the admissible set at some rung is
disconnected from the current point, the rung fails and every rung above it
fails with it, and the run says so rather than reporting the last success as
the answer.  Section 6.8 reports which rungs moved.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import GLOBAL_BOUNDS, Bounds, Design
from .robustness import failure_probability
from .synthesis import RangeFit, TargetMotion, maximise_range_from_target

DEFAULT_SCHEDULE: tuple[float, ...] = (0.0, 1.0, 2.0, 2.5, 3.0)
"""Reliability targets to climb, in order.

Coarse at the bottom, where the constraint is slack and each solve is a plain
range maximisation, and fine at the top, where it binds and the admissible set
is thin.  The first rung is deliberately below the target most starts already
meet, so that the first solve is a warm-up rather than a test.
"""

MOVED_TOLERANCE = 1.0e-9
"""Relative design change below which a rung counts as not having moved."""

CONSTRAINT_TOLERANCE = 1.0e-6
"""Largest constraint violation a rung may carry and still count as feasible.

SLSQP converges *onto* its active constraints, so a solution routinely sits a
few parts in ten million outside one of them.  Calling that a failure would
report every successful rung as a failure."""


@dataclass(frozen=True)
class ContinuationStep:
    """One rung of the ladder."""

    beta_target: float
    design: Design
    km_per_litre: float

    steered_beta: float
    """``min_i beta_i``, the quantity the solve actually holds (§3.10)."""

    system_beta: float
    """The system index, reported rather than steered."""

    system_pf: float
    worst_constraint: float
    """Largest constraint value; ``<= 0`` means every constraint holds."""

    feasible: bool
    """Whether every constraint holds *at the bands this solve is using*.

    Not the same as :attr:`strictly_feasible`, and the distinction matters:
    §6.2 settles the study on a widened specification, and a rung solved
    against a band of 0.15 must be judged against 0.15.  Scoring it against
    the specification as written would mark every rung a failure and say
    nothing about whether the continuation works.
    """

    strictly_feasible: bool
    """Whether it would also pass the specification as written."""

    moved: bool
    """Whether this rung changed the design at all.

    A rung that does not move is the failure mode this module exists to
    diagnose: it means the line search found no admissible direction, which is
    what happened to the single-shot solve at every target.
    """

    evaluations: int
    converged: bool


@dataclass(frozen=True)
class Continuation:
    """A ladder of reliability targets, climbed or not."""

    start: Design
    steps: tuple[ContinuationStep, ...]

    @property
    def reached(self) -> bool:
        """Whether the last rung was met, feasibly."""
        if not self.steps:
            return False
        last = self.steps[-1]
        return last.feasible and last.steered_beta >= last.beta_target - 1.0e-6

    @property
    def best(self) -> ContinuationStep | None:
        """The highest rung that was met, feasibly."""
        met = [
            step
            for step in self.steps
            if step.feasible and step.steered_beta >= step.beta_target - 1.0e-6
        ]
        return met[-1] if met else None

    @property
    def evaluations(self) -> int:
        """Total analyses across every rung."""
        return sum(step.evaluations for step in self.steps)


def reliability_continuation(
    target: TargetMotion,
    start: Design,
    schedule: tuple[float, ...] = DEFAULT_SCHEDULE,
    speed_rpm: float = 1000.0,
    bounds: Bounds = GLOBAL_BOUNDS,
    max_iterations: int = 60,
    band: float = 0.15,
    grade: int = 8,
    module: float | None = None,
    teeth: int | None = None,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
) -> Continuation:
    """Climb a ladder of reliability targets, warm-starting each rung.

    Args:
        target: The fallback motion, as §3.10's ladder needs.
        start: Where to begin -- typically a deterministic optimum.
        schedule: Reliability targets, ascending.
        speed_rpm: Analysis speed [rev/min].
        bounds: The design box.
        max_iterations: SLSQP budget *per rung*, not in total.
        band: Half-width of the two relaxed equalities.
        grade: ISO 286 grade the probabilities assume.
        module: Gear module to pin [mm].
        teeth: Tooth count to pin.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        Every rung attempted, in order, whether or not it succeeded.  A rung
        that fails does not stop the climb: the next is attempted from the same
        point, because a target that is unreachable in one step is sometimes
        reachable after the one above it has been tried and the design has
        drifted.  What it does mean is that :attr:`Continuation.reached` is
        false unless the *last* rung was met.
    """
    current = start
    steps: list[ContinuationStep] = []
    for level in schedule:
        fit = maximise_range_from_target(
            target,
            current,
            speed_rpm=speed_rpm,
            bounds=bounds,
            max_iterations=max_iterations,
            band=band,
            beta_target=level,
            grade=grade,
            module=module,
            teeth=teeth,
            targets=targets,
            spec=spec,
        )
        if fit is None:
            steps.append(_failed_step(level, current))
            continue
        step = _step(level, fit, current, band, grade, targets, spec)
        steps.append(step)
        # Only a rung that actually improved is worth carrying forward: a rung
        # that walked out of the feasible set would poison every rung above it.
        if step.feasible:
            current = step.design
    return Continuation(start=start, steps=tuple(steps))


def _step(
    level: float,
    fit: RangeFit,
    previous: Design,
    band: float,
    grade: int,
    targets: DesignTargets,
    spec: EngineSpec,
) -> ContinuationStep:
    """Assemble one rung's record, measuring the reliability it reached."""
    reliability = failure_probability(
        fit.design,
        grade=grade,
        targets=targets,
        spec=spec,
        band={"expansion_stroke": band, "compression_ratio": band},
    )
    steered = float("-inf") if reliability is None else float(np.min(reliability.moments.beta))
    moved = not np.allclose(
        fit.design.to_array(), previous.to_array(), rtol=MOVED_TOLERANCE, atol=0.0
    )
    return ContinuationStep(
        beta_target=float(level),
        design=fit.design,
        km_per_litre=fit.km_per_litre,
        steered_beta=steered,
        system_beta=float("-inf") if reliability is None else reliability.system_beta,
        system_pf=1.0 if reliability is None else reliability.system,
        worst_constraint=fit.worst_constraint,
        feasible=fit.worst_constraint <= CONSTRAINT_TOLERANCE,
        strictly_feasible=fit.feasible,
        moved=moved,
        evaluations=fit.evaluations,
        converged=fit.converged,
    )


def _failed_step(level: float, current: Design) -> ContinuationStep:
    """A rung whose solve produced nothing analysable."""
    return ContinuationStep(
        beta_target=float(level),
        design=current,
        km_per_litre=0.0,
        steered_beta=float("-inf"),
        system_beta=float("-inf"),
        system_pf=1.0,
        worst_constraint=float("inf"),
        feasible=False,
        strictly_feasible=False,
        moved=False,
        evaluations=0,
        converged=False,
    )


def _verdict(step: ContinuationStep) -> str:
    """One word for how a rung ended."""
    if not step.feasible:
        return "infeasible"
    if step.steered_beta < step.beta_target - 1.0e-6:
        return "short"
    return "met, strict" if step.strictly_feasible else "met"


def format_continuation(continuation: Continuation) -> str:
    """Render a :class:`Continuation` as an aligned table."""
    lines = [
        "reliability continuation",
        "=" * 72,
        "",
        f"  {'target':>7}{'km/L':>10}{'min beta':>11}{'sys beta':>10}"
        f"{'P_f':>11}{'worst g':>11}{'evals':>7}  moved  verdict",
    ]
    for step in continuation.steps:
        lines.append(
            f"  {step.beta_target:>7.2f}{step.km_per_litre:>10.1f}"
            f"{step.steered_beta:>11.3f}{step.system_beta:>10.2f}"
            f"{step.system_pf:>11.3e}{step.worst_constraint:>11.2e}"
            f"{step.evaluations:>7}  {'yes' if step.moved else ' NO':>3}"
            f"  {_verdict(step)}"
        )
    lines.append("")
    best = continuation.best
    if best is None:
        lines.append("  no rung was met feasibly")
    else:
        lines.append(
            f"  highest rung met: beta >= {best.beta_target:.2f} at "
            f"{best.km_per_litre:.1f} km/L, P_f = {best.system_pf:.3e}"
        )
    lines.append(f"  total evaluations {continuation.evaluations}")
    return "\n".join(lines)
