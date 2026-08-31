"""How strong is the coupling, and does the formulation choice matter?

Two questions a reviewer of any MDO study asks, and neither is answered by
asserting that a problem is "strongly coupled".

How strong is the coupling?
---------------------------
The sizing/dynamics loop is solved by Gauss-Seidel, and a Gauss-Seidel
iteration converges linearly with a factor that *is* the coupling strength:

.. math:: \\rho = \\lim_{k \\to \\infty}
    \\frac{\\|r_{k+1}\\|}{\\|r_k\\|}

``rho`` near 0 means the disciplines barely see each other and a single
sequential pass would have done; ``rho`` near 1 means each discipline
substantially rewrites the other's input and the fixed point is doing real
work; ``rho >= 1`` means no fixed point exists at all -- the mechanism cannot
be built to survive its own inertia.

:func:`coupling_strength` measures it from the residual history the solve
already records, at no extra cost.  Because the loop gain scales with the
inertia forces, ``rho`` is a function of engine speed, and reporting
``rho(omega)`` is the honest way to say how coupled this problem is: not very,
at rest; decisively, at speed.

Does the formulation matter?
-----------------------------
MDF converges the MDA at every optimizer iteration, so every point it evaluates
is physically consistent, and the optimizer sees a small design space (the 11
linkage variables).  IDF hands the coupling variables to the optimizer as extra
design variables with consistency constraints, so no inner iteration runs at
all, but the design space grows by the 7 diameters and 7 equality constraints
come with them.

The textbook trade is that IDF wins when the MDA is expensive and loses when
the coupling variables are numerous.  Here the MDA is expensive *and* the
coupling variables are numerous, so the answer is not obvious from theory and
has to be measured.  :func:`compare_formulations` measures it.

One asymmetry is worth stating in advance, because it is a property of this
problem rather than of the formulations in general: intermediate IDF iterates
are *not* physically consistent, so any quantity read from one -- the range,
the mass budget, the friction -- is meaningless until the consistency
constraints are satisfied.  With MDF, every iterate is a real engine.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any

import numpy as np

from .coupled import CoupledResult, solve_for_design
from .design import Design
from .dynamics import DEFAULT_SPEED_RPM


@dataclass(frozen=True)
class CouplingStrength:
    """Observed contraction of the sizing/dynamics fixed point."""

    speed_rpm: float
    rho: float
    """Asymptotic residual ratio per sweep; the coupling strength."""

    sweeps: int
    """Sweeps the solve took."""

    converged: bool
    history: list[float] = field(default_factory=list)

    @property
    def sweeps_per_decade(self) -> float:
        """How many sweeps the loop needs to gain one digit.

        A compact way to read ``rho``: 1 sweep per decade is a nearly decoupled
        problem, 20 is a strongly coupled one.
        """
        if not 0.0 < self.rho < 1.0:
            return float("inf")
        return -1.0 / np.log10(self.rho)

    @property
    def descriptor(self) -> str:
        """A word for how coupled this is, on the usual reading of ``rho``."""
        if self.rho >= 1.0 or not self.converged:
            return "divergent"
        if self.rho < 0.1:
            return "weak"
        if self.rho < 0.5:
            return "moderate"
        return "strong"


def coupling_strength(
    design: Design,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    tail: int = 5,
    **kwargs: Any,
) -> CouplingStrength:
    """Measure the coupling strength from the fixed point's own residuals.

    The ratio is taken over the last few sweeps, where the iteration has
    settled into its asymptotic linear rate; the early sweeps are contaminated
    by the starting guess and would understate it.

    Args:
        design: The mechanism dimensions.
        speed_rpm: Crankshaft speed [rev/min].
        tail: Sweeps at the end of the history to average the ratio over.
        **kwargs: Forwarded to :func:`~exlink.coupled.solve_for_design`.

    Returns:
        The measurement.  A design whose loop runs away comes back with
        ``rho >= 1`` rather than raising.
    """
    try:
        result: CoupledResult = solve_for_design(design, speed_rpm=speed_rpm, **kwargs)
    except ValueError:
        return CouplingStrength(
            speed_rpm=speed_rpm, rho=float("inf"), sweeps=0, converged=False
        )

    history = [value for value in result.history if np.isfinite(value) and value > 0.0]
    if len(history) < 3:
        # Converged before a rate could be observed: the disciplines are
        # effectively decoupled at this speed.
        return CouplingStrength(
            speed_rpm=speed_rpm,
            rho=0.0,
            sweeps=result.iterations,
            converged=result.converged,
            history=list(result.history),
        )

    window = history[-min(tail + 1, len(history)) :]
    ratios = [later / earlier for earlier, later in pairwise(window)]
    rho = float(np.exp(np.mean(np.log(np.clip(ratios, 1.0e-16, None)))))
    return CouplingStrength(
        speed_rpm=speed_rpm,
        rho=rho,
        sweeps=result.iterations,
        converged=result.converged and not result.saturated,
        history=list(result.history),
    )


def coupling_curve(
    design: Design,
    speeds: tuple[float, ...] = (0.0, 250.0, 500.0, 750.0, 1000.0, 1250.0, 1500.0),
    **kwargs: Any,
) -> list[CouplingStrength]:
    """Coupling strength against engine speed.

    The point of the curve is that the answer to "is this problem coupled?"
    is *it depends on the operating point*, and the dependence is steep.  At
    rest the loop gain is zero -- there are no inertia forces, so the sections
    do not feed back into the loads at all, and the quasi-static problem of the
    geometric formulation is recovered exactly.  The gain grows with the
    inertia, i.e. with ``omega^2``.

    Args:
        design: The mechanism dimensions.
        speeds: Crankshaft speeds [rev/min].
        **kwargs: Forwarded to :func:`coupling_strength`.

    Returns:
        One measurement per speed.
    """
    return [coupling_strength(design, speed_rpm=speed, **kwargs) for speed in speeds]


@dataclass(frozen=True)
class FormulationResult:
    """What one formulation cost and what it found."""

    name: str
    objective: float
    """Best objective value reached."""

    design: Design | None
    evaluations: int
    """Objective evaluations the optimizer performed."""

    disciplinary_calls: int
    """Total discipline executions, MDA sweeps included.

    The fair cost measure: MDF hides its work inside the MDA, and counting only
    optimizer iterations would flatter it.
    """

    seconds: float
    feasible: bool
    error: str = ""


def compare_formulations(
    initial: Design | None = None,
    objective: str = "total_mass",
    speed_rpm: float = DEFAULT_SPEED_RPM,
    max_iter: int = 25,
    relative: float = 0.20,
    algorithm: str = "SLSQP",
    names: tuple[str, ...] = ("MDF", "IDF"),
    **kwargs: Any,
) -> list[FormulationResult]:
    """Run the same coupled problem under each formulation and compare.

    Args:
        initial: Starting design; defaults to the refined reference.
        objective: Output to minimise.
        speed_rpm: Crankshaft speed [rev/min].
        max_iter: Evaluation budget for each run.
        relative: Half-width of the design box, as a fraction of ``|X_0|``.
        algorithm: A single-objective GEMSEO optimizer.
        names: Formulations to compare.
        **kwargs: Forwarded to the scenario builder.

    Returns:
        One result per formulation, in the order given.  A formulation that
        fails outright is reported as a result carrying its error, not raised:
        a failure is a comparison outcome.
    """
    from .design import Bounds
    from .model import analyse
    from .reference import REFINED_DESIGN
    from .scenarios import COUPLED_SAMPLES, _best_design, build_coupled_scenario, is_feasible

    start = REFINED_DESIGN if initial is None else initial
    box = Bounds.around(start, relative=relative)

    results: list[FormulationResult] = []
    for name in names:
        began = time.perf_counter()
        try:
            scenario = build_coupled_scenario(
                objective,
                bounds=box,
                initial=start,
                speed_rpm=speed_rpm,
                formulation_name=name,
                **kwargs,
            )
            scenario.execute(algo_name=algorithm, max_iter=max_iter)
            design = _best_design(scenario, objective=objective)
            problem = scenario.formulation.optimization_problem
            calls = sum(
                discipline.execution_statistics.n_executions or 0
                for discipline in scenario.disciplines
            )
            results.append(
                FormulationResult(
                    name=name,
                    objective=float(np.min(problem.database.get_function_history(objective))),
                    design=design,
                    evaluations=len(problem.database),
                    disciplinary_calls=int(calls),
                    seconds=time.perf_counter() - began,
                    feasible=is_feasible(analyse(design, samples=COUPLED_SAMPLES)),
                )
            )
        except Exception as error:
            results.append(
                FormulationResult(
                    name=name,
                    objective=float("nan"),
                    design=None,
                    evaluations=0,
                    disciplinary_calls=0,
                    seconds=time.perf_counter() - began,
                    feasible=False,
                    error=f"{type(error).__name__}: {error}",
                )
            )
    return results


def format_coupling(measurements: list[CouplingStrength]) -> str:
    """Render a coupling curve as an aligned table."""
    lines = ["coupling strength against speed", "=" * 31, ""]
    lines.append(f"  {'rpm':>7}{'rho':>10}{'sweeps':>9}{'per decade':>13}   verdict")
    for item in measurements:
        per_decade = item.sweeps_per_decade
        rendered = "inf" if not np.isfinite(per_decade) else f"{per_decade:.1f}"
        lines.append(
            f"  {item.speed_rpm:>7.0f}{item.rho:>10.4f}{item.sweeps:>9}"
            f"{rendered:>13}   {item.descriptor}"
        )
    return "\n".join(lines)
