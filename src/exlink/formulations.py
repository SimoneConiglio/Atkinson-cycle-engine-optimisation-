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
the coupling variables are numerous.  Here the MDA is expensive, so the trade
looks open -- until the coupling is counted.

It is not close.  The strong couplings are ``diameters`` and ``piston_mass``,
which are small, and ``member_axial`` and ``member_bending``, which are not
scalars at all: they are the internal load history of every member, at every
crank angle, at every station along it.  That is **45 367 coupling scalars
against 11 design variables**.  IDF would carry all of them in the design space
with a matching consistency constraint each, to optimise eleven real degrees of
freedom.

So the formulation question has a decisive answer on this problem, and it is
structural rather than a matter of timing: IDF is not slower here, it is
unavailable.  :func:`coupling_dimension` reports the count and
:func:`compare_formulations` runs both anyway, because an attempt that fails
for a stated reason is a better answer than an assertion.

The general lesson is worth separating from this mechanism.  IDF's cost scales
with the *dimension* of the coupling, and a discipline pair that exchanges
distributed fields -- load histories, pressure distributions, temperature
fields -- rather than a handful of scalars will always sit on the wrong side of
that trade, however expensive its MDA.

One asymmetry is worth stating in advance, because it is a property of this
problem rather than of the formulations in general: intermediate IDF iterates
are *not* physically consistent, so any quantity read from one -- the range,
the mass budget, the friction -- is meaningless until the consistency
constraints are satisfied.  With MDF, every iterate is a real engine.
"""

from __future__ import annotations

import time
from collections.abc import Sequence
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


def coupling_dimension(
    samples: int | None = None, stations: int | None = None
) -> dict[str, int]:
    """How many scalars the two coupled disciplines exchange.

    This is the number that decides the formulation question, and it is not
    close.  The strong couplings here are ``diameters``, ``piston_mass`` and --
    dominating everything -- ``member_axial`` and ``member_bending``, which are
    not scalars but the *internal load history of every member, at every crank
    angle, at every station along it*.

    IDF puts every coupling variable into the design space with a matching
    consistency constraint.  With 11 real degrees of freedom and tens of
    thousands of coupling scalars, that trade is not merely unfavourable, it is
    unavailable: the optimizer would carry four orders of magnitude more
    variables than the problem has.

    Args:
        samples: Crank angles per revolution; the package default if omitted.
        stations: Sections along each member; the package default if omitted.

    Returns:
        ``{variable: size}`` plus ``"total"`` and ``"design_variables"``.
    """
    from .design import VARIABLE_NAMES
    from .disciplines import COUPLED_SAMPLES
    from .dynamics import MEMBER_NAMES
    from .sizing import STATIONS

    n_angles = COUPLED_SAMPLES if samples is None else samples
    n_stations = STATIONS if stations is None else stations
    history = len(MEMBER_NAMES) * n_angles * n_stations
    sizes = {
        "member_axial": history,
        "member_bending": history,
        "diameters": len(MEMBER_NAMES),
        "piston_mass": 1,
    }
    sizes["total"] = sum(sizes.values())
    sizes["design_variables"] = len(VARIABLE_NAMES)
    return sizes


def motion_harmonics(
    design: Design,
    tolerances: Sequence[float] = (0.1, 0.01, 0.001),
    samples: int = 720,
) -> dict[float, int]:
    """How many Fourier harmonics reproduce the piston motion, per RMS tolerance.

    Why this is the interesting number
    ----------------------------------
    :func:`coupling_dimension` counts the coupling *pointwise* and concludes
    that IDF is unavailable.  That conclusion is a property of the
    parameterisation, not of the physics, and this function is the measurement
    that shows it.

    The linkage's eleven dimensions reach the rest of the model through exactly
    one quantity: the piston motion ``lam(theta)``.  The cycle turns it into a
    volume and a pressure, the dynamics differentiate it twice for the inertia
    loads, the vehicle sees only the work that results.  So ``lam`` *is* the
    coupling variable of the whole problem -- and it is smooth and periodic, so
    counting it at 45 368 grid points overstates it by orders of magnitude.

    In a Fourier basis the count is what this returns: of order twenty
    coefficients to reproduce the motion more tightly than the part can be
    machined.  A decomposition on that basis would carry tens of consistency
    variables where the pointwise one carries tens of thousands.  See §7.4 of
    the documentation.

    Args:
        design: The mechanism.
        tolerances: RMS errors to report, in millimetres, largest first.
        samples: Crank angles per revolution to analyse.

    Returns:
        ``{tolerance: harmonics}``: the fewest leading harmonics whose truncation
        error is below that RMS, one entry per requested tolerance.

    Raises:
        ValueError: If the design cannot be analysed.
    """
    from .model import analyse

    analysis = analyse(design, samples=samples)
    if not analysis.valid:
        msg = f"cannot analyse the motion of an unanalysable design: {analysis.metrics.reason}"
        raise ValueError(msg)
    lam = np.asarray(analysis.require_solved().kinematics.lam, dtype=float)
    # Parseval: the power in harmonic k, with the one-sided factor of two on
    # everything but the mean, sums to the mean square of ``lam``.
    power = np.abs(np.fft.rfft(lam) / lam.size) ** 2
    power[1:] *= 2.0
    # Truncation error after keeping the first k harmonics.
    residual = np.sqrt(np.maximum(power.sum() - np.cumsum(power), 0.0))
    counts: dict[float, int] = {}
    for tolerance in tolerances:
        below = np.flatnonzero(residual < float(tolerance))
        counts[float(tolerance)] = int(below[0]) if below.size else int(residual.size)
    return counts


def format_formulations(rows: list[FormulationResult]) -> str:
    """Render a formulation comparison, with the coupling dimension that explains it."""
    sizes = coupling_dimension()
    lines = ["formulation comparison", "=" * 22, ""]
    lines.append(f"  {'name':<6}{'objective':>12}{'evals':>8}{'seconds':>10}{'feasible':>10}")
    for row in rows:
        objective = "n/a" if not np.isfinite(row.objective) else f"{row.objective:.4f}"
        lines.append(
            f"  {row.name:<6}{objective:>12}{row.evaluations:>8}"
            f"{row.seconds:>10.1f}{row.feasible!s:>10}"
        )
        if row.error:
            lines.append(f"      {row.error[:100]}")
    lines.append("")
    lines.append("  why, in one number:")
    lines.append(f"    design variables      {sizes['design_variables']:>8}")
    lines.append(f"    coupling variables    {sizes['total']:>8}")
    lines.append(
        f"      of which load histories {2 * sizes['member_axial']:>6}"
        f"  ({sizes['member_axial'] // 1} each)"
    )
    return "\n".join(lines)
