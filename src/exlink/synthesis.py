"""Prescribed-motion synthesis, and the range problem it grew into.

The formulation of :mod:`exlink.scenarios` asks for the best mechanism and lets
the piston motion fall out of it.  Kinematic synthesis conventionally does the
opposite -- it *prescribes* a target motion and fits the linkage to it:

.. math::
    \\min_X\\; J(X) = \\sum_k \\bigl(\\lambda_k(X) - \\lambda^\\star_k\\bigr)^2

This module started there and did not stay there, because each round of the
exercise found the same thing: whatever constraint had been left out of the fit
was the one the solve then violated.  The entry points below are that argument
in order.

Why a target helps at all
-------------------------
Both equality requirements are functionals of :math:`\\lambda` alone --
:math:`STE` is the span from top dead centre to the deeper bottom dead centre,
and :math:`\\varepsilon` follows from the shallower one through the swept
volume.  So a target built to have the right two strokes satisfies both
requirements before any linkage exists.  :func:`target_motion` solves for one
and :func:`describe_target` confirms it with the code the constraints use.

That makes the fit a sampling device.  Uniform sampling of the design box finds
a design on the equality manifold in 0 of 12 000 attempts, because the manifold
has measure zero; sampling *and then fitting* reaches it in most attempts.

The entry points
----------------
:func:`fit_to_target`
    Unconstrained least squares.  Enough only when the target is reachable:
    against one that is not, it walks ten to fifteen units outside the geometric
    constraint set chasing a motion no buildable linkage produces.

:func:`fit_within_constraints`
    The same fit subject to ``g(X) <= 0``, and with ``hold_bands`` also subject
    to the tolerance bands.  Absorbing the equalities into the target removed
    the *measure-zero* part of the feasible set, which is the only part the
    reformulation earns the right to drop; the inequalities define a
    full-dimensional set and cost one SQP instead of one least-squares solve.

:func:`maximise_range_from_target`
    The endpoint, and no longer a synthesis problem: maximise the quantity the
    application scores, subject to all fourteen constraints, with the motion
    residual standing in wherever the range does not exist.  Its docstring has
    the objective ladder and why a constant penalty will not do.

What the target does *not* buy
------------------------------
Reachability is a property of the mechanism, not of the target.  A target can
sit exactly on the equality manifold and still be unfittable, which is why
:func:`target_from_design` exists: the obvious two-harmonic construction is one
of those.  And fitting to a single target returns a single linkage whatever
start it is given, so diversity has to come from varying the target rather than
the start.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import fsolve, least_squares

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .cycle import PhaseError, find_phases
from .design import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design
from .materials import FloatArray
from .model import INEQUALITY_NAMES, analyse, inequality_constraints

if TYPE_CHECKING:
    from .performance import Performance

DEFAULT_TARGET_SAMPLES = 360
"""Crank angles used to represent a target motion."""


def _two_harmonic(amplitude: float, skew: float, samples: int, offset: float) -> FloatArray:
    """A motion with two equal maxima and two unequal minima per revolution.

    The smallest basis that can produce an Atkinson motion at all.  The second
    harmonic supplies the two up-and-downs; the first, taken as a sine, is zero
    at both maxima and opposite at the two minima, so it deepens one bottom dead
    centre and lifts the other without separating the tops.  The top-dead-centre
    gap is therefore identically zero, which is what a target should ask for.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    return offset + amplitude * np.cos(2.0 * theta) + skew * np.sin(theta)


@dataclass(frozen=True)
class TargetMotion:
    """A piston motion prescribed by the two strokes it must realise."""

    lam: FloatArray
    """The target :math:`\\lambda^\\star(\\theta)` [mm]."""

    expansion_stroke: float
    """``STE`` the target was solved for [mm]."""

    compression_ratio: float
    """``epsilon`` the target was solved for."""

    amplitude: float
    skew: float

    @property
    def samples(self) -> int:
        """Crank angles in the target."""
        return int(self.lam.size)


def target_motion(
    expansion_stroke: float = 74.0,
    compression_ratio: float = 16.0,
    samples: int = DEFAULT_TARGET_SAMPLES,
    offset: float = 150.0,
    spec: EngineSpec = DEFAULT_SPEC,
) -> TargetMotion:
    """Build a motion that realises the two equality requirements exactly.

    The compression ratio fixes the compression stroke through the swept volume,
    ``STC = (epsilon - 1) V_0 / A_p``; the two strokes then fix the two
    harmonic coefficients, and the solve is two equations in two unknowns.

    Args:
        expansion_stroke: Required ``STE`` [mm].
        compression_ratio: Required ``epsilon``.
        samples: Crank angles to represent the target on.
        offset: Mean height; irrelevant to both requirements, and removed
            before any comparison with a real mechanism.
        spec: Fixed engine data.

    Returns:
        The target.

    Raises:
        ValueError: If no two-harmonic motion realises the pair.
    """
    stroke = (float(compression_ratio) - 1.0) * spec.dead_volume / spec.piston_area

    def residual(params: FloatArray) -> list[float]:
        phases = find_phases(_two_harmonic(params[0], params[1], samples, offset))
        return [
            phases.expansion_stroke - float(expansion_stroke),
            phases.compression_stroke - stroke,
        ]

    guess = np.array(
        [0.25 * (float(expansion_stroke) + stroke), 0.5 * (float(expansion_stroke) - stroke)]
    )
    solution, _info, status, message = fsolve(residual, guess, full_output=True)
    if status != 1:
        msg = (
            f"no two-harmonic motion realises STE={expansion_stroke}, "
            f"eps={compression_ratio}: {message}"
        )
        raise ValueError(msg)
    return TargetMotion(
        lam=_two_harmonic(float(solution[0]), float(solution[1]), samples, offset),
        expansion_stroke=float(expansion_stroke),
        compression_ratio=float(compression_ratio),
        amplitude=float(solution[0]),
        skew=float(solution[1]),
    )


def target_from_design(
    design: Design,
    expansion_stroke: float = 74.0,
    compression_ratio: float = 16.0,
    samples: int = DEFAULT_TARGET_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
) -> TargetMotion:
    """Build a target by correcting a real mechanism's motion onto the manifold.

    Why this exists, and why :func:`target_motion` is not enough
    ------------------------------------------------------------
    A two-harmonic target satisfies both equalities exactly but is not
    *reachable*: measured, the closest a seven-bar linkage gets to one is
    1.16 mm RMS, and that residual carries the fitted design outside both
    tolerance bands.  That is the attainability problem in its sharpest form --
    a target can be perfectly on the manifold and still be useless, because
    being on the manifold is a property of the target and being fittable is a
    property of the mechanism.

    So this starts from a motion that is reachable *by definition* -- one an
    actual design produces -- and adds only the two harmonics needed to move
    its strokes onto the requirement.  Every other harmonic is left as the
    mechanism made it, so the correction is small and the target stays near the
    reachable set.

    Args:
        design: A design whose motion seeds the target.
        expansion_stroke: Required ``STE`` [mm].
        compression_ratio: Required ``epsilon``.
        samples: Crank angles to represent the target on.
        spec: Fixed engine data.

    Returns:
        The target.

    Raises:
        ValueError: If the seed design is unanalysable, or no correction in
            these two harmonics realises the pair.
    """
    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        msg = f"cannot seed a target from an unanalysable design: {analysis.metrics.reason}"
        raise ValueError(msg)
    seed = np.asarray(analysis.require_solved().kinematics.lam, dtype=float)

    stroke = (float(compression_ratio) - 1.0) * spec.dead_volume / spec.piston_area
    theta = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    second, first = np.cos(2.0 * theta), np.sin(theta)

    def corrected(params: FloatArray) -> FloatArray:
        return seed + params[0] * second + params[1] * first

    def residual(params: FloatArray) -> list[float]:
        phases = find_phases(corrected(params))
        return [
            phases.expansion_stroke - float(expansion_stroke),
            phases.compression_stroke - stroke,
        ]

    solution, _info, status, message = fsolve(residual, np.zeros(2), full_output=True)
    if status != 1:
        msg = (
            f"no two-harmonic correction of this motion realises "
            f"STE={expansion_stroke}, eps={compression_ratio}: {message}"
        )
        raise ValueError(msg)
    return TargetMotion(
        lam=corrected(np.asarray(solution, dtype=float)),
        expansion_stroke=float(expansion_stroke),
        compression_ratio=float(compression_ratio),
        amplitude=float(solution[0]),
        skew=float(solution[1]),
    )


def describe_target(target: TargetMotion, spec: EngineSpec = DEFAULT_SPEC) -> dict[str, float]:
    """Measure a target with the same code the constraints use.

    The point of the exercise is that the target satisfies the equalities
    *as the model measures them*, so it is measured by :func:`find_phases` and
    not by the formula it was constructed from.

    Args:
        target: The target motion.
        spec: Fixed engine data.

    Returns:
        ``expansion_stroke``, ``compression_ratio`` and ``tdc_gap``.
    """
    phases = find_phases(target.lam)
    ratio = (spec.dead_volume + spec.piston_area * phases.compression_stroke) / spec.dead_volume
    return {
        "expansion_stroke": float(phases.expansion_stroke),
        "compression_ratio": float(ratio),
        "tdc_gap": float(phases.tdc_gap),
    }


def _residual(design: Design, target: TargetMotion, spec: EngineSpec) -> FloatArray | None:
    """Mean-removed difference between a mechanism's motion and the target.

    The mean is removed because the absolute height of the piston above the
    crank axis is set by the cylinder position, not by the motion the cycle
    cares about; both requirements are differences, and so is this.
    """
    analysis = analyse(design, samples=target.samples, spec=spec)
    if not analysis.valid:
        return None
    lam = np.asarray(analysis.require_solved().kinematics.lam, dtype=float)
    if lam.size != target.lam.size:
        return None
    return (lam - lam.mean()) - (target.lam - target.lam.mean())


@dataclass(frozen=True)
class FitResult:
    """One prescribed-motion fit, and how close it landed."""

    design: Design
    rms: float
    """RMS of the motion residual [mm]."""

    expansion_stroke: float
    compression_ratio: float
    stroke_error: float
    """``|STE - target|`` [mm]."""

    ratio_error: float
    """``|epsilon - target|``."""

    converged: bool
    evaluations: int


def fit_to_target(
    target: TargetMotion,
    start: Design,
    bounds: Bounds = GLOBAL_BOUNDS,
    max_evaluations: int = 300,
    spec: EngineSpec = DEFAULT_SPEC,
) -> FitResult | None:
    """Fit the linkage to a prescribed motion by least squares.

    Uses Levenberg-Marquardt through a trust region (SciPy's ``least_squares``),
    which is the natural method once the objective is a sum of squares: the
    Gauss-Newton Hessian approximation comes free from the residual Jacobian.

    Args:
        target: The motion to fit.
        start: Initial design.
        bounds: Box to search in.
        max_evaluations: Residual evaluation budget.
        spec: Fixed engine data.

    Returns:
        The fit, or ``None`` if the start is unanalysable.
    """
    if _residual(start, target, spec) is None:
        return None
    calls = 0
    penalty = np.full(target.samples, 1.0e3)

    def residual(vector: FloatArray) -> FloatArray:
        nonlocal calls
        calls += 1
        value = _residual(Design.from_array(vector), target, spec)
        # An unanalysable design is pushed away rather than raising: the
        # trust region shrinks and tries again.
        return penalty if value is None else value

    outcome = least_squares(
        residual,
        np.clip(start.to_array(), bounds.lower, bounds.upper),
        bounds=(bounds.lower, bounds.upper),
        max_nfev=max_evaluations,
        xtol=1.0e-10,
        ftol=1.0e-10,
    )
    design = Design.from_array(np.asarray(outcome.x, dtype=float))
    final = _residual(design, target, spec)
    if final is None:
        return None
    try:
        measured = describe_target(
            TargetMotion(
                lam=np.asarray(
                    analyse(design, samples=target.samples, spec=spec)
                    .require_solved()
                    .kinematics.lam,
                    dtype=float,
                ),
                expansion_stroke=target.expansion_stroke,
                compression_ratio=target.compression_ratio,
                amplitude=0.0,
                skew=0.0,
            ),
            spec=spec,
        )
    except PhaseError:
        return None
    return FitResult(
        design=design,
        rms=float(np.sqrt(np.mean(final**2))),
        expansion_stroke=measured["expansion_stroke"],
        compression_ratio=measured["compression_ratio"],
        stroke_error=abs(measured["expansion_stroke"] - target.expansion_stroke),
        ratio_error=abs(measured["compression_ratio"] - target.compression_ratio),
        converged=bool(outcome.success),
        evaluations=calls,
    )


def fit_within_constraints(
    target: TargetMotion,
    start: Design,
    bounds: Bounds = GLOBAL_BOUNDS,
    max_iterations: int = 120,
    hold_bands: bool = False,
    band: float = 0.05,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
) -> FitResult | None:
    """Fit to a prescribed motion **subject to** the inequality constraints.

    Why this is the right form of the reformulated problem
    ------------------------------------------------------
    :func:`fit_to_target` drops the constraints entirely, which is defensible
    only because the equalities were absorbed into the target.  The
    *inequalities* were never the difficulty -- they define a full-dimensional
    set -- so discarding them buys nothing and can hand back a design that
    tracks the motion beautifully and is not buildable.

    Keeping them costs almost nothing.  The reformulation has already removed
    the measure-zero part of the feasible set, so what remains is an ordinary
    box-and-inequality NLP that SQP handles directly:

    .. math::

        \\min_X \\; \\lVert \\lambda(X) - \\lambda^\\star \\rVert^2
        \\quad \\text{s.t.} \\quad g(X) \\le 0, \\; X \\in [X_{lb}, X_{ub}]

    That is the difference between a generator of points *near the manifold*
    and a generator of points that are actually usable as starts.

    ``hold_bands`` carries the same argument one step further.  Against a target
    the mechanism cannot reach, holding ``g <= 0`` and chasing the motion trade
    against each other, and the strokes drift out of the tolerance bands of
    §3.4.  Those bands are constraints too, so they belong in the problem rather
    than being checked afterwards: with them added the solve either returns a
    design that is feasible *and* in band, or reports that it found none, which
    is the useful answer in both cases.

    The same reasoning does not stop there.  A fit that satisfies the geometric
    inequalities and the bands can still fail the coupled and vehicle
    constraints, which are not included here because each evaluation would then
    carry an MDA.  Where that matters, they are the next constraints to add.

    Args:
        target: The motion to fit.
        start: Initial design.
        bounds: Box to search in.
        max_iterations: SLSQP iteration budget.
        hold_bands: Also require ``|STE - target| <= band`` and
            ``|epsilon - target| <= band`` rather than checking them after.
        band: Half-width for those two, when ``hold_bands`` is set.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        The fit, or ``None`` if the start is unanalysable or SLSQP fails to
        return a design that is analysable at the end.
    """
    from scipy.optimize import minimize

    if _residual(start, target, spec) is None:
        return None

    def objective(vector: FloatArray) -> float:
        value = _residual(Design.from_array(vector), target, spec)
        return 1.0e6 if value is None else float(np.sum(value**2))

    width = len(INEQUALITY_NAMES) + (4 if hold_bands else 0)

    def constraint(vector: FloatArray) -> FloatArray:
        analysis = analyse(Design.from_array(vector), samples=target.samples, spec=spec)
        if not analysis.valid:
            # Unanalysable reads as deeply infeasible, which steers SLSQP back
            # rather than letting it wander off the analysable set.
            return np.full(width, -1.0e3)
        # SciPy wants ``>= 0`` where the package states ``<= 0``.
        rows = -inequality_constraints(analysis, targets)
        if not hold_bands:
            return rows
        metrics = analysis.metrics
        stroke = float(metrics.expansion_stroke) - target.expansion_stroke
        ratio = float(metrics.compression_ratio) - target.compression_ratio
        return np.concatenate(
            [rows, [band - stroke, band + stroke, band - ratio, band + ratio]]
        )

    outcome = minimize(
        objective,
        np.clip(start.to_array(), bounds.lower, bounds.upper),
        method="SLSQP",
        bounds=list(zip(bounds.lower, bounds.upper, strict=True)),
        constraints=[{"type": "ineq", "fun": constraint}],
        options={"maxiter": int(max_iterations), "ftol": 1.0e-12},
    )
    design = Design.from_array(np.asarray(outcome.x, dtype=float))
    final = _residual(design, target, spec)
    if final is None:
        return None
    analysis = analyse(design, samples=target.samples, spec=spec)
    if not analysis.valid:
        return None
    metrics = analysis.metrics
    return FitResult(
        design=design,
        rms=float(np.sqrt(np.mean(final**2))),
        expansion_stroke=float(metrics.expansion_stroke),
        compression_ratio=float(metrics.compression_ratio),
        stroke_error=abs(float(metrics.expansion_stroke) - target.expansion_stroke),
        ratio_error=abs(float(metrics.compression_ratio) - target.compression_ratio),
        converged=bool(outcome.success),
        evaluations=int(outcome.nfev),
    )


CYCLE_PENALTY = 1.0e6
"""Objective value for a design whose motion is not an Atkinson cycle at all.

Worse than any range-unavailable value, so the search always climbs back to a
mechanism that at least completes four strokes before it starts caring about
how far the car goes.
"""

RANGE_UNAVAILABLE = 1.0e3
"""Objective floor for a design that is analysable but produces no range.

Any computable range scores negative, so this is worse than every real design
and better than a broken cycle -- the middle rung of the ladder.
"""

MOTION_WEIGHT = 1.0e2
"""Weight on the motion residual when it stands in for the range [1/mm^2]."""


@dataclass(frozen=True)
class RangeFit:
    """One range-maximising solve guided by a target motion."""

    design: Design
    km_per_litre: float
    rms: float
    """Motion residual against the target [mm]."""

    stroke_error: float
    ratio_error: float
    worst_constraint: float
    """Largest constraint value; ``<= 0`` means every constraint holds."""

    feasible: bool
    """Feasible against the whole model, by :func:`exlink.performance.evaluate`."""

    fell_back: int
    """Evaluations that scored on the target because the range was unavailable."""

    broken_cycle: int
    """Evaluations whose motion was not a four-stroke cycle."""

    evaluations: int
    converged: bool


def _range_constraints(
    performance: Performance,
    target: TargetMotion,
    band: float,
    targets: DesignTargets,
) -> FloatArray:
    """Every constraint of the full problem, as ``>= 0`` when satisfied.

    Geometric, tolerance bands, coupled and vehicle, in that order -- the same
    twelve §3.10 states, minus nothing.
    """
    from .disciplines import MEMBER_IS_SLENDER
    from .dynamics import MEMBER_NAMES
    from .gears import MAX_WIDTH_FACTOR
    from .sizing import MAX_DIAMETER, member_lengths

    analysis = performance.analysis
    if not analysis.valid:
        return np.full(NUMBER_OF_CONSTRAINTS, -1.0e3)

    metrics = analysis.metrics
    stroke = float(metrics.expansion_stroke) - target.expansion_stroke
    ratio = float(metrics.compression_ratio) - target.compression_ratio
    rows = [
        -inequality_constraints(analysis, targets),
        np.array([band - stroke, band + stroke, band - ratio, band + ratio]),
    ]

    coupled = performance.coupled
    if coupled is None:
        rows.append(np.full(3, -1.0e3))
    else:
        # Both arrays are indexed by MEMBER_NAMES, which is the order the
        # sizing and the slenderness mask are both written against; sorting
        # the dict would silently mis-pair diameters with lengths.
        diameters = np.array([float(coupled.diameters[name]) for name in MEMBER_NAMES])
        lengths = member_lengths(performance.design)
        slender = float(np.max((diameters / lengths)[MEMBER_IS_SLENDER], initial=0.0))
        rows.append(
            np.array(
                [
                    0.98 * MAX_DIAMETER - float(diameters.max()),
                    0.34 - slender,
                    1.0 - coupled.peak_bearing_load / targets.max_bearing_load,
                ]
            )
        )

    friction = performance.friction
    gears = performance.budget.gears
    runs = (
        friction.brake_work / max(friction.indicated_work, 1.0)
        if friction is not None
        else -1.0
    )
    gear = MAX_WIDTH_FACTOR - gears.width_factor if gears is not None else -1.0
    rows.append(np.array([runs, gear]))
    return np.concatenate(rows)


NUMBER_OF_CONSTRAINTS = 14
"""Five geometric, four band sides, three coupled, two vehicle."""


def maximise_range_from_target(
    target: TargetMotion,
    start: Design,
    speed_rpm: float = 1000.0,
    bounds: Bounds = GLOBAL_BOUNDS,
    max_iterations: int = 120,
    band: float = 0.05,
    beta_target: float | None = None,
    grade: int = 8,
    module: float | None = None,
    teeth: int | None = None,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
) -> RangeFit | None:
    """Maximise range under every constraint, with the target as a fallback.

    The formulation
    ---------------
    This is the whole problem of §3.10 with the prescribed motion carried
    along, rather than a synthesis problem that hands its answer to a separate
    optimization:

    .. math::

        \\max_X \\; R(X) \\quad \\text{s.t.} \\quad g(X) \\le 0, \\;
        |STE - 74| \\le \\delta, \\; |\\varepsilon - 16| \\le \\delta, \\;
        X \\in [X_{lb}, X_{ub}]

    where :math:`g` now carries *all* of it -- the five geometric constraints,
    the three coupled ones and the two vehicle ones -- because each previous
    round of this exercise found that whatever was left out was what the solve
    then violated.

    Why the target does not simply disappear
    ----------------------------------------
    The range is not computable everywhere.  A design whose kinematics closes
    can still fail to size, fail to run, or fail to produce a motion that is a
    four-stroke cycle at all, and at such a point :math:`R` has no value for a
    line search to descend.  Deleting the target and penalising those points
    with a constant leaves the optimizer a flat, uninformative region to cross.

    So the objective is a ladder, and the target holds the middle rung:

    ==========================================  ==========================================
    the design                                  scores
    ==========================================  ==========================================
    range computable                            :math:`-R(X)`, the real objective
    analysable, no range                        :data:`RANGE_UNAVAILABLE` plus the
                                                motion residual -- the target takes over
    motion is not a four-stroke cycle           :data:`CYCLE_PENALTY`
    ==========================================  ==========================================

    Every rung is worse than the one above it, so the search is always pushed
    back towards designs that run; and on the middle rung it still has a
    gradient to follow, because tracking :math:`\\lambda^\\star` is a proxy for
    getting back to a cycle that works.

    Args:
        target: The motion to fall back on.
        start: Initial design.
        speed_rpm: Operating point for the range.
        bounds: Box to search in.
        max_iterations: SLSQP iteration budget.
        band: Half-width on the two stroke requirements.
        beta_target: When given, the system reliability index is constrained to
            ``>= beta_target``, making this the reliability-based form of the
            problem.  ``None`` leaves it deterministic.
        grade: ISO 286 IT grade assumed for the dimensions when ``beta_target``
            is set; ignored otherwise.
        module: Gear module to pin, or ``None`` to let the sizer choose.
        teeth: Teeth on the small gear, or ``None``.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        The solve, or ``None`` if the final design is unanalysable.
    """
    from scipy.optimize import minimize

    from .performance import evaluate
    from .robustness import failure_probability

    calls = {"n": 0, "fallback": 0, "broken": 0}
    # SLSQP asks for the objective and the constraints at the same point, one
    # after the other, and each evaluation here is a converged MDA.  Without
    # this the solve pays for every point twice.
    cache: dict[bytes, Performance] = {}

    def score(vector: FloatArray) -> Performance:
        key = np.ascontiguousarray(vector, dtype=float).tobytes()
        hit = cache.get(key)
        if hit is not None:
            return hit
        # ``targets`` is deliberately not forwarded: ``evaluate`` passes its
        # extra keywords down to the sizing solve, which does not take them,
        # and the constraint bounds are applied by ``_range_constraints``
        # against the analysis rather than inside it.
        performance = evaluate(
            Design.from_array(vector),
            speed_rpm=speed_rpm,
            samples=target.samples,
            module=module,
            teeth=teeth,
            spec=spec,
        )
        # Only the current point is worth keeping: SLSQP moves on and never
        # asks again, so an unbounded cache would only hold memory.
        cache.clear()
        cache[key] = performance
        return performance

    def objective(vector: FloatArray) -> float:
        calls["n"] += 1
        performance = score(vector)
        analysis = performance.analysis
        if not analysis.valid:
            calls["broken"] += 1
            return CYCLE_PENALTY
        usable = (
            performance.friction is not None
            and performance.coupled is not None
            and performance.km_per_litre > 0.0
            and np.isfinite(performance.km_per_litre)
        )
        if usable:
            return -float(performance.km_per_litre)
        calls["fallback"] += 1
        residual = _residual(Design.from_array(vector), target, spec)
        motion = 0.0 if residual is None else float(np.mean(residual**2))
        return RANGE_UNAVAILABLE + MOTION_WEIGHT * motion

    def constraint(vector: FloatArray) -> FloatArray:
        rows = _range_constraints(score(vector), target, band, targets)
        if beta_target is None:
            return rows
        # The reliability index as one more row, on the same >= 0 convention.
        # It is evaluated against the *relaxed* bands this solve is using, not
        # the specified ones: a probability is only meaningful against the
        # bound actually being held.
        reliability = failure_probability(
            Design.from_array(vector),
            grade=grade,
            samples=target.samples,
            targets=targets,
            spec=spec,
            band={"expansion_stroke": band, "compression_ratio": band},
        )
        index = -10.0 if reliability is None else reliability.system_beta
        return np.concatenate([rows, [index - beta_target]])

    outcome = minimize(
        objective,
        np.clip(start.to_array(), bounds.lower, bounds.upper),
        method="SLSQP",
        bounds=list(zip(bounds.lower, bounds.upper, strict=True)),
        constraints=[{"type": "ineq", "fun": constraint}],
        options={"maxiter": int(max_iterations), "ftol": 1.0e-8},
    )
    design = Design.from_array(np.asarray(outcome.x, dtype=float))
    final = score(design.to_array())
    if not final.analysis.valid:
        return None
    residual = _residual(design, target, spec)
    metrics = final.analysis.metrics
    return RangeFit(
        design=design,
        km_per_litre=float(final.km_per_litre),
        rms=0.0 if residual is None else float(np.sqrt(np.mean(residual**2))),
        stroke_error=abs(float(metrics.expansion_stroke) - target.expansion_stroke),
        ratio_error=abs(float(metrics.compression_ratio) - target.compression_ratio),
        worst_constraint=-float(np.min(_range_constraints(final, target, band, targets))),
        feasible=bool(final.feasible),
        fell_back=calls["fallback"],
        broken_cycle=calls["broken"],
        evaluations=calls["n"],
        converged=bool(outcome.success),
    )


def feasible_starts(
    target: TargetMotion,
    attempts: int = 20,
    seed: int = 0,
    bounds: Bounds = GLOBAL_BOUNDS,
    reference: Design | None = None,
    scatter: float = 0.25,
    max_evaluations: int = 300,
    spec: EngineSpec = DEFAULT_SPEC,
) -> list[FitResult]:
    """Generate starting points by fitting random starts to the target.

    This is the answer to §3.4's measurement that uniform sampling finds
    0 feasible designs in 12 000 draws.  Sampling *and then fitting* is a
    different operation: the random draw only has to be analysable, and the fit
    carries it towards a motion that already satisfies both equalities.

    Args:
        target: The motion to fit to.
        attempts: Random starts to try.
        seed: Seed for the draws.
        bounds: Box to search in.
        reference: Centre of the draws; the draws are uniform over ``bounds``
            when omitted, and log-scattered around this design when given.
        scatter: Relative scatter when ``reference`` is given.
        max_evaluations: Per-fit residual budget.
        spec: Fixed engine data.

    Returns:
        The fits that converged, best (smallest stroke error) first.
    """
    rng = np.random.default_rng(seed)
    results: list[FitResult] = []
    for _ in range(int(attempts)):
        if reference is None:
            vector = rng.uniform(bounds.lower, bounds.upper)
        else:
            vector = reference.to_array() * (
                1.0 + scatter * rng.normal(size=len(VARIABLE_NAMES))
            )
            vector = np.clip(vector, bounds.lower, bounds.upper)
        fit = fit_to_target(
            target,
            Design.from_array(vector),
            bounds=bounds,
            max_evaluations=max_evaluations,
            spec=spec,
        )
        if fit is not None:
            results.append(fit)
    results.sort(key=lambda item: item.stroke_error)
    return results


def fit_report(
    results: list[FitResult], stroke_band: float = 0.05, ratio_band: float = 0.05
) -> str:
    """Render what a batch of fits achieved against the tolerance bands.

    Args:
        results: Fits, as returned by :func:`feasible_starts`.
        stroke_band: Half-width on ``STE`` [mm].
        ratio_band: Half-width on ``epsilon``.

    Returns:
        A table, one row per fit, with the in-band count.
    """
    lines = ["prescribed-motion fits", "=" * 22, ""]
    lines.append(f"  {'rms [mm]':>10} {'|dSTE|':>10} {'|deps|':>10}  in band")
    inside = 0
    for fit in results:
        ok = fit.stroke_error <= stroke_band and fit.ratio_error <= ratio_band
        inside += int(ok)
        lines.append(
            f"  {fit.rms:>10.4f} {fit.stroke_error:>10.4f} {fit.ratio_error:>10.4f}"
            f"  {'yes' if ok else 'no'}"
        )
    lines.append("")
    lines.append(f"  {inside} of {len(results)} fits land inside both bands")
    return "\n".join(lines)
