"""Prescribed-motion synthesis: fit the linkage to a target piston motion.

The formulation of :mod:`exlink.scenarios` asks for the best mechanism and lets
the piston motion fall out of it.  Kinematic synthesis conventionally does the
opposite -- it *prescribes* a target motion and fits the linkage to it:

.. math::
    \\min_X\\; J(X) = \\sum_k \\bigl(\\lambda_k(X) - \\lambda^\\star_k\\bigr)^2

This module implements that, and exists to answer one question the main
formulation cannot: **where do feasible starting points come from?**

Why it matters here
-------------------
The two equality requirements are functionals of :math:`\\lambda` alone --
:math:`STE` is the span from top dead centre to the deeper bottom dead centre,
and :math:`\\varepsilon` follows from the shallower one through the swept
volume.  So a target built to have the right two strokes satisfies both
equalities *as the model measures them*, before any linkage is involved:
:func:`target_motion` solves for one and :func:`describe_target` confirms it to
machine precision.

That makes the fit a sampling device.  Uniform sampling of the design box finds
a design on the equality manifold with probability zero -- measured at 0 in
12 000 attempts -- because the manifold has measure zero.  Fitting to a target
that is *already* on the manifold is instead an unconstrained box problem, so it
can be run from arbitrary random starts, and each converged fit lands near the
manifold rather than nowhere near it.  :func:`feasible_starts` does exactly
that.

What this does *not* claim
--------------------------
The fit is over-determined -- eleven variables against several hundred
residuals -- so :math:`\\lambda(X) = \\lambda^\\star` is not attainable and the
equalities are *not* satisfied exactly.  Whether the residual is small enough to
land inside the tolerance bands of §3.4 is an empirical question, and
:func:`fit_report` is what measures it.  A least-squares residual also has no
exchange rate with mass, friction or range: this is a generator of feasible
points, not a substitute for the objective.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import fsolve, least_squares

from .constants import DEFAULT_SPEC, EngineSpec
from .cycle import PhaseError, find_phases
from .design import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design
from .materials import FloatArray
from .model import analyse

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
