"""Checking the first-order reliability estimate against sampling.

Why the check is needed
-----------------------
Every probability in this study comes from FORM: the constraint surfaces are
linearised at the nominal design and the failure probability is read off a
multivariate-normal orthant (§3.8).  That is what makes the estimate affordable
enough to sit inside an optimizer -- one exact Jacobian per design point -- and
it is an approximation whose error is not bounded by anything the estimate
itself reports.  §7.2 lists it as a limitation and §7.3 asks for a sampling
check.

This module is that check.  It evaluates the *exact* constraints on random
builds drawn from the same covariance and counts how many miss, so the number
it produces makes no linearity assumption at all.

Why crude Monte Carlo is not enough on its own
----------------------------------------------
The study's headline is :math:`P_f = 1.3\\times10^{-3}` and §6.2 reaches
:math:`1.9\\times10^{-5}`.  Crude sampling resolves a probability :math:`p` to a
relative standard error of :math:`\\sqrt{(1-p)/(np)}`, so ten per cent on
:math:`10^{-5}` needs :math:`10^{7}` analyses -- six hours, to check a number
that took milliseconds.  Crude sampling is therefore offered for the range
where it works, and the tail is reached by **importance sampling**.

How the importance sampling is set up
-------------------------------------
Work in the standard normal space :math:`u`, where a build is
:math:`X = X_0 + \\Sigma^{1/2} u` and :math:`\\Sigma` is diagonal.  For a
linearised constraint, the nearest failure point to the origin -- the *most
probable point* -- is

.. math:: u^\\star_i = \\beta_i \\alpha_i, \\qquad
          \\alpha_i = -\\frac{\\Sigma^{1/2}\\nabla g_i}{\\sigma_i},

and the failure region lies beyond it.  Sampling instead from
:math:`\\mathcal{N}(u^\\star, I)` puts most of the draws where failures actually
are, and unbiasing the count with the likelihood ratio

.. math:: w(u) = \\frac{\\varphi(u)}{\\varphi(u - u^\\star)}
                = \\exp\\bigl(-u^\\star\\!\\cdot u + \\tfrac12 |u^\\star|^2\\bigr)

recovers the true probability.  The shift is chosen at the constraint with the
*smallest* :math:`\\beta`, which is where the failures are; the estimator stays
unbiased whatever shift is used, so a poor choice costs variance rather than
correctness -- and the reported standard error says when that has happened.

Note what is and is not checked.  The design point comes from the *linearised*
constraint, but every draw is evaluated exactly, so a curved constraint surface
shows up in the answer.  What is assumed is the input distribution: normal
errors of the stated width on each dimension, independent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import Design
from .model import analyse
from .robustness import (
    ANGULAR_TOLERANCE,
    DEFAULT_GRADE,
    RELIABILITY_NAMES,
    SIGMA_PER_HALF_WIDTH,
    _band_widths,
    constraint_jacobian,
    constraint_moments,
    tolerance_half_widths,
)

FloatArray = NDArray[np.float64]

DEFAULT_DRAWS = 20_000
"""Draws for the importance-sampled estimate.

Enough for a few per cent on a probability anywhere from 1 to 1e-8, because the
shifted density puts a large fraction of the draws in the failure region
whatever its size.
"""


@dataclass(frozen=True)
class SampledReliability:
    """A failure probability estimated by sampling the exact constraints."""

    names: tuple[str, ...]
    per_constraint: dict[str, float]
    """Probability that each constraint alone is missed."""

    system: float
    """Probability that *any* constraint is missed."""

    standard_error: float
    """Standard error of :attr:`system`."""

    draws: int
    unbuildable: int
    """Draws whose kinematics did not close; counted as failures."""

    importance: bool
    """Whether the estimate is importance sampled or crude."""

    effective_draws: float
    """Kish effective sample size, ``(sum w)^2 / sum w^2``.

    For a crude estimate this equals :attr:`draws`.  For an importance-sampled
    one it is the number of *independent* draws the weighted sample is worth,
    and a value far below ``draws`` means the shift was badly chosen and the
    estimate is noisier than its count suggests.
    """

    @property
    def interval(self) -> tuple[float, float]:
        """Approximate 95 % interval on :attr:`system`, clipped to [0, 1]."""
        half = 1.96 * self.standard_error
        return (max(self.system - half, 0.0), min(self.system + half, 1.0))

    def agrees_with(self, form: float, factor: float = 3.0) -> bool:
        """Whether a FORM estimate sits within ``factor`` of this one.

        A factor rather than an interval, because the two are compared over
        orders of magnitude and the question asked of FORM is whether it has
        the right *size*, not the right third digit.

        Args:
            form: The first-order estimate.
            factor: Ratio counted as agreement.

        Returns:
            Whether they agree.
        """
        if self.system <= 0.0:
            return form <= 1.0 / max(self.draws, 1)
        ratio = form / self.system
        return 1.0 / factor <= ratio <= factor


def _exact_constraints(
    design: Design,
    samples: int,
    targets: DesignTargets,
    spec: EngineSpec,
    band: dict[str, float] | None,
) -> FloatArray | None:
    """The constraints of :data:`RELIABILITY_NAMES`, in that order, not linearised.

    The two strokes appear once per top dead centre, because that is how
    :mod:`exlink.robustness` now prices them and a check that scored a
    different constraint set would be checking the wrong thing.
    """
    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        return None
    metrics = analysis.metrics
    phases = analysis.require_solved().thermodynamics.phases
    band_stroke, band_ratio = (float(width) for width in _band_widths(band))

    stroke = tuple(value - targets.expansion_stroke for value in phases.expansion_strokes)
    ratio = tuple(
        1.0 + spec.piston_area * value / spec.dead_volume - targets.compression_ratio
        for value in phases.compression_strokes
    )
    return np.array(
        [
            metrics.rod_angle - targets.max_rod_angle,
            metrics.compatibility - targets.max_transmission,
            metrics.tdc_gap - targets.max_tdc_gap,
            metrics.side_load_ratio - targets.max_side_load,
            stroke[0] - band_stroke,
            stroke[1] - band_stroke,
            -stroke[0] - band_stroke,
            -stroke[1] - band_stroke,
            ratio[0] - band_ratio,
            ratio[1] - band_ratio,
            -ratio[0] - band_ratio,
            -ratio[1] - band_ratio,
        ]
    )


def sampled_reliability(
    design: Design,
    draws: int = DEFAULT_DRAWS,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    band: dict[str, float] | None = None,
    seed: int = 0,
    importance: bool = True,
) -> SampledReliability | None:
    """Estimate the system failure probability by sampling exact builds.

    Args:
        design: The nominal design.
        draws: Number of builds to draw.
        grade: ISO 286 IT grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        samples: Crank angles per revolution in each analysis.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        band: Half-widths of the two relaxed equalities.
        seed: Random seed, so the estimate is reproducible.
        importance: Shift the sampling density to the most probable failure
            point.  Leave it on unless the point of the run is to check the
            importance sampling itself against crude sampling.

    Returns:
        The estimate, or ``None`` if the nominal design cannot be analysed.
    """
    moments = constraint_moments(design, grade, angular, samples, targets, spec, band)
    if moments is None:
        return None

    sigma_x = tolerance_half_widths(design, grade, angular) / SIGMA_PER_HALF_WIDTH
    base = design.to_array()
    rng = np.random.default_rng(seed)

    shift = np.zeros(base.size)
    if importance:
        exact = constraint_jacobian(design, samples, targets, spec, band)
        if exact is not None:
            shift = _design_point(moments.beta, exact[1], sigma_x)
    unit = rng.normal(size=(int(draws), base.size)) + shift
    # phi(u) / phi(u - shift), in logs for the tail.
    log_weight = -unit @ shift + 0.5 * float(shift @ shift)
    weight = np.exp(log_weight)

    failed = np.zeros(int(draws), dtype=bool)
    collected = np.zeros((int(draws), len(RELIABILITY_NAMES)))
    breached = np.zeros((int(draws), len(RELIABILITY_NAMES)), dtype=bool)
    unbuildable = 0
    for index, row in enumerate(unit):
        values = _exact_constraints(
            Design.from_array(base + sigma_x * row), samples, targets, spec, band
        )
        if values is None:
            # A build whose kinematics does not close has missed every
            # requirement there is; counting it as a success would be the one
            # way to make this estimate optimistic.
            unbuildable += 1
            failed[index] = True
            collected[index] = 1.0
            breached[index] = True
            continue
        collected[index] = values
        breached[index] = values > 0.0
        failed[index] = bool(breached[index].any())

    total = float(weight.sum())
    system = float((weight * failed).sum() / total) if total > 0.0 else 1.0
    # Var of a weighted mean of a Bernoulli, by the delta method.
    residual = weight * (failed.astype(float) - system)
    error = float(np.sqrt((residual**2).sum())) / max(total, 1.0e-300)
    effective = total**2 / float((weight**2).sum()) if total > 0.0 else 0.0

    if importance:
        rates = [
            float((weight * breached[:, column]).sum() / total) if total > 0.0 else 1.0
            for column in range(len(RELIABILITY_NAMES))
        ]
    else:
        rates = list(np.atleast_1d(_exceedance_rate(collected)))

    return SampledReliability(
        names=RELIABILITY_NAMES,
        per_constraint=dict(zip(RELIABILITY_NAMES, map(float, rates), strict=True)),
        system=system,
        standard_error=error,
        draws=int(draws),
        unbuildable=unbuildable,
        importance=bool(importance),
        effective_draws=effective,
    )


def _exceedance_rate(values: FloatArray) -> FloatArray:
    """Fraction of unweighted draws in which each constraint is violated.

    Routed through ``gemseo-umdo``'s own ``Probability`` estimator when that
    plugin is installed -- it is the statistic a sampling-based U-MDO
    formulation would attach to each constraint, and going through it keeps the
    check on the same footing as the rest of the GEMSEO stack.  The plugin is
    the project's optional ``uq`` extra, so the one-line equivalent is used
    when it is absent; the two agree exactly, and a test pins that they do.

    Args:
        values: ``(n_draws, n_constraints)`` constraint values.

    Returns:
        One rate per constraint.
    """
    try:
        from gemseo_umdo.formulations._statistics.sampling.probability import (
            Probability,
        )
    except ImportError:
        return np.asarray((values >= 0.0).mean(axis=0), dtype=float)
    return np.asarray(
        Probability(threshold=0.0, greater=True).estimate_statistic(values), dtype=float
    )


def _design_point(beta: FloatArray, jacobian: FloatArray, sigma_x: FloatArray) -> FloatArray:
    """The most probable failure point, in standard normal coordinates.

    Taken at the constraint with the smallest reliability index, which is where
    the failures are.  The estimator is unbiased for any shift, so this choice
    controls the variance and not the answer.

    Args:
        beta: Reliability index of each constraint.
        jacobian: ``(n_constraints, 11)`` exact gradients.
        sigma_x: Standard deviation of each design variable.

    Returns:
        ``beta alpha`` for the binding constraint, ``(11,)``.
    """
    finite = np.isfinite(beta)
    if not finite.any():
        return np.zeros(sigma_x.size)
    order = np.where(finite, beta, np.inf)
    binding = int(np.argmin(order))

    scaled = np.asarray(jacobian[binding], dtype=float) * sigma_x
    norm = float(np.linalg.norm(scaled))
    if norm <= 0.0:
        return np.zeros(sigma_x.size)
    # alpha points *towards* failure, which is the direction g increases.
    return float(order[binding]) * scaled / norm


def format_sampled(sampled: SampledReliability, form: float | None = None) -> str:
    """Render a :class:`SampledReliability` as an aligned table."""
    lines = ["sampled reliability", "=" * 52, ""]
    lines.append(f"  {'constraint':<15}{'P(fail), sampled':>18}")
    for name in sampled.names:
        lines.append(f"  {name:<15}{sampled.per_constraint[name]:>18.4e}")
    low, high = sampled.interval
    lines.append("")
    kind = "importance sampled" if sampled.importance else "crude Monte Carlo"
    lines.append(f"  {kind}, {sampled.draws} draws")
    lines.append(f"  effective sample size   {sampled.effective_draws:.0f}")
    lines.append(f"  builds that did not close {sampled.unbuildable}")
    lines.append(f"  system P(fail)          {sampled.system:.4e}")
    lines.append(f"  95 % interval           [{low:.4e}, {high:.4e}]")
    if form is not None:
        lines.append(f"  FORM says               {form:.4e}")
        ratio = form / sampled.system if sampled.system > 0.0 else float("inf")
        lines.append(f"  FORM / sampled          {ratio:.3f}")
    return "\n".join(lines)
