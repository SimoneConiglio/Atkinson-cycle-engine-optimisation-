"""What manufacturing tolerance does to a design that sits on a knife edge.

The central finding of the coupled study is that the quasi-statically optimal
mechanism sits at the transmission-angle singularity.  That is a statement
about *conditioning*, and a design chosen for its nominal performance in a
badly conditioned region is exactly the design a tolerance study exists to
catch.  Presenting a deterministic optimum here without one would be
negligent: the whole point of the finding is that small perturbations matter
enormously.

Two propagation methods, and why both
--------------------------------------
**First order, from the exact Jacobians.**  For a constraint ``g(X)`` and
independent dimensional errors with covariance ``Sigma``,

.. math:: \\sigma_g^2 = \\nabla g^T \\Sigma \\nabla g

Because :mod:`exlink.jacobian` already supplies ``\\nabla g`` analytically, a
full robustness assessment costs one extra Jacobian evaluation.  The exact
derivatives built to make the optimizer work turn out to pay for the
uncertainty analysis as well -- the same object serves both.

**Monte Carlo, to check the first order is not lying.**  Linearisation is
exactly what should be distrusted near a singularity, where the constraint
surfaces curve sharply.  A sample confirms or refutes the linear estimate, and
the comparison between them is itself a measure of how nonlinear the design
point is.

Where the tolerances come from
-------------------------------
ISO 286 IT grades, not invented numbers.  The standard tolerance unit is

.. math:: i = 0.45 \\sqrt[3]{D} + 0.001 D \\quad [\\mu m]

for a nominal size ``D`` in millimetres, and grade ``IT_n`` is a fixed multiple
of it (16i at IT7, 25i at IT8, 40i at IT9).  A machined linkage member is an
IT8 part; holding IT7 costs real money and IT9 is loose for a bearing centre.

Angular variables are assembly rather than machining quantities -- ``theta_f``
and ``theta_r`` are set by how the gears are clocked -- so they get a fixed
angular band instead.

What the answer turns out to be
--------------------------------
The top-dead-centre gap constraint is ``g <= 0.01 mm``.  The dimensions that
produce ``g`` are held to about ``+/- 0.02 mm`` each at IT8.  The constraint is
therefore **tighter than the tolerance of the parts that produce it**, and no
amount of optimization fixes that: it is a specification defect, not a design
one.  :func:`tolerance_report` quantifies it, and it is the strongest argument
in this package for treating ``g`` as a quantity to be *minimised* -- or
adjusted at assembly with a shim -- rather than bounded.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from gemseo.core.discipline import Discipline
from gemseo.typing import StrKeyMapping

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import ANGULAR_VARIABLES, VARIABLE_NAMES, Design
from .jacobian import kinematic_jacobian, metric_jacobian
from .materials import FloatArray
from .model import analyse, equality_constraints, inequality_constraints

IT_FACTORS: dict[int, float] = {6: 10.0, 7: 16.0, 8: 25.0, 9: 40.0, 10: 64.0, 11: 100.0}
"""Multiples of the standard tolerance unit ``i`` for each ISO 286 IT grade."""

DEFAULT_GRADE = 8
"""IT grade a machined linkage member is normally held to."""

ANGULAR_TOLERANCE = 0.05
"""Half-width of the assembly tolerance on the two clocking angles [deg]."""

EQUALITY_BAND: dict[str, float] = {
    "expansion_stroke": 0.05,
    "compression_ratio": 0.05,
}
"""Design band each equality is really held to.

``STE = 74`` and ``epsilon = 16`` are written as equalities, but no part is
made to an exact dimension, so a perturbation analysis that scores ``|residual|
> 0`` as a violation would report 100 % failure for any design whatsoever and
say nothing.  The equalities are scored against the same half-widths the
optimizer relaxes them to, which is what a drawing would actually call out.
"""

SIGMA_PER_HALF_WIDTH = 3.0
"""Standard deviations spanned by a tolerance half-width.

A ``+/- t`` drawing callout is read as ``3 sigma``, i.e. a process capability of
``Cp = 1``.  Reading it as ``1 sigma`` would be pessimistic and as ``6 sigma``
optimistic; this is the conventional middle.
"""

#: Names of the constraints, in the order the vectors below use.
CONSTRAINT_NAMES: tuple[str, ...] = (
    "expansion_stroke",
    "compression_ratio",
    "rod_angle",
    "compatibility",
    "tdc_gap",
    "clearance",
    "side_load",
)


def tolerance_unit(size: float) -> float:
    """ISO 286 standard tolerance unit ``i`` [mm].

    Args:
        size: Nominal size ``D`` [mm].

    Returns:
        ``i = (0.45 D^(1/3) + 0.001 D)`` in millimetres.
    """
    d = max(abs(float(size)), 1.0)
    return (0.45 * d ** (1.0 / 3.0) + 0.001 * d) * 1.0e-3


def tolerance_half_widths(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
) -> FloatArray:
    """Tolerance half-width on each design variable.

    Args:
        design: The nominal design.
        grade: ISO 286 IT grade for the dimensional variables.
        angular: Half-width on the two clocking angles [deg].

    Returns:
        One half-width per entry of :data:`~exlink.design.VARIABLE_NAMES`,
        in that variable's own unit (mm, or degrees for the angles).
    """
    factor = IT_FACTORS[grade]
    widths = []
    for name in VARIABLE_NAMES:
        if name in ANGULAR_VARIABLES:
            widths.append(float(angular))
        else:
            widths.append(0.5 * factor * tolerance_unit(getattr(design, name)))
    return np.array(widths)


def covariance(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
) -> FloatArray:
    """Diagonal covariance of the manufacturing errors.

    Errors on separate features are taken independent, which is the right
    default for parts made on separate operations and is conservative for the
    two angles, whose errors would partly cancel if they were clocked together.

    Args:
        design: The nominal design.
        grade: ISO 286 IT grade.
        angular: Angular half-width [deg].

    Returns:
        An ``(11, 11)`` diagonal covariance matrix.
    """
    sigma = tolerance_half_widths(design, grade, angular) / SIGMA_PER_HALF_WIDTH
    return np.diag(sigma**2)


def _band_widths(band: dict[str, float] | None) -> FloatArray:
    """Half-widths of the two relaxed equalities, defaulting to the specified ones."""
    widths = dict(EQUALITY_BAND)
    if band:
        widths.update(band)
    return np.array([widths["expansion_stroke"], widths["compression_ratio"]])


def _constraint_vector(
    design: Design,
    samples: int,
    targets: DesignTargets,
    spec: EngineSpec,
    band: dict[str, float] | None = None,
) -> FloatArray | None:
    """All seven constraints as one vector, negative meaning satisfied.

    ``band`` widens the two relaxed equalities.  It has to be a parameter
    rather than a constant, because §6.2's question -- what bound would this
    design need in order to be reliable? -- cannot be asked of a formulation
    whose bounds are fixed.
    """
    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        return None
    equality = equality_constraints(analysis, targets)
    inequality = inequality_constraints(analysis, targets)
    banded = np.abs(equality) - _band_widths(band)
    return np.concatenate([banded, inequality])


@dataclass(frozen=True)
class ToleranceReport:
    """How a design stands up to the tolerances of the parts that make it."""

    design: Design
    grade: int
    half_widths: FloatArray
    """Tolerance half-width on each design variable."""

    nominal: dict[str, float]
    """Constraint value at the nominal design; negative is satisfied."""

    linear_sigma: dict[str, float]
    """First-order standard deviation of each constraint."""

    monte_carlo_sigma: dict[str, float]
    """Sampled standard deviation of each constraint."""

    violation_rate: dict[str, float]
    """Fraction of samples violating each constraint."""

    any_violation_rate: float
    """Fraction of samples violating at least one constraint."""

    unbuildable_rate: float
    """Fraction of samples the kinematics could not close at all."""

    samples: int

    def capability(self) -> dict[str, float]:
        """One-sided process capability ``Cpk`` of each constraint.

        ``Cpk = -g_nominal / (3 sigma)``: how many three-sigma widths of margin
        the nominal design holds.  Below 1 the constraint is violated by
        ordinary manufacturing variation; below 0 the nominal design is already
        outside.  A value of 1.33 is the usual industrial target.
        """
        result = {}
        for name in CONSTRAINT_NAMES:
            sigma = self.monte_carlo_sigma.get(name, 0.0)
            if sigma <= 0.0:
                result[name] = float("inf")
                continue
            result[name] = -self.nominal[name] / (3.0 * sigma)
        return result

    def linearity_error(self) -> dict[str, float]:
        """Ratio of first-order to sampled sigma, per constraint.

        Far from 1 means the design point is strongly nonlinear over its own
        tolerance band -- which is precisely the symptom of sitting near a
        singularity, and a warning that a first-order robust formulation would
        be optimistic there.
        """
        result = {}
        for name in CONSTRAINT_NAMES:
            sampled = self.monte_carlo_sigma.get(name, 0.0)
            if sampled <= 0.0:
                result[name] = float("nan")
                continue
            result[name] = self.linear_sigma.get(name, 0.0) / sampled
        return result


def tolerance_report(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 2000,
    crank_samples: int = 720,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    seed: int = 0,
) -> ToleranceReport:
    """Propagate manufacturing tolerance through the constraints, two ways.

    Args:
        design: The nominal design.
        grade: ISO 286 IT grade for the dimensional variables.
        angular: Angular half-width [deg].
        samples: Monte Carlo sample count.
        crank_samples: Crank angles per revolution in each analysis.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        seed: Random seed, so the report is reproducible.

    Returns:
        The report.

    Raises:
        ValueError: If the nominal design cannot be analysed.
    """
    nominal_vector = _constraint_vector(design, crank_samples, targets, spec)
    if nominal_vector is None:
        msg = "cannot assess the tolerance of a design that does not close"
        raise ValueError(msg)

    sigma_matrix = covariance(design, grade, angular)
    widths = tolerance_half_widths(design, grade, angular)

    # -- first order, from the exact Jacobian --------------------------------------
    rows = _constraint_rows(design, crank_samples, targets, spec)
    linear = {
        name: float(np.sqrt(row @ sigma_matrix @ row))
        for name, row in zip(CONSTRAINT_NAMES, rows, strict=True)
    }

    # -- Monte Carlo ---------------------------------------------------------------
    rng = np.random.default_rng(seed)
    base = design.to_array()
    deviations = rng.normal(0.0, widths / SIGMA_PER_HALF_WIDTH, size=(samples, base.size))
    collected: list[FloatArray] = []
    unbuildable = 0
    for row in deviations:
        vector = _constraint_vector(Design.from_array(base + row), crank_samples, targets, spec)
        if vector is None:
            unbuildable += 1
            continue
        collected.append(vector)

    if collected:
        stacked = np.stack(collected)
        sampled = {
            name: float(np.std(stacked[:, index]))
            for index, name in enumerate(CONSTRAINT_NAMES)
        }
        violated = stacked > 0.0
        rate = {
            name: float(np.mean(violated[:, index]))
            for index, name in enumerate(CONSTRAINT_NAMES)
        }
        any_rate = float(np.mean(np.any(violated, axis=1)))
    else:
        sampled = dict.fromkeys(CONSTRAINT_NAMES, 0.0)
        rate = dict.fromkeys(CONSTRAINT_NAMES, 1.0)
        any_rate = 1.0

    built = max(samples - unbuildable, 1)
    scale = built / samples
    return ToleranceReport(
        design=design,
        grade=grade,
        half_widths=widths,
        nominal=dict(zip(CONSTRAINT_NAMES, nominal_vector.tolist(), strict=True)),
        linear_sigma=linear,
        monte_carlo_sigma=sampled,
        violation_rate={name: value * scale for name, value in rate.items()},
        any_violation_rate=any_rate * scale + (1.0 - scale),
        unbuildable_rate=unbuildable / samples,
        samples=samples,
    )


def _constraint_rows(
    design: Design,
    crank_samples: int,
    targets: DesignTargets,
    spec: EngineSpec,
) -> list[FloatArray]:
    """Gradient rows of the seven constraints, in :data:`CONSTRAINT_NAMES` order.

    Six come straight from :func:`~exlink.jacobian.metric_jacobian`, exactly.
    The seventh, the cylinder clearance, has no closed-form derivative in this
    package -- it is a minimum over both the crank revolution and the three
    edges of the trigonal link, and the edge attaining it can change -- so it
    falls back to a central difference.  That is honest about which parts are
    exact and which are not, and the clearance is far from active anyway.

    The two equalities enter as ``|residual|``, whose gradient is the residual
    gradient up to a sign; the sign does not affect ``sigma``.
    """
    analysis = analyse(design, samples=crank_samples, spec=spec)
    kinematic = kinematic_jacobian(design, analysis.require_solved().kinematics, spec)
    exact = metric_jacobian(design, analysis, kinematic, spec)

    base = design.to_array()
    clearance = np.zeros(base.size)
    for index in range(base.size):
        step = max(abs(base[index]), 1.0) * 1.0e-6
        rows = []
        for sign in (1.0, -1.0):
            shifted = base.copy()
            shifted[index] += sign * step
            probe = analyse(Design.from_array(shifted), samples=crank_samples, spec=spec)
            rows.append(probe.metrics.clearance if probe.valid else float("nan"))
        clearance[index] = (rows[0] - rows[1]) / (2.0 * step)

    return [
        exact["stroke_error"],
        exact["compression_ratio_error"],
        exact["rod_angle_margin"],
        exact["compatibility_margin"],
        exact["tdc_gap_margin"],
        -clearance,
        exact["side_load_margin"],
    ]


def format_report(report: ToleranceReport) -> str:
    """Render a :class:`ToleranceReport` as an aligned table."""
    title = f"tolerance study at IT{report.grade}, {report.samples} samples"
    lines = [title, "=" * len(title), ""]
    lines.append(
        f"  {'constraint':<20}{'nominal':>12}{'sigma_1st':>12}"
        f"{'sigma_MC':>12}{'Cpk':>8}{'violated':>11}"
    )
    capability = report.capability()
    for name in CONSTRAINT_NAMES:
        lines.append(
            f"  {name:<20}{report.nominal[name]:>12.4g}"
            f"{report.linear_sigma[name]:>12.4g}"
            f"{report.monte_carlo_sigma[name]:>12.4g}"
            f"{capability[name]:>8.2f}"
            f"{100.0 * report.violation_rate[name]:>10.1f}%"
        )
    lines.append("")
    lines.append(f"  any constraint violated: {100.0 * report.any_violation_rate:.1f}%")
    lines.append(f"  does not close at all:   {100.0 * report.unbuildable_rate:.1f}%")
    return "\n".join(lines)


def required_grade(
    design: Design,
    constraint: str = "tdc_gap",
    target_capability: float = 1.33,
    crank_samples: int = 720,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
) -> tuple[int | None, float]:
    """The loosest IT grade at which a constraint is actually holdable.

    Constraint sigma is linear in the tolerance factor, so the whole grade
    ladder can be scanned from one Jacobian.  Answering "which grade would fix
    this?" is what turns a robustness finding into an engineering decision: if
    the answer is a grade a machine shop can hold, tighten the drawing; if it
    is off the bottom of the ladder, the constraint itself is the problem.

    Args:
        design: The nominal design.
        constraint: Which of :data:`CONSTRAINT_NAMES` to assess.
        target_capability: ``Cpk`` to reach; 1.33 is the usual target.
        crank_samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        ``(grade, factor_needed)``.  The grade is the loosest listed one that
        reaches the target, or ``None`` when even the tightest does not -- in
        which case ``factor_needed`` says how far below the ladder the
        requirement sits, and the constraint has to be met by adjustment at
        assembly rather than by machining.

    Raises:
        ValueError: If the nominal design does not satisfy the constraint, or
            cannot be analysed.
    """
    index = CONSTRAINT_NAMES.index(constraint)
    nominal = _constraint_vector(design, crank_samples, targets, spec)
    if nominal is None:
        msg = "cannot assess the tolerance of a design that does not close"
        raise ValueError(msg)
    margin = -float(nominal[index])
    if margin <= 0.0:
        msg = f"the nominal design already violates {constraint}"
        raise ValueError(msg)

    row = _constraint_rows(design, crank_samples, targets, spec)[index]
    reference = DEFAULT_GRADE
    sigma_matrix = covariance(design, reference)
    sigma_reference = float(np.sqrt(row @ sigma_matrix @ row))
    if sigma_reference <= 0.0:
        return min(IT_FACTORS), 0.0

    # sigma scales linearly with the IT factor.
    allowed = margin / (3.0 * target_capability)
    factor_needed = IT_FACTORS[reference] * allowed / sigma_reference
    workable = [g for g, f in IT_FACTORS.items() if f <= factor_needed]
    return (max(workable) if workable else None), factor_needed


#: The constraints whose reliability is assessed, and the Jacobian row for each.
#: ``clearance`` is left out deliberately: its capability is above 400, so its
#: failure probability is numerically zero, and it is the one constraint here
#: without an analytic gradient.
CONSTRAINT_ROWS_FOR_RELIABILITY: dict[str, str] = {
    "rod_angle": "rod_angle_margin",
    "compatibility": "compatibility_margin",
    "tdc_gap": "tdc_gap_margin",
    "side_load": "side_load_margin",
    "stroke_upper": "stroke_error",
    "stroke_lower": "stroke_error",
    "ratio_upper": "compression_ratio_error",
    "ratio_lower": "compression_ratio_error",
}

RELIABILITY_NAMES: tuple[str, ...] = tuple(CONSTRAINT_ROWS_FOR_RELIABILITY)

TARGET_FAILURE_PROBABILITY = 1.0e-3
"""Default system probability of failure to design to.

A design target, not a fact about the parts: it says how often a
nominally-conforming build may miss *any* of its requirements.  1e-3 corresponds
to a reliability index of about 3.1 on a single constraint, which is the usual
structural-design level.
"""


@dataclass(frozen=True)
class ConstraintMoments:
    """First-order moments of every constraint under manufacturing scatter."""

    names: tuple[str, ...]
    value: FloatArray
    """Nominal constraint value; negative is satisfied."""

    sigma: FloatArray
    """Standard deviation induced by the tolerances."""

    correlation: FloatArray
    """Correlation between constraints, ``(n, n)``.

    Non-trivial and the whole reason a system probability differs from the
    product of the individual ones: every constraint is a function of the same
    eleven dimensions, so their errors are strongly dependent.  The stroke and
    compression-ratio residuals in particular move together.
    """

    @property
    def beta(self) -> FloatArray:
        """Reliability index ``-g / sigma`` of each constraint.

        Large and positive is safe.  This is the first-order (FORM) index: the
        distance, in standard deviations, from the nominal design to the
        constraint surface along the direction the gradient says is steepest.
        """
        safe = np.where(self.sigma > 0.0, self.sigma, np.inf)
        return -self.value / safe


def constraint_moments(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    band: dict[str, float] | None = None,
) -> ConstraintMoments | None:
    """Constraint values, standard deviations and correlations, from the exact Jacobians.

    First-order propagation: with ``Sigma`` the covariance of the dimensional
    errors,

    .. math::
        \\sigma_i = \\sqrt{\\nabla g_i^\\top \\Sigma \\nabla g_i}, \\qquad
        \\rho_{ij} = \\frac{\\nabla g_i^\\top \\Sigma \\nabla g_j}
                          {\\sigma_i \\sigma_j}

    The correlation matters.  Every constraint here is a function of the same
    eleven dimensions, so their scatter is strongly dependent, and treating them
    as independent -- which is what applying a fixed margin to each one
    separately amounts to -- overstates the chance that *some* constraint fails.

    Args:
        design: The design to assess.
        grade: ISO 286 IT grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        The moments, or ``None`` if the design cannot be analysed.
    """
    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        return None

    kinematic = kinematic_jacobian(design, analysis.require_solved().kinematics, spec)
    rows = metric_jacobian(design, analysis, kinematic, spec)
    sigma_matrix = covariance(design, grade, angular)

    metrics = analysis.metrics
    stroke = metrics.expansion_stroke - targets.expansion_stroke
    ratio = metrics.compression_ratio - targets.compression_ratio
    band_stroke, band_ratio = (float(width) for width in _band_widths(band))

    values = {
        "rod_angle": metrics.rod_angle - targets.max_rod_angle,
        "compatibility": metrics.compatibility - targets.max_transmission,
        "tdc_gap": metrics.tdc_gap - targets.max_tdc_gap,
        "side_load": metrics.side_load_ratio - targets.max_side_load,
        "stroke_upper": stroke - band_stroke,
        "stroke_lower": -stroke - band_stroke,
        "ratio_upper": ratio - band_ratio,
        "ratio_lower": -ratio - band_ratio,
    }
    gradients = []
    for name, source in CONSTRAINT_ROWS_FOR_RELIABILITY.items():
        gradient = np.asarray(rows[source], dtype=float)
        gradients.append(-gradient if name.endswith("_lower") else gradient)
    jacobian = np.stack(gradients)

    covariances = jacobian @ sigma_matrix @ jacobian.T
    sigma = np.sqrt(np.clip(np.diag(covariances), 0.0, None))
    outer = np.outer(sigma, sigma)
    correlation = np.divide(covariances, outer, out=np.eye(sigma.size), where=outer > 0.0)
    return ConstraintMoments(
        names=RELIABILITY_NAMES,
        value=np.array([values[name] for name in RELIABILITY_NAMES]),
        sigma=sigma,
        correlation=np.clip(correlation, -1.0, 1.0),
    )


@dataclass(frozen=True)
class Reliability:
    """Probability that a nominally-conforming build misses its requirements."""

    moments: ConstraintMoments
    per_constraint: dict[str, float]
    """First-order failure probability of each constraint, ``Phi(-beta)``."""

    system: float
    """Probability that *any* constraint fails, with the correlation kept."""

    independent_bound: float
    """The same quantity assuming the constraints are independent.

    Reported beside :attr:`system` because the gap between them is the price of
    the approximation a per-constraint margin makes.
    """

    @property
    def system_beta(self) -> float:
        """The system reliability index, ``-Phi^-1(P_f)``."""
        from scipy.stats import norm

        return float(norm.isf(min(max(self.system, 1e-16), 1.0 - 1e-16)))

    def binding(self) -> str:
        """The constraint contributing most of the failure probability."""
        return max(self.per_constraint, key=lambda name: self.per_constraint[name])


def failure_probability(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    band: dict[str, float] | None = None,
) -> Reliability | None:
    """Probability of failure by FORM, with the constraint correlation kept.

    This replaces the fixed-margin formulation ``g + k sigma <= 0`` that an
    earlier version of this module used, and the reason is worth stating.

    A margin of ``k`` standard deviations on *each* constraint separately is a
    reliability statement only if the constraints are independent.  They are
    emphatically not: every one of them is a function of the same eleven
    dimensions, and :attr:`ConstraintMoments.correlation` measures pairs at
    ``0.94`` and at exactly ``-1``.  Requiring all eight to hold at three sigma
    simultaneously is therefore much stronger than requiring the *system* to be
    reliable at three sigma, and it buys that strength by refusing designs that
    are in fact acceptable.

    What is computed instead is the thing actually wanted:

    .. math::
        P_f = 1 - \\Phi_n(\\beta; \\rho), \\qquad \\beta_i = -g_i / \\sigma_i

    the probability that *any* constraint is missed, from the multivariate
    normal orthant with the correlation in place.  ``Phi_n`` is evaluated by
    SciPy's implementation of Genz's algorithm, which is randomised
    quasi-Monte Carlo: repeated calls agree to about seven significant figures
    rather than exactly, so anything comparing two evaluations should allow for
    that.

    First order means the constraint surfaces are linearised, which is the same
    approximation the exact Jacobians already make cheap, and the same one
    :func:`tolerance_report` checks against sampling.

    Args:
        design: The design to assess.
        grade: ISO 286 IT grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        The reliability, or ``None`` if the design cannot be analysed.
    """
    from scipy.stats import multivariate_normal, norm

    moments = constraint_moments(design, grade, angular, samples, targets, spec, band)
    if moments is None:
        return None

    beta = moments.beta
    per_constraint = {
        name: float(norm.sf(value)) for name, value in zip(moments.names, beta, strict=True)
    }
    independent = float(1.0 - np.prod([1.0 - value for value in per_constraint.values()]))

    finite = np.isfinite(beta)
    if not np.any(finite):
        return Reliability(moments, per_constraint, 0.0, independent)

    # The correlation is singular by construction -- the two sides of each
    # relaxed equality are perfectly anti-correlated -- so a ridge is added
    # before the orthant integral.  It is small enough not to move the answer
    # and large enough to keep the factorisation well posed.
    reduced = moments.correlation[np.ix_(finite, finite)]
    ridge = reduced + 1.0e-8 * np.eye(reduced.shape[0])
    try:
        safe = float(
            multivariate_normal(mean=np.zeros(ridge.shape[0]), cov=ridge).cdf(beta[finite])
        )
    except Exception:  # a degenerate correlation falls back to the bound
        safe = float(np.prod([1.0 - value for value in per_constraint.values()]))
    system = float(min(max(1.0 - safe, 0.0), 1.0))
    return Reliability(moments, per_constraint, system, independent)


def format_reliability(
    reliability: Reliability, target: float = TARGET_FAILURE_PROBABILITY
) -> str:
    """Render a :class:`Reliability` as an aligned table."""
    lines = ["reliability", "===========", ""]
    lines.append(f"  {'constraint':<16}{'g':>11}{'sigma':>11}{'beta':>8}{'P(fail)':>12}")
    moments = reliability.moments
    for name, value, sigma in zip(moments.names, moments.value, moments.sigma, strict=True):
        beta = -value / sigma if sigma > 0 else float("inf")
        lines.append(
            f"  {name:<16}{value:>11.5f}{sigma:>11.5f}{beta:>8.2f}"
            f"{reliability.per_constraint[name]:>12.3e}"
        )
    lines.append("")
    lines.append(f"  system P(fail), correlation kept   {reliability.system:.4e}")
    lines.append(f"  the same assuming independence     {reliability.independent_bound:.4e}")
    lines.append(f"  system reliability index beta      {reliability.system_beta:.2f}")
    lines.append(f"  target                             {target:.1e}")
    lines.append(f"  binding constraint                 {reliability.binding()}")
    return "\n".join(lines)


def required_bound(
    design: Design,
    constraint: str = "tdc_gap",
    target: float = TARGET_FAILURE_PROBABILITY,
    grade: int = DEFAULT_GRADE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
) -> float:
    """The bound a constraint needs so its failure probability meets the target.

    Inverts the first-order reliability relation: for a target ``p``, the
    nominal value must sit ``beta = -Phi^-1(p)`` standard deviations inside the
    bound, so the bound has to be at least ``current + beta * sigma`` above the
    nominal quantity.

    This answers the practical question the tolerance study raises and cannot
    settle -- *how much* would the specification have to give? -- with a number
    rather than "more".

    Args:
        design: The design to assess.
        constraint: Which of :data:`RELIABILITY_NAMES` to invert.
        target: Failure probability to design to.
        grade: ISO 286 IT grade.
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        The required right-hand side, in the constraint's own units.

    Raises:
        ValueError: If the design cannot be analysed.
    """
    from scipy.stats import norm

    moments = constraint_moments(design, grade, samples=samples, targets=targets, spec=spec)
    if moments is None:
        msg = "cannot assess a design that does not close"
        raise ValueError(msg)

    index = moments.names.index(constraint)
    beta = float(norm.isf(target))
    current = {
        "tdc_gap": targets.max_tdc_gap,
        "rod_angle": targets.max_rod_angle,
        "compatibility": targets.max_transmission,
        "side_load": targets.max_side_load,
    }.get(constraint, 0.0)
    # g = quantity - bound, so quantity = g + bound; the new bound must exceed
    # that quantity by beta sigma.
    quantity = moments.value[index] + current
    return float(quantity + beta * moments.sigma[index])


class FailureProbabilityDiscipline(Discipline):
    """System probability of failure, as a GEMSEO constraint output.

    Publishes the FORM estimate so an optimizer can be given
    ``P_f <= target`` -- a single, correlation-aware reliability requirement --
    instead of a fixed margin on each constraint separately.

    First order is what makes this affordable inside an optimization: it costs
    one Jacobian evaluation of a geometric analysis, so it can be evaluated at
    every iteration.  Its error against sampling is real and is documented on
    :func:`failure_probability`; the sampling estimate is the reference, not
    this.

    Args:
        grade: ISO 286 IT grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        step: Relative finite-difference step for the Jacobian.
        name: Discipline name.
    """

    auto_detect_grammar_files: ClassVar[bool] = False

    def __init__(
        self,
        grade: int = DEFAULT_GRADE,
        angular: float = ANGULAR_TOLERANCE,
        samples: int = 360,
        targets: DesignTargets = DEFAULT_TARGETS,
        spec: EngineSpec = DEFAULT_SPEC,
        step: float = 1.0e-5,
        name: str = "",
    ) -> None:
        from .reference import PUBLISHED_DESIGN

        super().__init__(name=name)
        self.grade = grade
        self.angular = angular
        self.samples = samples
        self.targets = targets
        self.spec = spec
        self.step = step
        self.input_grammar.update_from_names(VARIABLE_NAMES)
        self.output_grammar.update_from_names(["failure_probability", "system_beta"])
        self.default_input_data = PUBLISHED_DESIGN.to_mapping()

    def _values(self, design: Design) -> tuple[float, float]:
        reliability = failure_probability(
            design, self.grade, self.angular, self.samples, self.targets, self.spec
        )
        if reliability is None:
            return 1.0, -10.0
        return reliability.system, reliability.system_beta

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        probability, beta = self._values(Design.from_mapping(dict(input_data)))
        return {
            "failure_probability": np.array([probability]),
            "system_beta": np.array([beta]),
        }

    def _compute_jacobian(
        self,
        input_names: Sequence[str] = (),
        output_names: Sequence[str] = (),
    ) -> None:
        self._init_jacobian(input_names, output_names)
        base = Design.from_mapping(dict(self.io.data)).to_array()
        rows: dict[str, list[float]] = {"failure_probability": [], "system_beta": []}
        for index, value in enumerate(base):
            delta = self.step * max(abs(float(value)), 1.0)
            ahead, behind = base.copy(), base.copy()
            ahead[index] += delta
            behind[index] -= delta
            forward = self._values(Design.from_array(ahead))
            backward = self._values(Design.from_array(behind))
            rows["failure_probability"].append((forward[0] - backward[0]) / (2.0 * delta))
            rows["system_beta"].append((forward[1] - backward[1]) / (2.0 * delta))
        for output, values in rows.items():
            if output not in self.jac:
                continue
            for index, variable in enumerate(VARIABLE_NAMES):
                if variable in self.jac[output]:
                    self.jac[output][variable] = np.array([[values[index]]])
