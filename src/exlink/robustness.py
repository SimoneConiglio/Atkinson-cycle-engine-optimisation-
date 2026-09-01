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


def _constraint_vector(
    design: Design,
    samples: int,
    targets: DesignTargets,
    spec: EngineSpec,
) -> FloatArray | None:
    """All seven constraints as one vector, negative meaning satisfied."""
    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        return None
    equality = equality_constraints(analysis, targets)
    inequality = inequality_constraints(analysis, targets)
    banded = np.abs(equality) - np.array(
        [EQUALITY_BAND["expansion_stroke"], EQUALITY_BAND["compression_ratio"]]
    )
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


#: The constraints worth robustifying, and the Jacobian row that carries each.
#: ``clearance`` is left out deliberately: its capability is above 400, so a
#: robust margin on it would only add cost, and it is the one constraint here
#: without an analytic gradient.
ROBUST_SOURCES: dict[str, str] = {
    "rod_angle_margin": "rod_angle_margin",
    "compatibility_margin": "compatibility_margin",
    "tdc_gap_margin": "tdc_gap_margin",
    "side_load_margin": "side_load_margin",
    "stroke_band_upper": "stroke_error",
    "stroke_band_lower": "stroke_error",
    "ratio_band_upper": "compression_ratio_error",
    "ratio_band_lower": "compression_ratio_error",
}

ROBUST_NAMES: tuple[str, ...] = tuple(f"{name}_robust" for name in ROBUST_SOURCES)
"""The robust counterparts, each ``<= 0`` when the design is robustly feasible."""

DEFAULT_SIGMA_LEVEL = 3.0
"""Standard deviations of margin the robust constraints hold.

``g + k sigma_g <= 0`` at ``k = 3`` puts the nominal design three standard
deviations inside the constraint, i.e. a one-sided process capability of 1.0.
The usual industrial target is 1.33, which is ``k = 4``; that is available and
is what a production drawing would ask for, but on this problem it is
unreachable for the top-dead-centre gap at any geometry -- see
:func:`required_grade`.
"""


def robust_margins(
    design: Design,
    sigma_level: float = DEFAULT_SIGMA_LEVEL,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
) -> dict[str, float]:
    """The robust counterpart of each geometric constraint.

    Replaces ``g(X) <= 0`` by

    .. math:: g(X) + k \\sqrt{\\nabla g^\\top \\Sigma \\nabla g} \\le 0

    so the optimizer is required to hold the constraint not at the nominal
    design but ``k`` standard deviations inside it, given the manufacturing
    covariance of the parts.  This is what
    :func:`tolerance_report` measures *after* the fact, moved into the
    formulation so the optimizer has to pay for it while choosing.

    Every gradient here is the exact one from :mod:`exlink.jacobian`, so the
    robust margin costs one Jacobian evaluation on top of the analysis -- the
    same observation that made the tolerance study cheap.

    Args:
        design: The design to assess.
        sigma_level: ``k``, the standard deviations of margin required.
        grade: ISO 286 IT grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.

    Returns:
        ``{name_robust: value}`` over :data:`ROBUST_NAMES`, each non-positive
        when the design is robustly feasible.  A design that cannot be analysed
        returns large positive values rather than raising, so an optimizer can
        walk through it.
    """
    from .scenarios import DEFAULT_EQUALITY_TOLERANCE

    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        return dict.fromkeys(ROBUST_NAMES, 1.0e3)

    kinematic = kinematic_jacobian(design, analysis.require_solved().kinematics, spec)
    rows = metric_jacobian(design, analysis, kinematic, spec)
    sigma_matrix = covariance(design, grade, angular)

    metrics = analysis.metrics
    stroke = metrics.expansion_stroke - targets.expansion_stroke
    ratio = metrics.compression_ratio - targets.compression_ratio
    nominal = {
        "rod_angle_margin": metrics.rod_angle - targets.max_rod_angle,
        "compatibility_margin": metrics.compatibility - targets.max_transmission,
        "tdc_gap_margin": metrics.tdc_gap - targets.max_tdc_gap,
        "side_load_margin": metrics.side_load_ratio - targets.max_side_load,
        "stroke_band_upper": stroke - DEFAULT_EQUALITY_TOLERANCE["stroke_error"],
        "stroke_band_lower": -stroke - DEFAULT_EQUALITY_TOLERANCE["stroke_error"],
        "ratio_band_upper": (ratio - DEFAULT_EQUALITY_TOLERANCE["compression_ratio_error"]),
        "ratio_band_lower": (-ratio - DEFAULT_EQUALITY_TOLERANCE["compression_ratio_error"]),
    }

    result: dict[str, float] = {}
    for name, source in ROBUST_SOURCES.items():
        gradient = np.asarray(rows[source], dtype=float)
        if name.endswith("_lower"):
            gradient = -gradient
        deviation = float(np.sqrt(max(gradient @ sigma_matrix @ gradient, 0.0)))
        result[f"{name}_robust"] = float(nominal[name]) + sigma_level * deviation
    return result


def robust_report(
    design: Design,
    sigma_level: float = DEFAULT_SIGMA_LEVEL,
    grade: int = DEFAULT_GRADE,
    samples: int = 360,
) -> str:
    """Render the robust margins beside the nominal ones."""
    from .model import analyse as _analyse

    robust = robust_margins(design, sigma_level, grade, samples=samples)
    analysis = _analyse(design, samples=samples)
    title = f"robust margins at IT{grade}, k = {sigma_level:g}"
    lines = [title, "=" * len(title), ""]
    lines.append(f"  {'constraint':<24}{'nominal':>12}{'robust':>12}   holds")
    nominal_rows = dict(
        zip(
            ("rod_angle_margin", "compatibility_margin", "tdc_gap_margin", "side_load_margin"),
            inequality_constraints(analysis)[[0, 1, 2, 4]],
            strict=True,
        )
    )
    for name in ROBUST_SOURCES:
        value = robust[f"{name}_robust"]
        shown = nominal_rows.get(name)
        nominal = f"{shown:12.5f}" if shown is not None else " " * 12
        lines.append(f"  {name:<24}{nominal}{value:>12.5f}   {'yes' if value <= 0 else 'NO'}")
    return "\n".join(lines)


class RobustMarginDiscipline(Discipline):
    """The robust counterparts of the geometric constraints, as a GEMSEO discipline.

    Publishes ``g + k sigma_g`` for each constraint in :data:`ROBUST_SOURCES`,
    so an optimizer can be *required* to hold every constraint ``k`` standard
    deviations inside its bound given the manufacturing covariance -- rather
    than being told afterwards that its answer does not survive tolerance.

    Derivatives are finite differences, deliberately.  The exact route would
    need second derivatives of ``g``: the robust margin contains
    ``sqrt(grad g' Sigma grad g)``, whose gradient carries a Hessian this
    package does not compute.  The usual first-order dodge is to hold ``sigma``
    locally constant and reuse ``grad g``, but that discards precisely the term
    that makes a robust optimum differ from a nominal one -- the pull towards
    regions where the constraint is *less* sensitive.  Since the underlying
    analysis is geometric and costs about ten milliseconds, differencing the
    whole robust margin is affordable and keeps that term.

    Args:
        sigma_level: ``k``, the standard deviations of margin required.
        grade: ISO 286 IT grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        step: Relative finite-difference step.
        name: Discipline name.
    """

    auto_detect_grammar_files: ClassVar[bool] = False

    def __init__(
        self,
        sigma_level: float = DEFAULT_SIGMA_LEVEL,
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
        self.sigma_level = sigma_level
        self.grade = grade
        self.angular = angular
        self.samples = samples
        self.targets = targets
        self.spec = spec
        self.step = step
        self.input_grammar.update_from_names(VARIABLE_NAMES)
        self.output_grammar.update_from_names(list(ROBUST_NAMES))
        self.default_input_data = PUBLISHED_DESIGN.to_mapping()

    def _margins(self, design: Design) -> dict[str, float]:
        return robust_margins(
            design,
            self.sigma_level,
            self.grade,
            self.angular,
            self.samples,
            self.targets,
            self.spec,
        )

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        values = self._margins(Design.from_mapping(dict(input_data)))
        return {name: np.array([values[name]]) for name in ROBUST_NAMES}

    def _compute_jacobian(
        self,
        input_names: Sequence[str] = (),
        output_names: Sequence[str] = (),
    ) -> None:
        self._init_jacobian(input_names, output_names)
        base = Design.from_mapping(dict(self.io.data)).to_array()
        gradients: dict[str, list[float]] = {name: [] for name in ROBUST_NAMES}
        for index, value in enumerate(base):
            delta = self.step * max(abs(float(value)), 1.0)
            ahead, behind = base.copy(), base.copy()
            ahead[index] += delta
            behind[index] -= delta
            forward = self._margins(Design.from_array(ahead))
            backward = self._margins(Design.from_array(behind))
            for name in ROBUST_NAMES:
                gradients[name].append((forward[name] - backward[name]) / (2.0 * delta))
        for output in ROBUST_NAMES:
            if output not in self.jac:
                continue
            for index, variable in enumerate(VARIABLE_NAMES):
                if variable in self.jac[output]:
                    self.jac[output][variable] = np.array([[gradients[output][index]]])
