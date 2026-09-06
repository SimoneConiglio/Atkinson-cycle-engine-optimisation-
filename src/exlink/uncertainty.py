"""Reliability against more than the dimensions.

:mod:`exlink.robustness` propagates ISO 286 dimensional tolerances through the
*geometric* constraints, exactly, using the analytic Jacobians.  That is the
right thing to have inside an optimizer -- it costs one Jacobian and no extra
analyses -- but it leaves five of the thirteen constraints without a
probability, and §7.2 records why: they are not functions of the eleven
dimensions alone.  ``saturation``, ``slenderness`` and ``bearing`` depend on the
material the members are made of; ``runs`` depends on how much of the indicated
work friction takes; and all of them depend on how hard the charge burns.

This module widens the uncertain vector to seventeen entries -- the eleven
dimensions plus six parameters that are not design variables and are not
constants either -- and prices every constraint against all of them.

What the six parameters are, and why these values
-------------------------------------------------
Coefficients of variation, taken as the JCSS *Probabilistic Model Code* and the
usual structural-reliability literature give them for hot-rolled steel, with two
that are this problem's rather than the material's:

==================  ======  ================================================
parameter           CoV     why
==================  ======  ================================================
yield strength      0.07    JCSS steel; the dominant material scatter
ultimate strength   0.05    correlated with yield in reality; see below
Young's modulus     0.03    JCSS; tight, and it only enters buckling
density             0.01    tighter still; a bar is a bar
friction coeff.     0.20    §7.2 puts the FMEP uncertainty near 30 %; this is
                            the part of it attributable to ``mu`` alone
explosion ratio     0.10    cycle-to-cycle variation of a single-cylinder
                            engine, plus calibration of an idealised burn
==================  ======  ================================================

Taking the two strengths independent is conservative: on a real heat they move
together, and a build with a low yield usually has a low ultimate, so the joint
probability of both being bad is smaller than the product used here.

Why finite differences here and analytic Jacobians there
--------------------------------------------------------
Four of the five newly covered constraints are outputs of the sizing/inertia
fixed point, and one is an output of the vehicle sub-problem beyond it.  There
is no closed form for ``d(bearing load)/d(yield strength)`` that does not amount
to differentiating the whole MDA, so this module differences it: one full
evaluation per uncertain entry, eighteen in all.  That is far too expensive to
put inside an optimizer, which is exactly why the fast analytic path in
:mod:`exlink.robustness` still exists and is still what the RBDO loop calls.
This is the wider, slower picture the fast one is checked against.

The step sizes matter and are not uniform: a dimensional step is absolute
(microns, matched to the tolerance being propagated) and a parameter step is
relative (a fraction of the nominal), because the two have no common scale.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import VARIABLE_NAMES, Design
from .dynamics import DEFAULT_SPEED_RPM
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors
from .performance import Performance, evaluate
from .robustness import (
    ANGULAR_TOLERANCE,
    DEFAULT_GRADE,
    RELIABILITY_NAMES,
    SIGMA_PER_HALF_WIDTH,
    ConstraintMoments,
    Reliability,
    _band_widths,
    constraint_jacobian,
    reliability_from_moments,
    tolerance_half_widths,
)
from .sizing import MAX_DIAMETER, MEMBER_IS_SLENDER, member_lengths

FloatArray = NDArray[np.float64]

MAX_SLENDERNESS = 0.34
"""Largest ``d / L`` a connecting link may reach before it is not a rod."""

SATURATION_FRACTION = 0.98
"""Fraction of :data:`~exlink.sizing.MAX_DIAMETER` counted as run away."""

MAX_WIDTH_FACTOR = 12.0
"""Largest gear face width as a multiple of the module."""

COUPLED_NAMES: tuple[str, ...] = ("saturation", "slenderness", "bearing", "runs", "gear")
"""The five constraints the dimensional-only model cannot price."""

WIDENED_NAMES: tuple[str, ...] = (*RELIABILITY_NAMES, *COUPLED_NAMES)
"""All thirteen, in the order this module's vectors use."""


@dataclass(frozen=True)
class ParameterScatter:
    """Coefficients of variation of the six non-dimensional uncertainties."""

    yield_strength: float = 0.07
    ultimate_strength: float = 0.05
    youngs_modulus: float = 0.03
    density: float = 0.01
    friction: float = 0.20
    explosion_ratio: float = 0.10

    def as_array(self) -> FloatArray:
        """The six, in the order :data:`PARAMETER_NAMES` uses."""
        return np.array(
            [
                self.yield_strength,
                self.ultimate_strength,
                self.youngs_modulus,
                self.density,
                self.friction,
                self.explosion_ratio,
            ]
        )


PARAMETER_NAMES: tuple[str, ...] = (
    "yield_strength",
    "ultimate_strength",
    "youngs_modulus",
    "density",
    "friction",
    "explosion_ratio",
)

DEFAULT_SCATTER = ParameterScatter()
"""The scatter the study reports against."""

UNCERTAIN_NAMES: tuple[str, ...] = (*VARIABLE_NAMES, *PARAMETER_NAMES)
"""The seventeen uncertain entries, dimensions first."""


def _perturbed(
    parameter: str,
    factor: float,
    material: Material,
    spec: EngineSpec,
) -> tuple[Material, EngineSpec, float]:
    """Apply a multiplicative perturbation to one parameter.

    Returns:
        The perturbed material, the perturbed engine data, and a multiplier on
        the friction coefficients (which live in neither).
    """
    friction = 1.0
    if parameter == "yield_strength":
        material = replace(material, yield_strength=material.yield_strength * factor)
    elif parameter == "ultimate_strength":
        material = replace(material, ultimate_strength=material.ultimate_strength * factor)
    elif parameter == "youngs_modulus":
        material = replace(material, youngs_modulus=material.youngs_modulus * factor)
    elif parameter == "density":
        material = replace(material, density=material.density * factor)
    elif parameter == "friction":
        friction = factor
    elif parameter == "explosion_ratio":
        spec = replace(spec, explosion_ratio=spec.explosion_ratio * factor)
    else:  # pragma: no cover - guarded by PARAMETER_NAMES
        msg = f"unknown uncertain parameter {parameter!r}"
        raise ValueError(msg)
    return material, spec, friction


def constraint_vector(
    performance: Performance,
    targets: DesignTargets = DEFAULT_TARGETS,
    band: dict[str, float] | None = None,
) -> FloatArray | None:
    """All thirteen constraints of one evaluation, negative meaning satisfied.

    The eight of :data:`~exlink.robustness.RELIABILITY_NAMES` in their usual
    form, then the five that only exist once the loads are dynamic and the car
    is in the problem.  The two sign conventions of :mod:`exlink.scenarios` are
    reconciled here: everything below is a violation when positive.

    Args:
        performance: A solved evaluation.
        targets: Constraint right-hand sides.
        band: Half-widths of the two relaxed equalities.

    Returns:
        Thirteen values, or ``None`` if the design was never sized.
    """
    coupled = performance.coupled
    friction = performance.friction
    if coupled is None or friction is None or not performance.analysis.valid:
        return None

    metrics = performance.metrics
    stroke = metrics.expansion_stroke - targets.expansion_stroke
    ratio = metrics.compression_ratio - targets.compression_ratio
    band_stroke, band_ratio = (float(width) for width in _band_widths(band))

    diameters = np.array(
        [coupled.diameters[name] for name in coupled.diameters],
        dtype=float,
    )
    lengths = member_lengths(performance.design)
    slenderness = float(np.max((diameters / lengths)[MEMBER_IS_SLENDER], initial=0.0))
    peak_bearing = float(np.max(np.linalg.norm(coupled.loads.reaction["R1"], axis=1)))
    pair = performance.budget.gears

    return np.array(
        [
            metrics.rod_angle - targets.max_rod_angle,
            metrics.compatibility - targets.max_transmission,
            metrics.tdc_gap - targets.max_tdc_gap,
            metrics.side_load_ratio - targets.max_side_load,
            stroke - band_stroke,
            -stroke - band_stroke,
            ratio - band_ratio,
            -ratio - band_ratio,
            float(diameters.max()) - SATURATION_FRACTION * MAX_DIAMETER,
            slenderness - MAX_SLENDERNESS,
            peak_bearing / targets.max_bearing_load - 1.0,
            # The two margins negated: positive means violated, as above.
            -friction.brake_work / max(friction.indicated_work, 1.0),
            (pair.width_factor - MAX_WIDTH_FACTOR if pair is not None else 1.0),
        ]
    )


def widened_moments(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    scatter: ParameterScatter = DEFAULT_SCATTER,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    band: dict[str, float] | None = None,
    relative_step: float = 1.0e-3,
    **kwargs: Any,
) -> tuple[ConstraintMoments, VarianceShares] | None:
    """First-order moments of all thirteen constraints, over all seventeen inputs.

    Differences one full evaluation per uncertain entry, so it costs eighteen
    analyses -- seconds, not milliseconds.  See the module docstring for why
    that is the right trade here and the wrong one inside an optimizer.

    Args:
        design: The design to assess.
        grade: ISO 286 IT grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        scatter: Coefficients of variation of the six parameters.
        speed_rpm: Analysis speed [rev/min].
        samples: Crank angles per revolution.  The study's 360, not less: the
            expansion stroke and the top-dead-centre gap are differences of
            nearby extrema, and their *analytic* gradients move by tens of per
            cent between 180 and 360 stations.  Eighteen analyses at 360 cost
            about a minute, which is the price of the two constraints that
            bind agreeing with section 6.2.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        material: The nominal material.
        safety: The design factors.
        band: Half-widths of the two relaxed equalities.
        relative_step: Differencing step on the six parameters, as a fraction.
        **kwargs: Forwarded to :func:`exlink.performance.evaluate`.

    Returns:
        The moments over thirteen constraints and the split of each one's
        variance by source, or ``None`` if the nominal design cannot be
        analysed or sized.
    """
    common = dict(
        speed_rpm=speed_rpm,
        samples=samples,
        spec=spec,
        material=material,
        safety=safety,
        **kwargs,
    )
    nominal = evaluate(design, **common)  # type: ignore[arg-type]
    base = constraint_vector(nominal, targets, band)
    if base is None:
        return None

    half_widths = tolerance_half_widths(design, grade, angular)
    sigma_design = half_widths / SIGMA_PER_HALF_WIDTH
    sigma_parameter = scatter.as_array()

    columns: list[FloatArray] = []

    # The eleven dimensions: an absolute step, one tenth of the tolerance being
    # propagated, so the difference is taken on the scale the answer is about.
    for index, name in enumerate(VARIABLE_NAMES):
        step = max(0.1 * float(half_widths[index]), 1.0e-6)
        moved = replace(design, **{name: getattr(design, name) + step})
        shifted = constraint_vector(evaluate(moved, **common), targets, band)  # type: ignore[arg-type]
        columns.append(np.zeros_like(base) if shifted is None else (shifted - base) / step)

    # The six parameters: a relative step, since they share no unit.
    for name in PARAMETER_NAMES:
        moved_material, moved_spec, friction = _perturbed(
            name, 1.0 + relative_step, material, spec
        )
        extra = dict(common)
        extra["material"] = moved_material
        extra["spec"] = moved_spec
        if friction != 1.0:
            extra["friction_scale"] = friction
        shifted = constraint_vector(evaluate(design, **extra), targets, band)  # type: ignore[arg-type]
        columns.append(
            np.zeros_like(base) if shifted is None else (shifted - base) / relative_step
        )

    jacobian = np.stack(columns, axis=1)

    # Where an exact derivative exists, use it.  The eight geometric
    # constraints have analytic gradients with respect to the eleven
    # dimensions, and differencing them instead would put a step-size error
    # into the one part of this calculation that does not need one -- the
    # stroke and the top-dead-centre gap are near-cancellations, and a
    # micron-scale difference of them carries tens of per cent.  The
    # widened model then contains the narrow one exactly rather than
    # approximately.
    exact = constraint_jacobian(design, samples, targets, spec, band)
    if exact is not None:
        _value, rows = exact
        jacobian[: rows.shape[0], : len(VARIABLE_NAMES)] = rows

    sigma_vector = np.concatenate([sigma_design, sigma_parameter])
    covariances = (jacobian * sigma_vector**2) @ jacobian.T
    sigma = np.sqrt(np.clip(np.diag(covariances), 0.0, None))
    outer = np.outer(sigma, sigma)
    correlation = np.divide(covariances, outer, out=np.eye(sigma.size), where=outer > 0.0)
    return ConstraintMoments(
        names=WIDENED_NAMES,
        value=base,
        sigma=sigma,
        correlation=np.clip(correlation, -1.0, 1.0),
    ), _shares(jacobian, sigma_vector)


def _shares(jacobian: FloatArray, sigma_vector: FloatArray) -> VarianceShares:
    """Split each constraint's first-order variance by source.

    First order and with the entries independent, the variance is a plain sum
    of ``(dg/du_j sigma_j)^2``, so the split is exact rather than an
    attribution -- there is no interaction term to allocate.

    Args:
        jacobian: ``(n_constraints, 17)`` derivatives.
        sigma_vector: The seventeen standard deviations.

    Returns:
        The pooled shares.
    """
    terms = (jacobian * sigma_vector) ** 2
    split = len(VARIABLE_NAMES)
    pooled = np.column_stack([terms[:, :split].sum(axis=1), terms[:, split:]])
    total = pooled.sum(axis=1, keepdims=True)
    return VarianceShares(
        names=WIDENED_NAMES,
        sources=("dimensions", *PARAMETER_NAMES),
        shares=np.divide(pooled, total, out=np.zeros_like(pooled), where=total > 0.0),
    )


@dataclass(frozen=True)
class VarianceShares:
    """Where each constraint's variance comes from."""

    names: tuple[str, ...]
    """Constraint names, as :data:`WIDENED_NAMES`."""

    sources: tuple[str, ...]
    """``("dimensions", *PARAMETER_NAMES)`` -- seven sources, not seventeen.

    The eleven dimensions are pooled, because they are one modelling decision
    (an ISO grade) and are reported as one in section 3.10.  The six parameters
    are kept apart, because the point of widening the model is to find out which
    of them matters.
    """

    shares: FloatArray
    """``(n_constraints, 7)``, each row summing to one where sigma is positive."""

    def dominant(self) -> dict[str, str]:
        """The largest source of each constraint's variance."""
        return {
            name: self.sources[int(np.argmax(row))]
            for name, row in zip(self.names, self.shares, strict=True)
        }


def format_shares(shares: VarianceShares) -> str:
    """Render a :class:`VarianceShares` as an aligned percentage table."""
    header = "".join(f"{source[:9]:>11}" for source in shares.sources)
    lines = ["variance shares", "=" * (17 + len(header)), "", f"  {'constraint':<15}{header}"]
    for name, row in zip(shares.names, shares.shares, strict=True):
        cells = "".join(f"{100.0 * value:>10.1f}%" for value in row)
        lines.append(f"  {name:<15}{cells}")
    return "\n".join(lines)


def widened_reliability(design: Design, **kwargs: Any) -> Reliability | None:
    """The system failure probability over all thirteen constraints.

    Args:
        design: The design to assess.
        **kwargs: Forwarded to :func:`widened_moments`.

    Returns:
        The reliability, or ``None`` if the design cannot be assessed.
    """
    outcome = widened_moments(design, **kwargs)
    if outcome is None:
        return None
    return reliability_from_moments(outcome[0])
