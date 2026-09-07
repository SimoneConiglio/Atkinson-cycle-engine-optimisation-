"""Second derivatives, and the measurement that they are the wrong instrument.

What this module was built for
------------------------------
§6.3 asked for :math:`\\nabla^2 g`, for two reasons that look like one problem.
The reliability estimate is first order, and where a constraint surface curves
towards the design that under-predicts; and the derivative of
:math:`\\beta_i = -g_i/\\sigma_i` is a finite difference, because
:math:`\\partial\\sigma_i/\\partial x` needs the Hessian.  §5.3 had just measured
the first-order estimate as optimistic by a factor of seven at the study's
result, so curvature was the obvious suspect.

What it found instead
---------------------
The Hessian obtained by differencing the *analytic* gradient -- which is the
well-conditioned route, since the first derivatives are exact -- does not exist.
Its norm scales as :math:`1/h` over two decades of step size and it is 100 %
asymmetric at every step, which is the signature of a **constant jump in the
gradient** rather than of a second derivative.  A finer crank-angle grid does
not change it, so it is not discretisation.

The jump is a branch switch.  The expansion stroke is a maximum over the two
top dead centres, and the design §5.5 arrives at has them 0.107 um
apart, so perturbing any dimension by a fraction of a micron swaps which one
attains the maximum and moves :math:`\\partial STE/\\partial a` by a factor of 56.
The surface is not curved, it is **kinked**, and no order of Taylor expansion
repairs a kink.  What was needed was not a second derivative but a first one, of
the other branch -- which is what :mod:`exlink.robustness` now carries.

So the contents here are a diagnostic rather than a correction:
:func:`differencing_scaling` produces the table that establishes the noise
floor, :func:`constraint_hessian` and :func:`hessian_asymmetry` are what it is
built from, and :func:`second_order_reliability` applies Breitung's curvature
correction and reports -- honestly -- that its expansion does not hold at this
design.  :func:`beta_gradient` is the closed-form reliability derivative the
entry also asked for; it is correct wherever the Hessian is, which is to say
away from the tie, and is offered on that condition.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import VARIABLE_NAMES, Design
from .robustness import (
    ANGULAR_TOLERANCE,
    DEFAULT_GRADE,
    RELIABILITY_NAMES,
    SIGMA_PER_HALF_WIDTH,
    ConstraintMoments,
    constraint_jacobian,
    covariance,
    reliability_from_moments,
    tolerance_half_widths,
)

FloatArray = NDArray[np.float64]

DEFAULT_STEP_FRACTION = 0.05
"""Differencing step, as a fraction of each variable's tolerance half-width.

Scaled to the tolerance rather than to the variable, because the derivative
wanted is the one that governs behaviour *over the scatter*: a step much
smaller than a micron differences numerical noise in the extremum refinement,
and one much larger than the tolerance measures the wrong neighbourhood.
"""

MIN_CURVATURE_PRODUCT = 1.0e-6
"""Floor on ``1 + beta kappa``, below which Breitung's formula has broken down.

The expansion is valid while every factor is positive.  A curvature that drives
one negative means the limit surface curves back on itself within a
:math:`\\beta` of the design point, and the answer there is a sampled one, not
a corrected first-order one.
"""


@dataclass(frozen=True)
class SecondOrder:
    """A reliability estimate corrected for constraint curvature."""

    names: tuple[str, ...]
    beta: FloatArray
    """First-order index of each constraint."""

    first_order: dict[str, float]
    """``Phi(-beta)`` for each constraint."""

    second_order: dict[str, float]
    """The same, with Breitung's curvature correction."""

    curvature: dict[str, float]
    """``prod (1 + beta kappa)^(-1/2)``: the factor the correction applies."""

    valid: dict[str, bool]
    """Whether Breitung's expansion held for that constraint."""

    system: float
    """System probability with every constraint corrected."""

    first_order_system: float
    """The uncorrected system probability, for comparison."""

    asymmetry: float
    """Largest relative asymmetry of any computed Hessian; the differencing error."""

    @property
    def ratio(self) -> float:
        """How much the correction moved the system probability."""
        if self.first_order_system <= 0.0:
            return float("inf")
        return self.system / self.first_order_system


def constraint_hessian(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    band: dict[str, float] | None = None,
    step_fraction: float = DEFAULT_STEP_FRACTION,
) -> tuple[FloatArray, FloatArray, FloatArray, float] | None:
    """Second derivatives of the eight geometric constraints.

    Args:
        design: The design to differentiate at.
        grade: ISO 286 grade, which sets the differencing step.
        angular: Angular half-width [deg], likewise.
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        band: Half-widths of the two relaxed equalities.
        step_fraction: Step as a fraction of each tolerance half-width.

    Returns:
        ``(value, jacobian, hessian, asymmetry)`` of shapes ``(8,)``,
        ``(8, 11)``, ``(8, 11, 11)`` and a scalar, or ``None`` if the design
        cannot be analysed.  The asymmetry is measured *before* symmetrising
        and is the honest error bar on the differencing.
    """
    centre = constraint_jacobian(design, samples, targets, spec, band)
    if centre is None:
        return None
    value, jacobian = centre

    widths = tolerance_half_widths(design, grade, angular)
    steps = np.maximum(step_fraction * widths, 1.0e-7)

    columns = []
    for index, name in enumerate(VARIABLE_NAMES):
        moved = replace(design, **{name: getattr(design, name) + float(steps[index])})
        shifted = constraint_jacobian(moved, samples, targets, spec, band)
        if shifted is None:
            columns.append(np.zeros_like(jacobian))
            continue
        columns.append((shifted[1] - jacobian) / steps[index])

    # (8, 11, 11): axis 1 is the gradient index, axis 2 the differencing index.
    raw = np.stack(columns, axis=2)
    hessian = 0.5 * (raw + np.swapaxes(raw, 1, 2))
    return value, jacobian, hessian, hessian_asymmetry(raw)


def hessian_asymmetry(raw: FloatArray) -> float:
    """Largest relative asymmetry in a computed Hessian.

    The exact Hessian is symmetric, so ``|H - H^T| / |H|`` measures nothing but
    the differencing error and is the honest error bar on it.

    Args:
        raw: ``(n, k, k)`` Hessians, before symmetrisation.

    Returns:
        The largest relative asymmetry over the stack; zero for an empty one.
    """
    worst = 0.0
    for block in np.atleast_3d(raw):
        scale = float(np.max(np.abs(block)))
        if scale <= 0.0:
            continue
        worst = max(worst, float(np.max(np.abs(block - block.T))) / scale)
    return worst


def _curvature_factor(
    beta: float, gradient: FloatArray, hessian: FloatArray, sigma_x: FloatArray
) -> tuple[float, bool]:
    """Breitung's factor for one constraint.

    Work in the standard normal space ``u`` where ``x = x_0 + diag(sigma) u``.
    There the gradient is ``sigma * grad g`` and the Hessian
    ``sigma_i sigma_j H_ij``.  Rotate so the last axis is the normal at the
    most probable point, and the principal curvatures are the eigenvalues of
    the leading ``(n-1)x(n-1)`` block of the scaled Hessian.

    Args:
        beta: The constraint's reliability index.
        gradient: ``dg/dx``.
        hessian: ``d2g/dx2``.
        sigma_x: Standard deviation of each design variable.

    Returns:
        ``(factor, valid)``.
    """
    scaled_gradient = gradient * sigma_x
    norm = float(np.linalg.norm(scaled_gradient))
    if norm <= 0.0 or not np.isfinite(beta):
        return 1.0, True
    scaled_hessian = hessian * np.outer(sigma_x, sigma_x)

    # An orthonormal basis whose last vector is the unit normal.
    normal = scaled_gradient / norm
    basis = np.linalg.qr(np.column_stack([normal, np.eye(normal.size)]))[0]
    if float(basis[:, 0] @ normal) < 0.0:
        basis = -basis
    tangent = basis[:, 1:]

    block = tangent.T @ scaled_hessian @ tangent / norm
    curvatures = np.linalg.eigvalsh(0.5 * (block + block.T))
    factors = 1.0 + beta * curvatures
    if np.any(factors <= MIN_CURVATURE_PRODUCT):
        return 1.0, False
    return float(np.prod(factors**-0.5)), True


def second_order_reliability(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    band: dict[str, float] | None = None,
    step_fraction: float = DEFAULT_STEP_FRACTION,
) -> SecondOrder | None:
    """Correct the first-order probabilities for constraint curvature.

    Args:
        design: The design to assess.
        grade: ISO 286 grade of the machined dimensions.
        angular: Angular assembly half-width [deg].
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        band: Half-widths of the two relaxed equalities.
        step_fraction: Differencing step, as a fraction of the tolerance.

    Returns:
        The corrected estimate, or ``None`` if the design cannot be analysed.
    """
    from scipy.stats import norm as normal

    outcome = constraint_hessian(
        design, grade, angular, samples, targets, spec, band, step_fraction
    )
    if outcome is None:
        return None
    value, jacobian, hessian, asymmetry = outcome

    sigma_matrix = covariance(design, grade, angular)
    sigma = np.sqrt(np.clip(np.diag(jacobian @ sigma_matrix @ jacobian.T), 0.0, None))
    safe = np.where(sigma > 0.0, sigma, np.inf)
    beta = -value / safe
    sigma_x = tolerance_half_widths(design, grade, angular) / SIGMA_PER_HALF_WIDTH

    first: dict[str, float] = {}
    second: dict[str, float] = {}
    factors: dict[str, float] = {}
    valid: dict[str, bool] = {}
    for index, name in enumerate(RELIABILITY_NAMES):
        base = float(normal.sf(beta[index]))
        factor, ok = _curvature_factor(
            float(beta[index]), jacobian[index], hessian[index], sigma_x
        )
        first[name] = base
        factors[name] = factor
        valid[name] = ok
        second[name] = float(min(base * factor, 1.0))

    # The system probability, with each constraint's corrected index carried
    # through the same orthant integral the first-order route uses -- so the
    # correlation is still kept and only the margins move.
    corrected = np.array(
        [
            float(normal.isf(min(max(second[name], 1.0e-300), 1.0 - 1.0e-16)))
            for name in RELIABILITY_NAMES
        ]
    )
    moments = ConstraintMoments(
        names=RELIABILITY_NAMES,
        value=-corrected * sigma,
        sigma=sigma,
        correlation=_correlation(jacobian, sigma_matrix, sigma),
    )
    plain = ConstraintMoments(
        names=RELIABILITY_NAMES,
        value=value,
        sigma=sigma,
        correlation=moments.correlation,
    )
    return SecondOrder(
        names=RELIABILITY_NAMES,
        beta=beta,
        first_order=first,
        second_order=second,
        curvature=factors,
        valid=valid,
        system=reliability_from_moments(moments).system,
        first_order_system=reliability_from_moments(plain).system,
        asymmetry=asymmetry,
    )


def _correlation(
    jacobian: FloatArray, sigma_matrix: FloatArray, sigma: FloatArray
) -> FloatArray:
    """Constraint correlation, as :func:`exlink.robustness.constraint_moments` forms it."""
    covariances = jacobian @ sigma_matrix @ jacobian.T
    outer = np.outer(sigma, sigma)
    correlation = np.divide(covariances, outer, out=np.eye(sigma.size), where=outer > 0.0)
    return np.clip(correlation, -1.0, 1.0)


def beta_gradient(
    design: Design,
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    band: dict[str, float] | None = None,
    step_fraction: float = DEFAULT_STEP_FRACTION,
) -> FloatArray | None:
    """Exact derivative of each reliability index, ``(8, 11)``.

    The quantity §4.7's search is steered on is ``min_i beta_i``, and until
    now its derivative was a finite difference of a difference.  With the
    Hessian in hand it is closed form:

    .. math:: \\frac{\\partial \\beta_i}{\\partial x}
              = -\\frac{\\nabla g_i}{\\sigma_i}
                + \\frac{g_i}{\\sigma_i^{3}} \\nabla^2 g_i \\Sigma \\nabla g_i

    Args:
        design: The design to differentiate at.
        grade: ISO 286 grade.
        angular: Angular half-width [deg].
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        band: Half-widths of the two relaxed equalities.
        step_fraction: Step for the Hessian.

    Returns:
        ``d beta_i / d x_j``, or ``None`` if the design cannot be analysed.
    """
    outcome = constraint_hessian(
        design, grade, angular, samples, targets, spec, band, step_fraction
    )
    if outcome is None:
        return None
    value, jacobian, hessian, _asymmetry = outcome
    sigma_matrix = covariance(design, grade, angular)
    sigma = np.sqrt(np.clip(np.diag(jacobian @ sigma_matrix @ jacobian.T), 0.0, None))

    rows = np.zeros_like(jacobian)
    for index in range(jacobian.shape[0]):
        if sigma[index] <= 0.0:
            continue
        d_sigma = hessian[index] @ sigma_matrix @ jacobian[index] / sigma[index]
        rows[index] = (
            -jacobian[index] / sigma[index] + value[index] * d_sigma / sigma[index] ** 2
        )
    return rows


def _metric_gradient(
    design: Design, row: str, samples: int, spec: EngineSpec
) -> FloatArray | None:
    """One row of :func:`exlink.jacobian.metric_jacobian`, by name."""
    from .jacobian import kinematic_jacobian, metric_jacobian
    from .model import analyse

    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        return None
    kinematic = kinematic_jacobian(design, analysis.require_solved().kinematics, spec)
    return np.asarray(metric_jacobian(design, analysis, kinematic, spec)[row], dtype=float)


def differencing_scaling(
    design: Design,
    rows: tuple[str, ...] = ("stroke_error", "expansion_stroke_1", "expansion_stroke_2"),
    fractions: tuple[float, ...] = (0.5, 0.2, 0.1, 0.05, 0.02, 0.01),
    grade: int = DEFAULT_GRADE,
    angular: float = ANGULAR_TOLERANCE,
    samples: int = 360,
    spec: EngineSpec = DEFAULT_SPEC,
) -> dict[str, list[tuple[float, float, float]]]:
    """Norm and asymmetry of the differenced Hessian, against step size.

    The diagnostic of §5.3, and the reason it is worth having as a function
    rather than as a one-off script: a genuine second derivative is insensitive
    to the step over a range of steps and its matrix is symmetric, while a
    discontinuous gradient gives a norm proportional to ``1/h`` and an
    asymmetry of one.  Run on the three rows it defaults to, it separates the
    two cases side by side -- ``stroke_error`` is the maximum over the two top
    dead centres and jumps, the two branch rows are smooth and do not.

    Args:
        design: The design to differentiate at.
        rows: Jacobian row names to test.
        fractions: Steps to try, as fractions of the tolerance half-width.
        grade: ISO 286 grade, which sets the tolerance.
        angular: Angular half-width [deg].
        samples: Crank angles per revolution.
        spec: Fixed engine data.

    Returns:
        ``{row: [(fraction, norm, asymmetry), ...]}``, scaled to the standard
        normal space so the norms are comparable across rows.
    """
    sigma_x = tolerance_half_widths(design, grade, angular) / SIGMA_PER_HALF_WIDTH
    widths = tolerance_half_widths(design, grade, angular)
    outcome: dict[str, list[tuple[float, float, float]]] = {}

    for row in rows:
        centre = _metric_gradient(design, row, samples, spec)
        if centre is None:
            continue
        measured = []
        for fraction in fractions:
            steps = np.maximum(fraction * widths, 1.0e-7)
            columns = []
            for index, name in enumerate(VARIABLE_NAMES):
                moved = replace(design, **{name: getattr(design, name) + float(steps[index])})
                shifted = _metric_gradient(moved, row, samples, spec)
                columns.append(
                    np.zeros_like(centre)
                    if shifted is None
                    else (shifted - centre) / steps[index]
                )
            raw = np.stack(columns, axis=1) * np.outer(sigma_x, sigma_x)
            measured.append(
                (
                    float(fraction),
                    float(np.max(np.abs(raw))),
                    hessian_asymmetry(raw[None]),
                )
            )
        outcome[row] = measured
    return outcome


def format_scaling(measured: dict[str, list[tuple[float, float, float]]]) -> str:
    """Render :func:`differencing_scaling` as a table.

    The last column is the tell: ``norm x step`` constant means the numerator
    is a fixed jump rather than a derivative, so there is no Hessian to find.
    """
    lines = ["Hessian by differencing the analytic gradient", "=" * 62, ""]
    for row, entries in measured.items():
        lines.append(f"  {row}")
        lines.append(f"  {'step':>8}{'norm':>13}{'asymmetry':>12}{'norm x step':>14}")
        for fraction, norm, asymmetry in entries:
            lines.append(
                f"  {fraction:>8.2f}{norm:>13.6f}{asymmetry:>12.2f}{norm * fraction:>14.7f}"
            )
        lines.append("")
    lines.append("  a constant last column means the gradient jumps; there is no Hessian")
    return "\n".join(lines)


def format_second_order(outcome: SecondOrder) -> str:
    """Render a :class:`SecondOrder` as an aligned table."""
    lines = ["second-order reliability", "=" * 62, ""]
    lines.append(
        f"  {'constraint':<15}{'beta':>8}{'first order':>14}{'curvature':>12}{'second':>13}"
    )
    for name, beta in zip(outcome.names, outcome.beta, strict=True):
        mark = "" if outcome.valid[name] else "  *"
        lines.append(
            f"  {name:<15}{beta:>8.2f}{outcome.first_order[name]:>14.4e}"
            f"{outcome.curvature[name]:>12.3f}{outcome.second_order[name]:>13.4e}{mark}"
        )
    lines.append("")
    lines.append(f"  system, first order     {outcome.first_order_system:.4e}")
    lines.append(f"  system, second order    {outcome.system:.4e}")
    lines.append(f"  the correction is       x{outcome.ratio:.2f}")
    lines.append(f"  worst Hessian asymmetry {outcome.asymmetry:.2e}")
    if not all(outcome.valid.values()):
        lines.append("  * Breitung's expansion did not hold; left uncorrected")
    return "\n".join(lines)
