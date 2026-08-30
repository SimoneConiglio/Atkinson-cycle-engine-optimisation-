"""Exact derivatives of the analysis chain with respect to the design vector.

Why this exists
---------------
The feasible set of this problem is a *sliver*.  At the reference design the
two equality constraints leave a band 0.1 mm wide on ``STE`` and 0.1 wide on
``epsilon``, inside a box whose sides are tens of millimetres, and ``W`` and
``gamma`` sit within 0.4 % and 7 % of their bounds.  A derivative-free method
cannot work in that: every trial step of a sensible size lands outside, so
COBYLA returns its starting point unchanged however large its budget.

Everything the report derives is closed form, so the derivatives are too.  This
module propagates them forward through the same chain
:mod:`exlink.kinematics` evaluates, carrying ``d/dX`` alongside each quantity.
That turns the problem over to SLSQP, which handles thin feasible sets by
construction.

Two ideas do the work
---------------------
**Forward-mode chaining.**  Each intermediate carries an array of shape
``(n_angles, 11)`` holding its derivative with respect to the design vector.
Differentiating the closed forms is mechanical; the result is exact, and one
pass produces all eleven components at once.

**The envelope theorem.**  Most metrics are extrema over the crank angle --
``W`` is a maximum of ``|cos T|``, ``STE`` a difference of extrema of
``lambda``, ``H`` and ``B`` extremes of a point cloud.  For a maximum attained
at ``theta*``,

.. math:: \\frac{d}{dX}\\max_\\theta f(X, \\theta) = \\frac{\\partial f}{\\partial X}
          \\Big|_{\\theta^*}

because the term through the moving maximiser carries ``partial f / partial
theta = 0``.  So no derivative of the *location* of the extremum is needed --
only the partial derivative evaluated there.

Angles in the design vector are stored in degrees, so their columns carry the
``pi/180`` factor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec
from .design import VARIABLE_NAMES, Design

FloatArray = NDArray[np.float64]

N_VARIABLES = len(VARIABLE_NAMES)
INDEX = {name: i for i, name in enumerate(VARIABLE_NAMES)}
DEGREES = np.pi / 180.0
"""Chain factor for the two design variables held in degrees."""


def _seed(n_angles: int, name: str) -> FloatArray:
    """A ``(n_angles, 11)`` array whose only non-zero column is ``name``."""
    seed = np.zeros((n_angles, N_VARIABLES))
    seed[:, INDEX[name]] = 1.0
    return seed


@dataclass(frozen=True)
class KinematicJacobian:
    """Derivatives of the kinematic solution with respect to ``X``.

    Every array is shaped ``(n_angles, 11)``: one row per crank angle, one
    column per design variable, in the order of
    :data:`exlink.design.VARIABLE_NAMES`.
    """

    theta_1: FloatArray
    cos_transmission: FloatArray
    """``d(cos T)/dX`` -- the argument of condition (4a)."""

    cos_rod: FloatArray
    """``d(cos theta_e)/dX`` -- the argument of condition (6a)."""

    theta_T: FloatArray
    theta_a: FloatArray
    theta_e: FloatArray
    lam: FloatArray
    """``d lambda / dX``, the piston height."""

    joints: dict[str, FloatArray]
    """``d r_J / dX`` for each joint, shaped ``(n_angles, 2, 11)``."""


def kinematic_jacobian(
    design: Design, kinematics: object, spec: EngineSpec = DEFAULT_SPEC
) -> KinematicJacobian:
    """Differentiate the closed-form kinematics, exactly.

    Follows :func:`exlink.kinematics.solve` step for step, carrying the
    derivative of each intermediate alongside it.

    Args:
        design: The mechanism dimensions the solution was computed at.
        kinematics: The corresponding :class:`~exlink.kinematics.Kinematics`.
        spec: Fixed engine data.

    Returns:
        The derivatives of every kinematic quantity with respect to ``X``.
    """
    k = kinematics
    theta_1 = k.theta_1  # type: ignore[attr-defined]
    n = theta_1.size
    zero = np.zeros((n, N_VARIABLES))

    a, c, e = design.a, design.c, design.e
    q_1, q_2, inter = design.q_1, design.q_2, design.I
    b, theta_b = design.b, design.theta_b
    theta_r = design.theta_r_rad

    d_a, d_c, d_e = _seed(n, "a"), _seed(n, "c"), _seed(n, "e")
    d_q1, d_q2 = _seed(n, "q_1"), _seed(n, "q_2")
    d_I, d_x1 = _seed(n, "I"), _seed(n, "x_1")
    d_theta_r = _seed(n, "theta_r") * DEGREES
    d_theta_f = _seed(n, "theta_f") * DEGREES

    # b = hypot(x_b, y_b) and theta_b = atan2(y_b, x_b): the reparametrisation
    # that lets the trigonal link be described without a triangle inequality.
    d_xb, d_yb = _seed(n, "x_b"), _seed(n, "y_b")
    d_b = (design.x_b * d_xb + design.y_b * d_yb) / b
    d_theta_b = (design.x_b * d_yb - design.y_b * d_xb) / b**2

    # (1) the gear relation.
    theta_2 = k.theta_2  # type: ignore[attr-defined]
    d_theta_2 = d_theta_f

    # (3a) the four-bar closure.
    sin_1, cos_1 = np.sin(theta_1)[:, None], np.cos(theta_1)[:, None]
    sin_2, cos_2 = np.sin(theta_2)[:, None], np.cos(theta_2)[:, None]
    proj_a = q_1 * np.sin(theta_1) - q_2 * np.sin(theta_2) + inter * np.cos(theta_r)
    proj_b = -q_1 * np.cos(theta_1) + q_2 * np.cos(theta_2) + inter * np.sin(theta_r)
    d_proj_a = (
        d_q1 * sin_1
        - d_q2 * sin_2
        - q_2 * cos_2 * d_theta_2
        + d_I * np.cos(theta_r)
        - inter * np.sin(theta_r) * d_theta_r
    )
    d_proj_b = (
        -d_q1 * cos_1
        + d_q2 * cos_2
        - q_2 * sin_2 * d_theta_2
        + d_I * np.sin(theta_r)
        + inter * np.cos(theta_r) * d_theta_r
    )

    # (4) the transmission angle.
    square = (proj_a**2 + proj_b**2)[:, None]
    d_square = 2.0 * proj_a[:, None] * d_proj_a + 2.0 * proj_b[:, None] * d_proj_b
    cos_T = ((proj_a**2 + proj_b**2 - a**2 - c**2) / (2.0 * a * c))[:, None]
    d_cos_T = (d_square - 2.0 * a * d_a - 2.0 * c * d_c) / (2.0 * a * c) - cos_T * (
        d_a / a + d_c / c
    )
    transmission = k.transmission_angle[:, None]  # type: ignore[attr-defined]
    sin_T = np.sin(transmission)
    # d(arccos u) = -du / sqrt(1 - u^2); guarded because a design at the
    # singularity has no derivative there either.
    safe = np.maximum(1.0 - cos_T**2, 1.0e-14)
    d_T = -d_cos_T / np.sqrt(safe)

    # (2a) the auxiliary angle q.
    upper = a * sin_T
    lower = a * np.cos(transmission) + c
    d_upper = d_a * sin_T + a * np.cos(transmission) * d_T
    d_lower = d_a * np.cos(transmission) - a * sin_T * d_T + d_c
    d_q = (lower * d_upper - upper * d_lower) / (upper**2 + lower**2)

    # (5) the trigonal and swing-rod orientations.
    d_atan = (proj_a[:, None] * d_proj_b - proj_b[:, None] * d_proj_a) / square
    d_theta_T = d_atan - d_q
    d_theta_a = d_theta_T + d_T

    # (6) the piston-rod angle.
    theta_a_v = k.theta_a[:, None]  # type: ignore[attr-defined]
    theta_T_v = k.theta_T[:, None]  # type: ignore[attr-defined]
    combined = theta_b + theta_T_v
    d_combined = d_theta_b + d_theta_T
    numerator = (
        q_1 * np.sin(theta_1)[:, None]
        - a * np.cos(theta_a_v)
        - b * np.cos(combined)
        + design.x_1
    )
    d_numerator = (
        d_q1 * sin_1
        - d_a * np.cos(theta_a_v)
        + a * np.sin(theta_a_v) * d_theta_a
        - d_b * np.cos(combined)
        + b * np.sin(combined) * d_combined
        + d_x1
    )
    cos_e = numerator / e
    d_cos_e = d_numerator / e - numerator * d_e / e**2
    theta_e_v = k.theta_e[:, None]  # type: ignore[attr-defined]
    safe_e = np.maximum(1.0 - cos_e**2, 1.0e-14)
    d_theta_e = -d_cos_e / np.sqrt(safe_e)

    # (7) the piston height.
    d_lam = (
        d_q1 * cos_1
        + d_a * np.sin(theta_a_v)
        + a * np.cos(theta_a_v) * d_theta_a
        + d_b * np.sin(combined)
        + b * np.cos(combined) * d_combined
        + d_e * np.sin(theta_e_v)
        + e * np.cos(theta_e_v) * d_theta_e
    )

    # -- joint positions ---------------------------------------------------------
    joints: dict[str, FloatArray] = {}
    joints["R1"] = np.zeros((n, 2, N_VARIABLES))
    joints["R2"] = np.stack(
        [
            d_I * np.cos(theta_r) - inter * np.sin(theta_r) * d_theta_r,
            d_I * np.sin(theta_r) + inter * np.cos(theta_r) * d_theta_r,
        ],
        axis=1,
    )
    joints["Q"] = np.stack([-d_q1 * sin_1, d_q1 * cos_1], axis=1)
    d_A = joints["Q"] + np.stack(
        [
            d_a * np.cos(theta_a_v) - a * np.sin(theta_a_v) * d_theta_a,
            d_a * np.sin(theta_a_v) + a * np.cos(theta_a_v) * d_theta_a,
        ],
        axis=1,
    )
    joints["A"] = d_A
    joints["D"] = d_A + np.stack(
        [
            d_c * np.cos(theta_T_v) - c * np.sin(theta_T_v) * d_theta_T,
            d_c * np.sin(theta_T_v) + c * np.cos(theta_T_v) * d_theta_T,
        ],
        axis=1,
    )
    joints["E"] = d_A + np.stack(
        [
            d_b * np.cos(combined) - b * np.sin(combined) * d_combined,
            d_b * np.sin(combined) + b * np.cos(combined) * d_combined,
        ],
        axis=1,
    )
    joints["P"] = np.stack([d_x1, d_lam - zero], axis=1)
    joints["P"][:, 1, :] = d_lam
    joints["H"] = np.stack([d_x1, d_lam], axis=1)
    del spec

    return KinematicJacobian(
        theta_1=np.zeros((n, N_VARIABLES)),
        cos_transmission=d_cos_T,
        cos_rod=d_cos_e,
        theta_T=d_theta_T,
        theta_a=d_theta_a,
        theta_e=d_theta_e,
        lam=d_lam,
        joints=joints,
    )


def _refined_extremum(values: FloatArray, index: int) -> tuple[float, float]:
    """Sub-sample position and value of an extremum, by parabolic fit.

    Returns ``(offset, value)`` where ``offset`` is in samples from ``index``.
    """
    n = values.size
    y0, y1, y2 = values[(index - 1) % n], values[index], values[(index + 1) % n]
    denominator = y0 - 2.0 * y1 + y2
    if denominator == 0.0:
        return 0.0, float(y1)
    offset = 0.5 * (y0 - y2) / denominator
    if abs(offset) > 1.0:
        return 0.0, float(y1)
    return float(offset), float(y1 - 0.25 * (y0 - y2) * offset)


def _at_extremum(derivative: FloatArray, index: int, offset: float) -> FloatArray:
    """Interpolate a derivative row to the refined location of an extremum.

    The envelope theorem needs the partial derivative *at* the extremum, and the
    extremum does not sit on a grid point.  Evaluating at the nearest sample
    instead leaves an error of order the grid spacing -- far too coarse against a
    0.05 mm band on ``STE`` -- so a three-point Lagrange quadratic through the
    samples either side is used, matching the order of the parabolic refinement
    that located the extremum in the first place.
    """
    n = derivative.shape[0]
    before = derivative[(index - 1) % n]
    here = derivative[index]
    after = derivative[(index + 1) % n]
    linear = 0.5 * (after - before)
    curvature = before - 2.0 * here + after
    return here + offset * linear + 0.5 * offset**2 * curvature


#: Outputs whose derivative this module supplies exactly.
ANALYTIC_OUTPUTS: tuple[str, ...] = (
    "expansion_stroke",
    "compression_ratio",
    "stroke_error",
    "compression_ratio_error",
    "tdc_gap_margin",
    "compatibility",
    "compatibility_margin",
    "rod_angle",
    "rod_angle_margin",
    "side_load_ratio",
    "side_load_margin",
)
"""Every constraint that is *tight* at a typical design.

These are the ones a derivative-free method cannot navigate and the ones whose
finite-difference gradients are worst, because each is an extremum over the
crank angle: the grid point attaining it switches as the design moves, so a
difference quotient straddles the switch and reports nonsense. The envelope
theorem gives them exactly instead.

``efficiency``, ``height``, ``width`` and ``clearance`` are left to finite
differences: they are smooth, none of them is tight, and the first carries
moving-boundary terms from the combustion pressure jump that would need the
crank angle of top dead centre differentiated as well.
"""


def metric_jacobian(
    design: Design,
    analysis: object,
    jacobian: KinematicJacobian,
    spec: EngineSpec = DEFAULT_SPEC,
) -> dict[str, FloatArray]:
    """Exact derivatives of the tight metrics, by the envelope theorem.

    Args:
        design: The design the analysis was computed at.
        analysis: A valid :class:`~exlink.model.Analysis`.
        jacobian: The kinematic derivatives from :func:`kinematic_jacobian`.
        spec: Fixed engine data.

    Returns:
        ``{output name: d(output)/dX}``, each shaped ``(11,)``, covering
        :data:`ANALYTIC_OUTPUTS`.
    """
    solved = analysis.require_solved()  # type: ignore[attr-defined]
    k, thermo = solved.kinematics, solved.thermodynamics
    phases = thermo.phases
    area = spec.piston_area

    # -- W: max |cos T| against max |cos theta_e| ---------------------------------
    cos_transmission = np.cos(k.transmission_angle)
    cos_rod = np.cos(k.theta_e)
    i_1 = int(np.argmax(np.abs(cos_transmission)))
    i_2 = int(np.argmax(np.abs(cos_rod)))
    delta_1 = abs(float(cos_transmission[i_1]))
    delta_2 = abs(float(cos_rod[i_2]))
    if delta_1 >= delta_2:
        d_compatibility = np.sign(cos_transmission[i_1]) * jacobian.cos_transmission[i_1]
    else:
        d_compatibility = np.sign(cos_rod[i_2]) * jacobian.cos_rod[i_2]

    # -- mra: max deviation of the piston rod from vertical -----------------------
    deviation = k.theta_e - np.pi / 2.0
    i_rod = int(np.argmax(np.abs(deviation)))
    d_rod_angle = np.degrees(1.0) * np.sign(deviation[i_rod]) * jacobian.theta_e[i_rod]

    # -- strokes: differences of extrema of lambda --------------------------------
    top_1, top_2 = phases.maxima_indices
    deep, shallow = phases.minima_indices
    offsets = {i: _refined_extremum(k.lam, i)[0] for i in (top_1, top_2, deep, shallow)}
    values = {i: _refined_extremum(k.lam, i)[1] for i in (top_1, top_2, deep, shallow)}
    d_lam_at = {i: _at_extremum(jacobian.lam, i, offsets[i]) for i in offsets}

    top = top_1 if values[top_1] >= values[top_2] else top_2
    d_expansion = d_lam_at[top] - d_lam_at[deep]
    d_compression = d_lam_at[top] - d_lam_at[shallow]
    d_ratio = area / spec.dead_volume * d_compression
    gap_sign = np.sign(values[top_1] - values[top_2])
    d_gap = gap_sign * (d_lam_at[top_1] - d_lam_at[top_2])

    # -- gamma: the piston side-load ratio ----------------------------------------
    # The gauge pressure depends on the design only through lambda, its value at
    # top dead centre, and the compression stroke; differentiating the two
    # adiabats gives the rest.
    gamma_exponent = spec.heat_capacity_ratio
    volume = thermo.volume
    d_volume = area * (d_lam_at[top][None, :] - jacobian.lam)
    v_1 = spec.dead_volume + area * phases.compression_stroke
    d_v1 = area * d_compression
    d_epsilon = d_v1 / spec.dead_volume

    from .cycle import Phase

    labels = phases.labels
    d_pressure = np.zeros_like(jacobian.lam)
    is_compression = labels == int(Phase.COMPRESSION)
    if np.any(is_compression):
        ratio = v_1 / volume[is_compression]
        d_pressure[is_compression] = (
            spec.p_intake
            * gamma_exponent
            * ratio[:, None] ** (gamma_exponent - 1.0)
            * (
                d_v1[None, :] / volume[is_compression][:, None]
                - v_1 * d_volume[is_compression] / volume[is_compression][:, None] ** 2
            )
        )
    is_expansion = labels == int(Phase.EXPANSION)
    if np.any(is_expansion):
        epsilon = thermo.compression_ratio
        d_p3 = (
            spec.explosion_ratio
            * spec.p_intake
            * gamma_exponent
            * epsilon ** (gamma_exponent - 1.0)
            * d_epsilon
        )
        ratio = spec.dead_volume / volume[is_expansion]
        d_pressure[is_expansion] = ratio[:, None] ** gamma_exponent * d_p3[None, :] + (
            thermo.p_combustion
            * gamma_exponent
            * ratio[:, None] ** (gamma_exponent - 1.0)
            * (-spec.dead_volume * d_volume[is_expansion] / volume[is_expansion][:, None] ** 2)
        )

    d_force = area * d_pressure
    force = thermo.piston_force
    cot = np.cos(k.theta_e) / np.sin(k.theta_e)
    side = force * cot
    d_side = (
        d_force * cot[:, None] - (force / np.sin(k.theta_e) ** 2)[:, None] * jacobian.theta_e
    )

    i_side = int(np.argmax(np.abs(side)))
    i_force = int(np.argmax(np.abs(force)))
    peak_side = abs(float(side[i_side]))
    peak_force = abs(float(force[i_force]))
    d_peak_side = np.sign(side[i_side]) * d_side[i_side]
    d_peak_force = np.sign(force[i_force]) * d_force[i_force]
    d_side_ratio = d_peak_side / peak_force - peak_side * d_peak_force / peak_force**2

    return {
        "expansion_stroke": d_expansion,
        "stroke_error": d_expansion,
        "compression_ratio": d_ratio,
        "compression_ratio_error": d_ratio,
        "tdc_gap_margin": d_gap,
        "compatibility": d_compatibility,
        "compatibility_margin": d_compatibility,
        "rod_angle": d_rod_angle,
        "rod_angle_margin": d_rod_angle,
        "side_load_ratio": d_side_ratio,
        "side_load_margin": d_side_ratio,
    }
