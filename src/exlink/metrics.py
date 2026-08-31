"""Objectives and constraint measures of the optimization problem.

Collects, from a solved mechanism, the quantities the optimizer works with:

============  ==========================================================
symbol        meaning
============  ==========================================================
``eta``       average mechanical efficiency -- **maximised**
``H``         envelope along the stroke -- **minimised**
``B``         envelope across the stroke -- **minimised**
``STE``       expansion stroke, must equal 74 mm
``epsilon``   compression ratio, must equal 16
``mra``       peak piston-rod tilt, at most 10 deg
``W``         ``max(delta_c1, delta_c2)``, at most 0.985
``g``         gap between the two top dead centres, at most 0.01 mm
``d``         trigonal-link to cylinder clearance, at least 10 mm
``gamma``     ``max(D) / max(P)`` side-load ratio, at most 0.02
============  ==========================================================
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec
from .kinematics import Kinematics
from .loads import Loads

FloatArray = NDArray[np.float64]

EDGE_SAMPLES = 33
"""Points per trigonal-link edge used for the clearance distance."""


def rod_angle_deviation(kinematics: Kinematics) -> float:
    """``mra``: peak deviation of the piston rod from vertical [deg].

    The brief caps this at 10 deg to keep the rod clear of the liner and to
    limit the side load and its friction losses.
    """
    return float(np.degrees(np.max(np.abs(kinematics.theta_e - np.pi / 2.0))))


def envelope(kinematics: Kinematics, spec: EngineSpec = DEFAULT_SPEC) -> tuple[float, float]:
    """``(H, B)``: the bounding box swept by the whole mechanism [mm].

    ``H`` is measured along the stroke (``y``) and ``B`` across it (``x``), as
    they are defined here.  The box covers every configuration of every body:
    the two gear primitives, all moving joints, and the piston envelope over
    its full travel.

    Args:
        kinematics: A solved mechanism.
        spec: Fixed engine data (bore and piston length).

    Returns:
        ``(H, B)`` in mm.
    """
    design = kinematics.design
    points = [
        kinematics.Q,
        kinematics.A,
        kinematics.D,
        kinematics.E,
        kinematics.P,
    ]
    cloud = np.concatenate(points, axis=0)
    x_min, y_min = cloud.min(axis=0)
    x_max, y_max = cloud.max(axis=0)

    # The two gear primitives are discs, not points.
    r_2 = design.r_2
    centre_2 = kinematics.R2[0]
    x_min = min(x_min, -design.r_1, centre_2[0] - r_2)
    x_max = max(x_max, design.r_1, centre_2[0] + r_2)
    y_min = min(y_min, -design.r_1, centre_2[1] - r_2)
    y_max = max(y_max, design.r_1, centre_2[1] + r_2)

    # The piston is a cylinder of diameter ``bore`` sweeping the whole stroke.
    half_bore = 0.5 * spec.bore
    x_min = min(x_min, design.x_1 - half_bore)
    x_max = max(x_max, design.x_1 + half_bore)
    y_min = min(y_min, float(kinematics.lam.min()) - spec.piston_length)
    y_max = max(y_max, float(kinematics.lam.max()))

    return float(y_max - y_min), float(x_max - x_min)


def cylinder_clearance(kinematics: Kinematics, spec: EngineSpec = DEFAULT_SPEC) -> float:
    """``d``: smallest distance from the trigonal link to the cylinder [mm].

    The liner is modelled as the half-strip
    ``x in [x_1 - phi/2, x_1 + phi/2]``, ``y >= y_bottom``, where ``y_bottom``
    is the lowest point the piston skirt reaches over the revolution -- the
    liner has to extend at least that far down.  The distance is minimised over
    the crank revolution and over the three edges ``AD``, ``DE``, ``EA``.

    The requirement is ``d >= 10 mm``.  The construction below is one specific
    reading of it: monotone in the right direction and vanishing exactly on
    contact, which is what the constraint needs.  Any other reasonable
    construction will differ in absolute value while keeping the same feasible
    region shape.
    """
    design = kinematics.design
    half_bore = 0.5 * spec.bore
    left = design.x_1 - half_bore
    right = design.x_1 + half_bore
    bottom = float(kinematics.lam.min()) - spec.piston_length

    weights = np.linspace(0.0, 1.0, EDGE_SAMPLES)[:, None, None]
    edges = [
        (kinematics.A, kinematics.D),
        (kinematics.D, kinematics.E),
        (kinematics.E, kinematics.A),
    ]
    best = np.inf
    for start, end in edges:
        # (EDGE_SAMPLES, n_angles, 2) points along the edge at every crank angle.
        samples = start[None] + weights * (end - start)[None]
        x, y = samples[..., 0], samples[..., 1]
        dx = np.maximum(np.maximum(left - x, x - right), 0.0)
        dy = np.maximum(bottom - y, 0.0)
        best = min(best, float(np.min(np.hypot(dx, dy))))
    return best


@dataclass(frozen=True)
class Metrics:
    """Every objective and constraint measure of one design."""

    efficiency: float
    """``eta``, the average mechanical efficiency [-]."""

    torque_pressure_ratio: float
    """``phi = M_r,ave / p_ave`` [mm^3]."""

    lever_arm: float
    """``M_r,ave / P_ave``, the effective crank lever arm [mm]."""

    height: float
    """``H``, envelope along the stroke [mm]."""

    width: float
    """``B``, envelope across the stroke [mm]."""

    expansion_stroke: float
    """``STE`` [mm]."""

    compression_stroke: float
    """``STC`` [mm]."""

    compression_ratio: float
    """``epsilon`` [-]."""

    rod_angle: float
    """``mra`` [deg]."""

    compatibility: float
    """``W = max(delta_c1, delta_c2)`` [-]."""

    tdc_gap: float
    """``g`` [mm]."""

    clearance: float
    """``d`` [mm]."""

    side_load_ratio: float
    """``gamma = max(D) / max(P)`` [-]."""

    mean_torque: float
    """``M_r,ave`` [N.mm]."""

    mean_piston_force: float
    """``P_ave`` [N]."""

    valid: bool = True
    """``False`` when the design was penalised instead of analysed."""

    reason: str = ""
    """Why the design was penalised, empty when :attr:`valid`."""


def efficiency(
    loads: Loads,
    expansion_stroke: float,
    compression_stroke: float,
    spec: EngineSpec = DEFAULT_SPEC,
) -> tuple[float, float, float]:
    """Return ``(eta, phi, lever_arm)``.

    The average mechanical efficiency is defined as

    .. math::
        \\eta = \\frac{\\int_0^{2\\pi} M_r \\, d\\theta_1}
                     {2 (STE + STC) \\frac{1}{2\\pi}\\int_0^{2\\pi} P \\, d\\theta_1}
             = \\frac{M_{r,ave}}{P_{ave}} \\frac{\\pi}{STE + STC},

    a ratio of two works: the torque's work on the crankshaft over the gas
    force's work on the piston.  Because it rewards a long lever arm, it grows
    without bound as the mechanism grows -- which is exactly why ``H`` and ``B``
    have to enter the problem as well.
    """
    mean_force = loads.mean_piston_force
    if mean_force <= 0.0:
        return 0.0, 0.0, 0.0
    lever_arm = loads.mean_torque / mean_force
    eta = lever_arm * np.pi / (expansion_stroke + compression_stroke)
    return float(eta), float(lever_arm * spec.piston_area), float(lever_arm)
