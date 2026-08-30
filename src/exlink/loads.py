"""Quasi-static force analysis of the linkage.

Inertia forces and torques are neglected, as in the report: this is a first
sizing iteration, and the masses are not known until the parts have a shape.
The chain runs from the piston down to the two shafts.

**Piston.** With ``P`` the gas force on the crown and ``theta_e`` the piston-rod
angle, ``C = P / sin(theta_e)`` is the rod load and ``D = P cot(theta_e)`` the
side load carried by the liner.

**Trigonal link.** Force balance at ``A`` (swing rod), ``D`` (crank ``q_2``) and
``E`` (piston rod), plus the moment about ``D``, give the swing-rod load ``A``
and the reaction ``Q``.  The moment equation carries ``sin(theta_a - theta_T)``
in its denominator, so ``A`` blows up exactly at the critical configurations
that condition (4a) is there to keep away -- a second, independent reason to
enforce it.

**Shafts.** The gear tooth load ``T`` follows from the moment balance of the
eccentric shaft, and the crankshaft moment balance then yields the output
torque

.. math:: M_r = q_1 A \\cos(\\theta_a - \\theta_1) + r_1 T \\cos\\alpha .
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec
from .kinematics import Kinematics

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Loads:
    """Internal forces and output torque over one crankshaft revolution."""

    piston_force: FloatArray
    """Gas force ``P`` on the piston crown [N]."""

    rod_force: FloatArray
    """Piston-rod load ``C`` [N]."""

    side_force: FloatArray
    """Piston side load ``D`` reacted by the liner [N]."""

    liner_moment: FloatArray
    """Reaction moment ``M_D`` about the wrist pin [N.mm]."""

    swing_force: FloatArray
    """Swing-rod load ``A`` [N]."""

    trigonal_reaction: FloatArray
    """Reaction ``Q = (Q_x, Q_y)`` at ``D``, shaped ``(n_angles, 2)`` [N]."""

    gear_force: FloatArray
    """Gear tooth load ``T`` [N]."""

    eccentric_reaction: FloatArray
    """Eccentric-shaft bearing reaction ``R_2``, shaped ``(n_angles, 2)`` [N]."""

    crank_reaction: FloatArray
    """Crankshaft bearing reaction ``R_1``, shaped ``(n_angles, 2)`` [N]."""

    torque: FloatArray
    """Output torque ``M_r`` on the crankshaft [N.mm]."""

    @property
    def mean_torque(self) -> float:
        """``M_r,ave``, averaged over one revolution [N.mm]."""
        return float(np.mean(self.torque))

    @property
    def mean_piston_force(self) -> float:
        """``P_ave``, averaged over one revolution [N]."""
        return float(np.mean(self.piston_force))

    @property
    def side_load_ratio(self) -> float:
        """``gamma = max(D) / max(P)``, the report's side-load limit."""
        peak = float(np.max(np.abs(self.piston_force)))
        if peak == 0.0:
            return float("inf")
        return float(np.max(np.abs(self.side_force))) / peak


def solve(
    kinematics: Kinematics,
    piston_force: FloatArray,
    spec: EngineSpec = DEFAULT_SPEC,
) -> Loads:
    """Propagate the gas load down to the crankshaft.

    Args:
        kinematics: Pose of the linkage over one revolution.
        piston_force: Gas force on the piston crown at each crank angle [N].
        spec: Fixed engine data (piston length, skirt, gear pressure angle).

    Returns:
        The internal forces and the output torque.
    """
    design = kinematics.design
    a, c, e = design.a, design.c, design.e  # noqa: F841 - e kept for symmetry
    b, theta_b = design.b, design.theta_b
    q_1, q_2 = design.q_1, design.q_2
    r_1, r_2 = design.r_1, design.r_2
    theta_r = design.theta_r_rad
    alpha = spec.pressure_angle

    theta_1 = kinematics.theta_1
    theta_2 = kinematics.theta_2
    theta_T = kinematics.theta_T
    theta_a = kinematics.theta_a
    theta_e = kinematics.theta_e

    force = np.asarray(piston_force, dtype=float)

    # -- piston ------------------------------------------------------------------
    sin_e, cos_e = np.sin(theta_e), np.cos(theta_e)
    rod_force = force / sin_e
    side_force = force * cos_e / sin_e
    liner_moment = side_force * (0.5 * spec.piston_skirt - spec.piston_length)

    # -- trigonal link: moment about D --------------------------------------------
    # DE = b u(theta_b + theta_T) - c u(theta_T); the bracket below is the moment
    # arm of the rod load about D, and sin(theta_a - theta_T) = sin(T).
    de_x = b * np.cos(theta_T + theta_b) - c * np.cos(theta_T)
    de_y = b * np.sin(theta_T + theta_b) - c * np.sin(theta_T)
    lever = de_x * sin_e - de_y * cos_e
    # NOTE the leading minus sign.  Expanding ``DA ^ F_A + DE ^ F_E = 0`` gives
    # ``-c A sin(theta_a - theta_T) - C (DE ^ u_e)_z = 0``, so the swing-rod load
    # is the *negative* of the expression printed in the report -- a sign slip
    # there.  With the sign below the chain reproduces the virtual-work torque
    # ``-P dlambda/dtheta_1`` to machine precision at every crank angle, which
    # is exactly the identity the report appeals to when it defines ``eta``;
    # with the report's sign it does not.  See ``tests/test_loads.py``.
    swing_force = -rod_force * lever / (c * np.sin(theta_a - theta_T))

    reaction_x = rod_force * cos_e - swing_force * np.cos(theta_a)
    reaction_y = rod_force * sin_e - swing_force * np.sin(theta_a)
    trigonal_reaction = np.stack([reaction_x, reaction_y], axis=-1)

    # -- eccentric shaft ----------------------------------------------------------
    gear_force = -(q_2 * reaction_y * np.sin(theta_2) + q_2 * reaction_x * np.cos(theta_2)) / (
        r_2 * np.cos(alpha)
    )
    eccentric_reaction = np.stack(
        [
            reaction_x - gear_force * np.sin(theta_r + alpha),
            reaction_y + gear_force * np.cos(theta_r + alpha),
        ],
        axis=-1,
    )

    # -- crankshaft ---------------------------------------------------------------
    torque = q_1 * swing_force * np.cos(theta_a - theta_1) + r_1 * gear_force * np.cos(alpha)
    crank_reaction = np.stack(
        [
            swing_force * np.cos(theta_a) + gear_force * np.sin(theta_r + alpha),
            swing_force * np.sin(theta_a) - gear_force * np.cos(theta_r + alpha),
        ],
        axis=-1,
    )

    return Loads(
        piston_force=force,
        rod_force=rod_force,
        side_force=side_force,
        liner_moment=liner_moment,
        swing_force=swing_force,
        trigonal_reaction=trigonal_reaction,
        gear_force=gear_force,
        eccentric_reaction=eccentric_reaction,
        crank_reaction=crank_reaction,
        torque=torque,
    )
