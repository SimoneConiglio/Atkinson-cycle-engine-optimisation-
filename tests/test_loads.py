"""The quasi-static force chain, checked against the principle of virtual work."""

from __future__ import annotations

import numpy as np
import pytest

from exlink import PUBLISHED_DESIGN, analyse
from exlink.reference import REFINED_DESIGN


@pytest.mark.parametrize(
    "design", [PUBLISHED_DESIGN, REFINED_DESIGN], ids=["published", "refined"]
)
def test_torque_equals_the_virtual_work_torque(design) -> None:
    """``M_r`` must equal ``-P dlambda/dtheta_1`` at every crank angle.

    In a massless, frictionless, quasi-static mechanism the instantaneous power
    in equals the power out, so the whole force chain -- piston, trigonal link,
    both shafts, the gear pair -- is pinned by this one identity.  It is what
    justifies calling ``eta`` an efficiency, and it is what exposed the sign
    slip in a careless inversion of the
    trigonal-link moment equation: with the sign as printed, this test fails by
    a factor of about -4.
    """
    analysis = analyse(design, samples=7200)
    kinematics, loads = analysis.kinematics, analysis.loads
    gradient = np.gradient(kinematics.lam, kinematics.theta_1, edge_order=2)
    expected = -analysis.thermodynamics.piston_force * gradient

    scale = np.max(np.abs(expected))
    assert np.max(np.abs(loads.torque - expected)) < 1e-4 * scale


def test_mean_torque_matches_the_indicated_cycle_work(refined_analysis) -> None:
    """``M_r,ave`` must equal the p-V loop area divided by ``2 pi``.

    A route to the same number that never touches the force chain: integrate
    the gauge pressure around the closed p-V loop to get the indicated work per
    cycle. It ties the thermodynamics to the mechanics.
    """
    thermo = refined_analysis.thermodynamics
    # Close the loop: the last sample is followed by the first.
    pressure = np.append(thermo.gauge_pressure, thermo.gauge_pressure[0])
    volume = np.append(thermo.volume, thermo.volume[0])
    work = np.trapezoid(pressure, volume)
    assert refined_analysis.loads.mean_torque == pytest.approx(work / (2.0 * np.pi), rel=2e-3)


def test_the_piston_rod_balances_the_gas_force(refined_analysis) -> None:
    """``C sin(theta_e) = P`` and ``D = P cot(theta_e)``."""
    loads = refined_analysis.loads
    theta_e = refined_analysis.kinematics.theta_e
    assert loads.rod_force * np.sin(theta_e) == pytest.approx(loads.piston_force)
    assert loads.side_force == pytest.approx(loads.piston_force / np.tan(theta_e))


def test_the_trigonal_link_is_in_force_balance(refined_analysis) -> None:
    """The three forces on the trigonal link must sum to zero."""
    loads = refined_analysis.loads
    kinematics = refined_analysis.kinematics
    at_e = -loads.rod_force[:, None] * np.stack(
        [np.cos(kinematics.theta_e), np.sin(kinematics.theta_e)], axis=-1
    )
    at_a = loads.swing_force[:, None] * np.stack(
        [np.cos(kinematics.theta_a), np.sin(kinematics.theta_a)], axis=-1
    )
    residual = at_e + at_a + loads.trigonal_reaction
    assert np.max(np.abs(residual)) < 1e-8 * np.max(np.abs(loads.rod_force))


def test_the_gear_pair_transmits_no_net_power(refined_analysis) -> None:
    """``r_1 omega_1 + r_2 omega_2 = 0``, so the tooth loads cancel in power.

    With ``omega_2 = -2 omega_1`` and ``r_1 = 2 r_2`` the two gear torques must
    cancel exactly; a gear pair that generated power would silently inflate the
    efficiency.
    """
    design = refined_analysis.design
    loads = refined_analysis.loads
    alpha = refined_analysis.spec.pressure_angle
    torque_1 = design.r_1 * loads.gear_force * np.cos(alpha)
    torque_2 = design.r_2 * loads.gear_force * np.cos(alpha)
    power = torque_1 * 1.0 + torque_2 * (-2.0)
    assert np.max(np.abs(power)) < 1e-9 * np.max(np.abs(torque_1))


def test_side_load_ratio_is_max_side_over_max_gas_force(refined_analysis) -> None:
    loads = refined_analysis.loads
    expected = np.max(np.abs(loads.side_force)) / np.max(np.abs(loads.piston_force))
    assert loads.side_load_ratio == pytest.approx(expected)
