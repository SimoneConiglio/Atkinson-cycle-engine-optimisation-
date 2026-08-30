"""Inertia in the load path, and the invariants that pin it."""

from __future__ import annotations

import numpy as np
import pytest

from exlink import analyse
from exlink.dynamics import (
    MEMBER_NAMES,
    MEMBERS,
    mass_properties,
    rpm_to_rad_per_s,
)
from exlink.dynamics import solve as solve_dynamics
from exlink.materials import DEFAULT_MATERIAL
from exlink.reference import REFINED_DESIGN
from exlink.sizing import piston_mass


@pytest.fixture(scope="module")
def solved():
    return analyse(REFINED_DESIGN, samples=720).require_solved()


@pytest.fixture(scope="module")
def properties(solved):
    _, piston = piston_mass(solved.thermodynamics)
    return mass_properties(
        solved.kinematics,
        dict.fromkeys(MEMBER_NAMES, 10.0),
        DEFAULT_MATERIAL.density,
        piston,
    )


def test_zero_speed_reproduces_the_quasi_static_chain(solved, properties) -> None:
    """The whole point of building the 18x18 system: it must contain the old one.

    With no inertia every rod becomes a two-force member again and the
    simultaneous solve has to return exactly what the report's sequential
    elimination gives -- torque, joint forces, gear load and all.
    """
    loads = solve_dynamics(
        solved.kinematics, solved.thermodynamics.piston_force, properties, speed=0.0
    )
    reference = solved.loads
    torque_scale = np.max(np.abs(reference.torque))
    assert np.max(np.abs(loads.torque - reference.torque)) < 1e-8 * torque_scale

    unit_e = np.stack(
        [np.cos(solved.kinematics.theta_e), np.sin(solved.kinematics.theta_e)], axis=-1
    )
    assert loads.reaction["P"] == pytest.approx(reference.rod_force[:, None] * unit_e, abs=1e-6)
    unit_a = np.stack(
        [np.cos(solved.kinematics.theta_a), np.sin(solved.kinematics.theta_a)], axis=-1
    )
    assert loads.reaction["A"] == pytest.approx(
        reference.swing_force[:, None] * unit_a, abs=1e-6
    )
    assert loads.gear_force == pytest.approx(reference.gear_force, abs=1e-6)


@pytest.mark.parametrize("rpm", [0.0, 750.0, 1500.0, 3000.0])
def test_mean_torque_does_not_depend_on_speed(solved, properties, rpm) -> None:
    """Inertia reshapes the torque curve but cannot move its mean.

    At constant crankshaft speed the mechanism returns to its starting state
    every revolution, so its kinetic energy is unchanged and the inertia forces
    do no net work. Any speed dependence in the mean torque would mean the
    equilibrium assembly is leaking or creating energy.
    """
    quasi_static = solve_dynamics(
        solved.kinematics, solved.thermodynamics.piston_force, properties, speed=0.0
    )
    loads = solve_dynamics(
        solved.kinematics,
        solved.thermodynamics.piston_force,
        properties,
        speed=rpm_to_rad_per_s(rpm),
    )
    assert loads.mean_torque == pytest.approx(quasi_static.mean_torque, rel=1e-6)


def test_inertia_does_change_the_instantaneous_loads(solved, properties) -> None:
    """The mean is invariant; the peaks are not, and that is what sizing sees."""
    at_rest = solve_dynamics(
        solved.kinematics, solved.thermodynamics.piston_force, properties, speed=0.0
    )
    at_speed = solve_dynamics(
        solved.kinematics,
        solved.thermodynamics.piston_force,
        properties,
        speed=rpm_to_rad_per_s(1500.0),
    )
    peak_at_rest = np.max(np.linalg.norm(at_rest.reaction["R1"], axis=1))
    peak_at_speed = np.max(np.linalg.norm(at_speed.reaction["R1"], axis=1))
    assert peak_at_speed > 1.2 * peak_at_rest


@pytest.mark.parametrize("joint", ["A", "R1", "P"])
def test_inertia_loads_scale_exactly_with_the_square_of_speed(
    solved, properties, joint
) -> None:
    """With the gas load removed, every reaction is exactly quadratic in speed.

    The whole load case is then linear in the inertia forces, which carry
    ``Omega^2`` as a common factor, so doubling the speed must multiply every
    reaction by exactly four. Superposing peaks of the *combined* case would not
    work -- the gas and inertia peaks fall at different crank angles -- which is
    why the gas force is switched off here rather than subtracted.
    """
    no_gas = np.zeros_like(solved.thermodynamics.piston_force)

    def peak(rpm: float) -> float:
        loads = solve_dynamics(
            solved.kinematics, no_gas, properties, speed=rpm_to_rad_per_s(rpm)
        )
        return float(np.max(np.linalg.norm(loads.reaction[joint], axis=1)))

    assert peak(6000.0) == pytest.approx(4.0 * peak(3000.0), rel=1e-9)


def test_the_shafts_have_no_angular_acceleration(solved, properties) -> None:
    """Constant crankshaft speed means both shafts do, via the 1:2 gear pair."""
    loads = solve_dynamics(
        solved.kinematics,
        solved.thermodynamics.piston_force,
        properties,
        speed=rpm_to_rad_per_s(1500.0),
    )
    assert loads.body_angular_acceleration["crank_1"] == pytest.approx(0.0)
    assert loads.body_angular_acceleration["crank_2"] == pytest.approx(0.0)


def test_joint_accelerations_match_finite_differences(solved, properties) -> None:
    """Guard against the spectral derivative amplifying noise into the loads."""
    loads = solve_dynamics(
        solved.kinematics,
        solved.thermodynamics.piston_force,
        properties,
        speed=rpm_to_rad_per_s(1000.0),
    )
    kinematics = solved.kinematics
    step = 2.0 * np.pi / kinematics.theta_1.size
    speed = rpm_to_rad_per_s(1000.0)
    for name in ("A", "E", "P"):
        points = getattr(kinematics, name)
        difference = (
            speed**2
            * (np.roll(points, -1, axis=0) - 2.0 * points + np.roll(points, 1, axis=0))
            / step**2
        )
        scale = np.max(np.abs(loads.joint_acceleration[name]))
        assert np.max(np.abs(loads.joint_acceleration[name] - difference)) < 5e-3 * scale


def test_the_crank_pin_traces_a_circle_of_the_crank_radius(solved, properties) -> None:
    """A pure centripetal check on the acceleration field."""
    speed = rpm_to_rad_per_s(1500.0)
    loads = solve_dynamics(
        solved.kinematics, solved.thermodynamics.piston_force, properties, speed=speed
    )
    magnitude = np.linalg.norm(loads.joint_acceleration["Q"], axis=1)
    assert magnitude == pytest.approx(speed**2 * REFINED_DESIGN.q_1, rel=1e-6)


def test_the_shafts_do_not_accelerate(solved, properties) -> None:
    """Both bearing centres are fixed points of the mechanism."""
    loads = solve_dynamics(
        solved.kinematics,
        solved.thermodynamics.piston_force,
        properties,
        speed=rpm_to_rad_per_s(1500.0),
    )
    # Compared against a real acceleration in the mechanism rather than to an
    # absolute floor: differentiating a constant leaves round-off, which the
    # Omega^2 factor then scales up.
    scale = np.max(np.linalg.norm(loads.joint_acceleration["Q"], axis=1))
    for name in ("R1", "R2"):
        assert np.max(np.abs(loads.joint_acceleration[name])) < 1e-9 * scale


def test_mass_properties_scale_with_the_section(solved) -> None:
    """Doubling every diameter must quadruple every member mass."""
    _, piston = piston_mass(solved.thermodynamics)
    thin = mass_properties(
        solved.kinematics, dict.fromkeys(MEMBER_NAMES, 10.0), DEFAULT_MATERIAL.density, piston
    )
    thick = mass_properties(
        solved.kinematics, dict.fromkeys(MEMBER_NAMES, 20.0), DEFAULT_MATERIAL.density, piston
    )
    for name in MEMBER_NAMES:
        assert thick.member_mass[name] == pytest.approx(4.0 * thin.member_mass[name])


def test_every_member_has_a_length(solved) -> None:
    for member in MEMBERS:
        assert abs(getattr(REFINED_DESIGN, member.length_attribute)) > 0.0


def test_the_equilibrium_matrix_is_well_conditioned(solved, properties) -> None:
    """Far from the singularity the 18x18 system should solve cleanly.

    The condition number is the same quantity the report's ``W`` constraint
    protects; a design inside that constraint must not be near-singular.
    """
    loads = solve_dynamics(
        solved.kinematics, solved.thermodynamics.piston_force, properties, speed=0.0
    )
    assert np.isfinite(loads.conditioning)
    assert loads.conditioning < 1e8
