"""Exact derivatives through the sizing / dynamics coupling.

Each piece is compared against a converged central difference.  As in
``tests/test_jacobian.py``, the step size matters: several quantities are
extrema over the crank angle, and at a step large enough to move the sample
attaining one, the difference quotient straddles the switch and is wrong.  That
is pinned below too, because it is the reason these are analytic.
"""

from __future__ import annotations

import numpy as np
import pytest

from exlink import Design, analyse
from exlink.design import VARIABLE_NAMES
from exlink.dynamics import (
    MEMBER_NAMES,
    mass_properties,
    rpm_to_rad_per_s,
)
from exlink.dynamics import solve as solve_dynamics
from exlink.dynamics_jacobian import (
    DESIGN_SLICE,
    N_DESIGN,
    N_PARAMETERS,
    acceleration_jacobian,
    coupled_jacobian,
    mass_property_jacobian,
    member_length_jacobian,
    sizing_jacobian,
)
from exlink.jacobian import kinematic_jacobian
from exlink.materials import DEFAULT_MATERIAL, DEFAULT_SAFETY
from exlink.reference import REFINED_DESIGN
from exlink.sizing import (
    STATIONS,
    member_lengths,
    member_loads,
    piston_mass,
    size_from_arrays,
)

SAMPLES = 360
RPM = 1000.0
DIAMETERS = {name: 10.0 + index for index, name in enumerate(MEMBER_NAMES)}


def _evaluate(design: Design, diameters: dict[str, float]):
    """Analysis, mass properties and load case at one point."""
    analysis = analyse(design, samples=SAMPLES)
    _, piston = piston_mass(analysis.thermodynamics)
    properties = mass_properties(
        analysis.kinematics, diameters, DEFAULT_MATERIAL.density, piston
    )
    loads = solve_dynamics(
        analysis.kinematics,
        analysis.thermodynamics.piston_force,
        properties,
        rpm_to_rad_per_s(RPM),
    )
    return analysis, properties, loads


def _perturb(index: int, step: float) -> tuple[Design, dict[str, float]]:
    """Move one entry of the combined ``(X, diameters)`` parameter vector."""
    values = REFINED_DESIGN.to_array()
    diameters = dict(DIAMETERS)
    if index < N_DESIGN:
        values[index] += step
    else:
        diameters[MEMBER_NAMES[index - N_DESIGN]] += step
    return Design.from_array(values), diameters


def _step(index: int) -> float:
    if index < N_DESIGN:
        return 1e-6 * max(abs(REFINED_DESIGN.to_array()[index]), 1.0)
    return 1e-6 * DIAMETERS[MEMBER_NAMES[index - N_DESIGN]]


@pytest.fixture(scope="module")
def state():
    analysis, properties, loads = _evaluate(REFINED_DESIGN, DIAMETERS)
    kinematic = kinematic_jacobian(REFINED_DESIGN, analysis.kinematics)
    mass = mass_property_jacobian(
        REFINED_DESIGN,
        analysis.kinematics,
        kinematic,
        properties,
        DIAMETERS,
        DEFAULT_MATERIAL.density,
    )
    acceleration = acceleration_jacobian(kinematic, mass, rpm_to_rad_per_s(RPM), SAMPLES)
    return analysis, properties, loads, kinematic, mass, acceleration


def test_member_length_derivatives(state) -> None:
    """The two derived trigonal sides carry the reparametrisation's chain rule."""
    rows = member_length_jacobian(REFINED_DESIGN)
    assert rows.shape == (len(MEMBER_NAMES), len(VARIABLE_NAMES))
    for index in range(N_DESIGN):
        step = _step(index)
        forward, _ = _perturb(index, step)
        backward, _ = _perturb(index, -step)
        numerical = (member_lengths(forward) - member_lengths(backward)) / (2.0 * step)
        assert rows[:, index] == pytest.approx(numerical, abs=1e-6)


def test_mass_property_derivatives(state) -> None:
    """Masses and centres of mass, against differences."""
    _, properties, _, _, mass, _ = state
    for index in range(N_PARAMETERS):
        step = _step(index)
        _, forward_properties, _ = _evaluate(*_perturb(index, step))
        _, backward_properties, _ = _evaluate(*_perturb(index, -step))
        for member in MEMBER_NAMES:
            numerical = (
                forward_properties.member_mass[member] - backward_properties.member_mass[member]
            ) / (2.0 * step)
            analytic = mass.member_mass[MEMBER_NAMES.index(member), index]
            assert analytic == pytest.approx(numerical, rel=1e-6, abs=1e-16)
        for body in properties.body_mass:
            numerical = (
                forward_properties.body_com[body] - backward_properties.body_com[body]
            ) / (2.0 * step)
            # Scaled by the centre-of-mass position rather than by its own
            # derivative: many entries are exactly zero, and a relative
            # comparison there measures only the difference's round-off.
            scale = float(np.max(np.abs(properties.body_com[body])))
            assert np.max(np.abs(mass.body_com[body][..., index] - numerical)) < 1e-6 * scale


def test_acceleration_derivatives_reuse_the_spectral_operator(state) -> None:
    """``a = Omega^2 D^2 r`` with ``D^2`` linear, so the derivative is the same map."""
    _, _properties, loads, _, _, acceleration = state
    for index in range(N_PARAMETERS):
        step = _step(index)
        _, _, forward = _evaluate(*_perturb(index, step))
        _, _, backward = _evaluate(*_perturb(index, -step))
        for joint in ("Q", "A", "D", "E", "P"):
            numerical = (
                forward.joint_acceleration[joint] - backward.joint_acceleration[joint]
            ) / (2.0 * step)
            scale = max(float(np.max(np.abs(loads.joint_acceleration[joint]))), 1.0)
            assert (
                np.max(np.abs(acceleration.joint[joint][..., index] - numerical)) < 1e-4 * scale
            )


def test_the_equilibrium_solve_is_differentiated_exactly(state) -> None:
    """``dx = A^-1 (db - dA x)``, against differences of the whole solve."""
    analysis, properties, loads, _, _, _ = state
    derivatives = coupled_jacobian(
        REFINED_DESIGN,
        analysis,
        DIAMETERS,
        properties,
        loads,
        STATIONS,
        DEFAULT_MATERIAL,
    )
    scale = float(np.max(np.abs(loads.torque)))
    for index in range(N_PARAMETERS):
        step = _step(index)
        _, _, forward = _evaluate(*_perturb(index, step))
        _, _, backward = _evaluate(*_perturb(index, -step))
        numerical = (forward.mean_torque - backward.mean_torque) / (2.0 * step)
        assert derivatives.mean_torque[index] == pytest.approx(numerical, abs=1e-5 * scale)


def test_internal_load_derivatives(state) -> None:
    """Axial force and bending moment, at every station of every member."""
    analysis, properties, loads, _, _, _ = state
    derivatives = coupled_jacobian(
        REFINED_DESIGN,
        analysis,
        DIAMETERS,
        properties,
        loads,
        STATIONS,
        DEFAULT_MATERIAL,
    )

    def stacked(case):
        per_member = member_loads(case, stations=STATIONS)
        return (
            np.stack([per_member[n][0] for n in MEMBER_NAMES]),
            np.stack([per_member[n][1] for n in MEMBER_NAMES]),
        )

    axial, bending = stacked(loads)
    axial_scale = float(np.max(np.abs(axial)))
    bending_scale = float(np.max(np.abs(bending)))
    for index in range(N_PARAMETERS):
        step = _step(index)
        _, _, forward = _evaluate(*_perturb(index, step))
        _, _, backward = _evaluate(*_perturb(index, -step))
        forward_axial, forward_bending = stacked(forward)
        backward_axial, backward_bending = stacked(backward)
        assert (
            np.max(
                np.abs(
                    derivatives.axial[..., index]
                    - (forward_axial - backward_axial) / (2.0 * step)
                )
            )
            < 1e-5 * axial_scale
        )
        assert (
            np.max(
                np.abs(
                    derivatives.bending[..., index]
                    - (forward_bending - backward_bending) / (2.0 * step)
                )
            )
            < 1e-5 * bending_scale
        )


def _sized(axial, bending, lengths):
    result = size_from_arrays(axial, bending, lengths, DEFAULT_MATERIAL, DEFAULT_SAFETY)
    return np.array([result[name].diameter for name in MEMBER_NAMES])


def test_the_sizing_bisection_is_differentiated_by_implicit_function(state) -> None:
    """``dd/dq = -(dU/dq)/(dU/dd)`` -- the bisection itself is never differentiated."""
    _, _, loads, _, _, _ = state
    per_member = member_loads(loads, stations=STATIONS)
    axial = np.stack([per_member[n][0] for n in MEMBER_NAMES])
    bending = np.stack([per_member[n][1] for n in MEMBER_NAMES])
    lengths = member_lengths(REFINED_DESIGN)
    diameters = _sized(axial, bending, lengths)
    jacobian = sizing_jacobian(
        axial, bending, diameters, lengths, DEFAULT_MATERIAL, DEFAULT_SAFETY
    )

    generator = np.random.default_rng(0)
    direction_axial = generator.normal(size=axial.shape) * np.max(np.abs(axial))
    direction_bending = generator.normal(size=bending.shape) * np.max(np.abs(bending))
    analytic = np.sum(jacobian.d_axial * direction_axial, axis=(1, 2)) + np.sum(
        jacobian.d_bending * direction_bending, axis=(1, 2)
    )
    # Small enough that the crank angle attaining each extremum does not move.
    scale = 1e-7
    numerical = (
        _sized(axial + scale * direction_axial, bending + scale * direction_bending, lengths)
        - _sized(axial - scale * direction_axial, bending - scale * direction_bending, lengths)
    ) / (2.0 * scale)
    assert analytic == pytest.approx(numerical, rel=1e-3)


def test_a_coarse_load_perturbation_breaks_the_difference_quotient(state) -> None:
    """Why the sizing derivative is analytic rather than differenced.

    Perturb the loads hard enough and the crank angle attaining the fatigue
    extremum hops across a near-flat minimum; the difference quotient then
    reports a derivative tens of percent away from the truth, while the implicit
    function theorem is unaffected.
    """
    _, _, loads, _, _, _ = state
    per_member = member_loads(loads, stations=STATIONS)
    axial = np.stack([per_member[n][0] for n in MEMBER_NAMES])
    bending = np.stack([per_member[n][1] for n in MEMBER_NAMES])
    lengths = member_lengths(REFINED_DESIGN)
    diameters = _sized(axial, bending, lengths)
    jacobian = sizing_jacobian(
        axial, bending, diameters, lengths, DEFAULT_MATERIAL, DEFAULT_SAFETY
    )

    generator = np.random.default_rng(0)
    direction_axial = generator.normal(size=axial.shape) * np.max(np.abs(axial))
    direction_bending = generator.normal(size=bending.shape) * np.max(np.abs(bending))
    analytic = np.sum(jacobian.d_axial * direction_axial, axis=(1, 2)) + np.sum(
        jacobian.d_bending * direction_bending, axis=(1, 2)
    )

    def difference(scale: float) -> np.ndarray:
        return (
            _sized(
                axial + scale * direction_axial,
                bending + scale * direction_bending,
                lengths,
            )
            - _sized(
                axial - scale * direction_axial,
                bending - scale * direction_bending,
                lengths,
            )
        ) / (2.0 * scale)

    assert analytic == pytest.approx(difference(1e-7), rel=1e-3)
    coarse = difference(1e-4)
    assert np.max(np.abs(coarse - analytic) / np.abs(analytic)) > 0.1


def test_the_design_half_of_the_parameter_vector_is_the_design_vector() -> None:
    assert len(VARIABLE_NAMES) == DESIGN_SLICE.stop
    assert len(VARIABLE_NAMES) + len(MEMBER_NAMES) == N_PARAMETERS
