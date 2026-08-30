"""Internal loads and sizing against yield, fatigue and buckling."""

from __future__ import annotations

import numpy as np
import pytest

from exlink import analyse
from exlink.dynamics import MEMBER_NAMES, mass_properties
from exlink.dynamics import solve as solve_dynamics
from exlink.materials import (
    DEFAULT_MATERIAL,
    Material,
    goodman_utilisation,
)
from exlink.reference import REFINED_DESIGN
from exlink.sizing import (
    MAX_DIAMETER,
    MIN_DIAMETER,
    internal_loads,
    member_loads,
    piston_mass,
    size_members,
)

ZERO = np.zeros((1, 2))
START = np.array([[0.0, 0.0]])
END = np.array([[100.0, 0.0]])


# -- the closed-form internal loads ------------------------------------------------


def test_a_massless_member_under_an_axial_load_carries_no_moment() -> None:
    axial, bending = internal_loads(START, END, np.array([[-500.0, 0.0]]), 0.0, ZERO, ZERO)
    assert axial == pytest.approx(500.0)
    assert bending == pytest.approx(0.0)


def test_a_massless_cantilever_gives_the_textbook_root_moment() -> None:
    """Tip load ``F`` on a cantilever of length ``L`` gives ``F L`` at the root."""
    axial, bending = internal_loads(START, END, np.array([[0.0, 500.0]]), 0.0, ZERO, ZERO)
    assert bending[0, -1] == pytest.approx(500.0 * 100.0)
    assert axial == pytest.approx(0.0)


def test_uniform_self_inertia_gives_the_textbook_moment() -> None:
    """A uniform bar under transverse acceleration: ``M = m a L / 2`` at the far end."""
    mass, acceleration = 1.0e-6, 1.0e6  # 1 N of total inertia force
    _, bending = internal_loads(
        START,
        END,
        ZERO,
        mass,
        np.array([[0.0, acceleration]]),
        np.array([[0.0, acceleration]]),
    )
    assert abs(bending[0, -1]) == pytest.approx(mass * acceleration * 100.0 / 2.0)


def test_internal_loads_vanish_at_the_loaded_end() -> None:
    _, bending = internal_loads(START, END, np.array([[0.0, 500.0]]), 0.0, ZERO, ZERO)
    assert bending[0, 0] == pytest.approx(0.0)


# -- the fatigue model -------------------------------------------------------------


def test_goodman_ignores_a_compressive_mean_stress() -> None:
    """The Goodman line is truncated at zero mean, which is the safe convention."""
    tensile = goodman_utilisation(
        np.array([100.0]), np.array([200.0]), np.array([250.0]), 900.0
    )
    compressive = goodman_utilisation(
        np.array([100.0]), np.array([-200.0]), np.array([250.0]), 900.0
    )
    assert tensile > compressive
    assert compressive == pytest.approx(100.0 / 250.0)


def test_the_size_factor_penalises_thick_sections() -> None:
    material = Material()
    assert material.size_factor(10.0) > material.size_factor(80.0)
    assert material.endurance_limit(10.0) > material.endurance_limit(80.0)


def test_the_endurance_limit_is_below_the_ultimate_strength() -> None:
    material = Material()
    assert 0.0 < float(material.endurance_limit(20.0)[0]) < material.ultimate_strength


# -- sizing ------------------------------------------------------------------------


@pytest.fixture(scope="module")
def loads_at_rest():
    solved = analyse(REFINED_DESIGN, samples=360).require_solved()
    _, piston = piston_mass(solved.thermodynamics)
    properties = mass_properties(
        solved.kinematics,
        dict.fromkeys(MEMBER_NAMES, 10.0),
        DEFAULT_MATERIAL.density,
        piston,
    )
    return solve_dynamics(
        solved.kinematics, solved.thermodynamics.piston_force, properties, speed=0.0
    )


def test_every_member_gets_a_buildable_diameter(loads_at_rest) -> None:
    sizing = size_members(loads_at_rest)
    assert set(sizing) == set(MEMBER_NAMES)
    for item in sizing.values():
        assert MIN_DIAMETER <= item.diameter < MAX_DIAMETER
        assert item.mass > 0.0


def test_sizing_drives_the_critical_mode_to_its_allowable(loads_at_rest) -> None:
    """The solve looks for the *smallest* safe section, so one mode must be tight."""
    for item in size_members(loads_at_rest).values():
        worst = max(
            item.static_utilisation,
            item.fatigue_utilisation,
            item.buckling_utilisation,
        )
        assert worst == pytest.approx(1.0, abs=1e-3)


def test_the_reported_critical_mode_is_the_binding_one(loads_at_rest) -> None:
    for item in size_members(loads_at_rest).values():
        utilisations = {
            "static": item.static_utilisation,
            "fatigue": item.fatigue_utilisation,
            "buckling": item.buckling_utilisation,
        }
        assert item.critical_mode == max(utilisations, key=utilisations.__getitem__)


def test_a_bigger_safety_factor_needs_more_material(loads_at_rest) -> None:
    from exlink.materials import SafetyFactors

    lean = size_members(loads_at_rest, safety=SafetyFactors(static=1.5, fatigue=2.0))
    cautious = size_members(loads_at_rest, safety=SafetyFactors(static=3.0, fatigue=4.0))
    assert sum(s.mass for s in cautious.values()) > sum(s.mass for s in lean.values())


def test_a_stronger_material_needs_less_material(loads_at_rest) -> None:
    strong = Material(yield_strength=1400.0, ultimate_strength=1800.0)
    lean = size_members(loads_at_rest, material=strong)
    standard = size_members(loads_at_rest)
    assert sum(s.mass for s in lean.values()) < sum(s.mass for s in standard.values())


def test_member_loads_cover_every_member(loads_at_rest) -> None:
    per_member = member_loads(loads_at_rest)
    assert set(per_member) == set(MEMBER_NAMES)
    for axial, bending in per_member.values():
        assert axial.shape == bending.shape
        assert np.all(np.isfinite(axial))
        assert np.all(np.isfinite(bending))


def test_trigonal_sides_carry_almost_no_bending_at_rest(loads_at_rest) -> None:
    """With no inertia the pin-jointed triangle is a pure truss."""
    per_member = member_loads(loads_at_rest)
    for name in ("trigonal_ad", "trigonal_ae", "trigonal_de"):
        axial, bending = per_member[name]
        assert np.max(np.abs(bending)) < 1e-6 * np.max(np.abs(axial)) * 100.0


def test_the_piston_crown_thickens_with_pressure() -> None:
    solved = analyse(REFINED_DESIGN, samples=180).require_solved()
    crown, mass = piston_mass(solved.thermodynamics)
    assert crown >= 2.0
    assert mass > 0.0
    weak = Material(yield_strength=200.0, ultimate_strength=400.0)
    thicker, heavier = piston_mass(solved.thermodynamics, material=weak)
    assert thicker >= crown
    assert heavier >= mass
