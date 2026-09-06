"""Hollow sections: the four shape factors, and what they do downstream."""

from __future__ import annotations

import math

import numpy as np
import pytest

from exlink.coupled import solve_for_design
from exlink.dynamics import MEMBERS, mass_properties
from exlink.model import analyse
from exlink.performance import evaluate
from exlink.reference import RELIABLE_DESIGN
from exlink.sections import (
    HOLLOW_KINDS,
    MAX_BORE_RATIO,
    SOLID,
    TUBULAR,
    Section,
    area_factor,
    gyration_factor,
    modulus_factor,
)
from exlink.sizing import MEMBER_KINDS

PINNED = {"module": 0.8, "teeth": 48}


def test_a_solid_bar_is_the_zero_of_every_factor():
    assert SOLID.solid
    assert area_factor(np.zeros(3)) == pytest.approx(1.0)
    assert modulus_factor(np.zeros(3)) == pytest.approx(1.0)
    assert gyration_factor(np.zeros(3)) == pytest.approx(1.0)


def test_the_factors_are_the_textbook_ones():
    k = np.array([0.5])
    outer, inner = 20.0, 10.0
    area = math.pi * (outer**2 - inner**2) / 4.0
    second = math.pi * (outer**4 - inner**4) / 64.0
    assert math.pi * outer**2 / 4.0 * area_factor(k)[0] == pytest.approx(area)
    assert math.pi * outer**4 / 64.0 * modulus_factor(k)[0] == pytest.approx(second)


def test_a_bore_takes_mass_faster_than_it_takes_stiffness():
    """The whole reason tubes exist."""
    k = np.array([0.75])
    assert area_factor(k)[0] < modulus_factor(k)[0]


def test_a_tube_carries_more_inertia_per_kilogram_than_a_bar():
    assert gyration_factor(np.array([0.75]))[0] > 1.0


def test_an_impossible_bore_ratio_is_rejected():
    with pytest.raises(ValueError, match="bore ratio"):
        Section(bore_ratio=MAX_BORE_RATIO + 0.01)
    with pytest.raises(ValueError, match="bore ratio"):
        Section(bore_ratio=-0.1)
    with pytest.raises(ValueError, match="wall floor"):
        Section(bore_ratio=0.5, min_wall=0.0)


def test_only_the_kinds_that_can_be_made_from_tube_are_bored():
    """A crank throw is a forging; a rod is tube."""
    ratios = TUBULAR.ratios(MEMBER_KINDS)
    for kind, ratio in zip(MEMBER_KINDS, ratios, strict=True):
        assert (ratio > 0.0) == (kind in HOLLOW_KINDS)
    assert "cantilever" not in HOLLOW_KINDS
    assert {member.kind for member in MEMBERS} >= HOLLOW_KINDS


def test_the_wall_floor_sets_a_smallest_outer_diameter():
    section = Section(bore_ratio=0.6, min_wall=1.5)
    floors = section.minimum_diameter(MEMBER_KINDS)
    hollow = floors > 0.0
    # A wall of exactly min_wall at the floor diameter.
    assert floors[hollow] == pytest.approx(2.0 * 1.5 / (1.0 - 0.6))
    # Solid members have no wall and so no floor.
    assert not hollow.all()
    assert SOLID.minimum_diameter(MEMBER_KINDS) == pytest.approx(0.0)


def test_the_sizing_respects_the_wall_floor():
    solid = solve_for_design(RELIABLE_DESIGN, speed_rpm=1000.0)
    section = Section(bore_ratio=0.6, min_wall=1.5)
    tube = solve_for_design(RELIABLE_DESIGN, speed_rpm=1000.0, section=section)
    floor = 2.0 * 1.5 / 0.4
    for member, kind in zip(MEMBERS, MEMBER_KINDS, strict=True):
        if kind in HOLLOW_KINDS:
            assert tube.diameters[member.name] >= floor - 1.0e-9
        else:
            # An unbored member is sized exactly as before, to within the
            # feedback the lighter tubes send round the fixed point.
            assert tube.diameters[member.name] == pytest.approx(
                solid.diameters[member.name], rel=0.02
            )


def test_a_bore_removes_member_mass_where_the_floor_does_not_bite():
    """A generous wall floor lets the bore do what it is for."""
    solved = analyse(RELIABLE_DESIGN).require_solved()
    diameters = dict.fromkeys([m.name for m in MEMBERS], 20.0)
    solid = mass_properties(solved.kinematics, diameters, 7.85e-9, 1.0e-4)
    tube = mass_properties(
        solved.kinematics,
        diameters,
        7.85e-9,
        1.0e-4,
        section=Section(bore_ratio=0.6, min_wall=0.5),
    )
    # The piston is unbored, so at least one body is untouched.
    lighter = [
        body
        for body in solid.body_mass
        if tube.body_mass[body] < solid.body_mass[body] - 1.0e-12
    ]
    assert lighter
    assert tube.body_mass["piston"] == pytest.approx(solid.body_mass["piston"])


def test_the_study_result_is_unchanged_when_the_sections_stay_solid():
    """Every number in sections 6.1 to 6.5 must survive the new parameter."""
    default = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **PINNED)
    explicit = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, section=SOLID, **PINNED)
    assert explicit.km_per_litre == pytest.approx(default.km_per_litre, rel=1.0e-12)
    assert explicit.engine_mass_kg == pytest.approx(default.engine_mass_kg, rel=1.0e-12)


def test_the_bore_buys_almost_nothing_at_the_design_speed():
    """Section 6.6's first finding: the members are 1.3 % of this engine."""
    solid = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **PINNED)
    assert solid.budget.shares()["linkage"] < 0.02

    tube = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, section=TUBULAR, **PINNED)
    assert abs(tube.km_per_litre / solid.km_per_litre - 1.0) < 0.01


def test_the_bore_buys_a_great_deal_where_inertia_sizes_the_members():
    """Section 6.6's second: it is a speed effect, not a mass effect."""
    speed = 1600.0
    solid = evaluate(RELIABLE_DESIGN, speed_rpm=speed, **PINNED)
    tube = evaluate(RELIABLE_DESIGN, speed_rpm=speed, section=TUBULAR, **PINNED)
    assert tube.km_per_litre > 1.2 * solid.km_per_litre
    assert tube.budget.items["linkage"] < solid.budget.items["linkage"]
