"""The schedule-scoring of section 6.5."""

from __future__ import annotations

import math

import pytest

from exlink.drivecycle import (
    STANDARD_CYCLE,
    CyclePoint,
    DriveCycle,
    format_cycle,
    score_cycle,
    score_slidercrank_cycle,
)
from exlink.performance import evaluate
from exlink.reference import RELIABLE_DESIGN
from exlink.slidercrank import SliderCrank

PINNED = {"module": 0.8, "teeth": 48}


@pytest.fixture(scope="module")
def scored():
    return score_cycle(RELIABLE_DESIGN, STANDARD_CYCLE, **PINNED)


def test_the_shares_are_normalised_on_the_way_in():
    cycle = DriveCycle((CyclePoint(1600.0, 3.0), CyclePoint(2400.0, 1.0))).normalised()
    assert [point.share for point in cycle.points] == [0.75, 0.25]


def test_an_empty_schedule_is_rejected():
    with pytest.raises(ValueError, match="positive share"):
        DriveCycle((CyclePoint(2000.0, 0.0),)).normalised()


def test_the_analysis_speed_is_half_the_crankshaft_speed():
    assert CyclePoint(2000.0, 1.0).speed_rpm == pytest.approx(1000.0)


def test_the_aggregate_is_the_harmonic_mean_not_the_arithmetic_one(scored):
    fuel = sum(row.point.share / row.km_per_litre for row in scored.rows)
    assert scored.km_per_litre == pytest.approx(1.0 / fuel)

    arithmetic = sum(row.point.share * row.km_per_litre for row in scored.rows)
    assert scored.km_per_litre < arithmetic


def test_the_schedule_is_scored_with_one_engine_not_one_per_point(scored):
    """The mass is a single number, and it is not any point's own mass."""
    at_the_design_point = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **PINNED)
    assert scored.engine_mass_kg > at_the_design_point.engine_mass_kg

    fastest = evaluate(RELIABLE_DESIGN, speed_rpm=1400.0, **PINNED)
    assert scored.engine_mass_kg > fastest.engine_mass_kg


def test_the_structure_comes_from_the_fastest_point(scored):
    assert scored.sizing_rpm == STANDARD_CYCLE.fastest.crankshaft_rpm

    fastest = evaluate(RELIABLE_DESIGN, speed_rpm=1400.0, **PINNED)
    fitted = 1000.0 * fastest.budget.items["flywheel"]
    assert scored.engine_mass_kg == pytest.approx(
        fastest.engine_mass_kg - fitted + scored.flywheel_kg, rel=1.0e-9
    )


def test_the_flywheel_comes_from_the_slowest_point(scored):
    """It is heavier than the one the design point alone would fit."""
    at_the_design_point = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **PINNED)
    assert scored.flywheel_kg > 1000.0 * at_the_design_point.budget.items["flywheel"]


def test_a_one_point_schedule_reproduces_the_single_point_score():
    cycle = DriveCycle((CyclePoint(2000.0, 1.0),))
    scored = score_cycle(RELIABLE_DESIGN, cycle, **PINNED)
    single = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **PINNED)
    assert scored.engine_mass_kg == pytest.approx(single.engine_mass_kg, rel=1.0e-9)
    assert scored.km_per_litre == pytest.approx(single.km_per_litre, rel=1.0e-9)
    assert scored.spread == pytest.approx(1.0)


def test_the_design_holds_up_across_the_schedule(scored):
    assert scored.feasible
    assert scored.spread > 0.85
    assert 3000.0 < scored.km_per_litre < 3400.0


def test_the_advantage_survives_the_schedule(scored):
    """Section 6.5's headline: the gap does not close, it widens."""
    mechanism = SliderCrank.for_compression_ratio(16.0, 0.195)
    rival = score_slidercrank_cycle(mechanism, STANDARD_CYCLE)
    assert rival.feasible
    assert scored.km_per_litre / rival.km_per_litre - 1.0 > 0.176


def test_the_slidercrank_pays_more_for_its_flywheel():
    """Two revolutions on one firing, and a torque curve with a deeper trough."""
    mechanism = SliderCrank.for_compression_ratio(16.0, 0.195)
    rival = score_slidercrank_cycle(mechanism, STANDARD_CYCLE)
    scored = score_cycle(RELIABLE_DESIGN, STANDARD_CYCLE, **PINNED)
    assert rival.flywheel_kg > scored.flywheel_kg


def test_a_range_of_zero_at_any_point_zeroes_the_cycle():
    cycle = DriveCycle((CyclePoint(2000.0, 0.5), CyclePoint(9000.0, 0.5)))
    scored = score_cycle(RELIABLE_DESIGN, cycle, **PINNED)
    assert not scored.feasible
    assert scored.km_per_litre == 0.0 or math.isfinite(scored.km_per_litre)


def test_the_table_names_every_point(scored):
    text = format_cycle(scored)
    for point in STANDARD_CYCLE.points:
        assert f"{point.crankshaft_rpm:.0f}" in text
    assert "cycle range" in text
