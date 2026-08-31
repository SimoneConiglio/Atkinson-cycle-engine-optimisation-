"""The second mechanism, and what the head-to-head comparison actually shows."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.constants import DEFAULT_SPEC
from exlink.performance import evaluate
from exlink.reference import COUPLED_DESIGN
from exlink.slidercrank import (
    SliderCrank,
    evaluate_slidercrank,
    firing_frequency_sensitivity,
    friction_work,
    kinematics,
    mass_budget,
    otto_cycle,
    solve,
)

SAMPLES = 360


@pytest.fixture(scope="module")
def mechanism() -> SliderCrank:
    return SliderCrank.for_compression_ratio(16.0)


@pytest.fixture(scope="module")
def solved(mechanism: SliderCrank) -> object:
    return solve(mechanism, 1500.0, samples=SAMPLES)


# -- the mechanism ------------------------------------------------------------


def test_compression_ratio_is_realised(mechanism: SliderCrank) -> None:
    cycle = otto_cycle(kinematics(mechanism, SAMPLES)["lam"])
    assert float(cycle["compression_ratio"]) == pytest.approx(16.0, rel=1.0e-9)


def test_piston_motion_is_the_textbook_slider_crank(mechanism: SliderCrank) -> None:
    """``lambda = r cos(theta) + sqrt(l^2 - r^2 sin^2(theta))``, spot-checked.

    At the two dead centres the closed form degenerates to ``l +/- r``, which
    is the cheapest possible check that the kinematics is the right one.
    """
    motion = kinematics(mechanism, SAMPLES)
    assert float(np.max(motion["lam"])) == pytest.approx(mechanism.rod + mechanism.crank)
    assert float(np.min(motion["lam"])) == pytest.approx(
        mechanism.rod - mechanism.crank, rel=1.0e-6
    )


def test_indicated_efficiency_matches_the_ideal_otto_formula(mechanism: SliderCrank) -> None:
    """``eta = 1 - epsilon^(1-gamma)``, exactly.

    The strongest available check on the cycle model: an independent
    closed-form result the numerical p-V loop has to reproduce.  It also
    isolates the one thermodynamic difference from the Atkinson cycle, since
    the compression ratio is identical.
    """
    cycle = otto_cycle(kinematics(mechanism, 1440)["lam"])
    quantity = (
        DEFAULT_SPEC.dead_volume
        * (float(cycle["p_combustion"]) - float(cycle["p_compression_end"]))
        / (DEFAULT_SPEC.heat_capacity_ratio - 1.0)
    )
    ideal = 1.0 - 16.0 ** (1.0 - DEFAULT_SPEC.heat_capacity_ratio)
    assert float(cycle["indicated_work"]) / quantity == pytest.approx(ideal, rel=2.0e-3)


def test_torque_integral_equals_the_pv_loop(solved: object) -> None:
    """Virtual work, on the second mechanism.

    The same identity that validates the EX-link force chain: whatever the
    linkage, the work the gas does on the piston must arrive at the crankshaft.
    Inertia reshapes the torque curve but does no net work over a closed cycle,
    so this holds at speed as well as at rest.
    """
    work = 4.0 * np.pi * solved.mean_torque
    assert work == pytest.approx(solved.indicated_work, rel=1.0e-3)


def test_the_sizing_loop_converges(solved: object) -> None:
    assert solved.converged
    assert all(2.0 < value < 100.0 for value in solved.diameters.values())


def test_a_well_conditioned_mechanism_gets_inertia_relief(mechanism: SliderCrank) -> None:
    """Speed *reduces* the peak main-bearing load here, and that is correct.

    The peak gas force lands near top dead centre, where the reciprocating
    masses are decelerating and their inertia force pulls the opposite way.  A
    well-proportioned slider-crank therefore sees its peak journal load fall as
    it speeds up -- the classic inertia relief of the gas load, and the reason
    high-speed engines do not need proportionally bigger main bearings.

    This is the exact opposite of what the near-singular EX-link does, and the
    contrast is the generalisation the second mechanism exists to establish:
    inertia relieves a well-conditioned linkage and amplifies an ill-conditioned
    one.  Same physics, opposite sign, and conditioning is what decides which.
    """
    at_rest = solve(mechanism, 0.0, samples=SAMPLES)
    at_speed = solve(mechanism, 3000.0, samples=SAMPLES)
    assert at_speed.peak_bearing_load < at_rest.peak_bearing_load
    assert sum(at_speed.member_mass.values()) < sum(at_rest.member_mass.values())


def test_the_near_singular_linkage_gets_the_opposite(mechanism: SliderCrank) -> None:
    """The EX-link cannot survive the speeds the slider-crank prefers.

    The slider-crank is still comfortable at 3000 rpm; the near-singular
    reference linkage has no feasible structure there at all.  That is the
    whole finding, in one comparison.
    """
    at_rest = solve(mechanism, 0.0, samples=SAMPLES)
    fast = solve(mechanism, 3000.0, samples=SAMPLES)
    assert fast.converged and at_rest.converged

    linkage = evaluate(COUPLED_DESIGN, speed_rpm=3000.0)
    assert not linkage.feasible


# -- parity of treatment ------------------------------------------------------


def test_both_mechanisms_are_sized_by_the_same_code() -> None:
    """The comparison is only meaningful if the structural model is identical.

    ``size_from_arrays`` takes the member list as a parameter precisely so that
    both mechanisms go through the same yield, fatigue and buckling checks
    rather than through two implementations that might differ.
    """
    from exlink.sizing import size_from_arrays

    axial = np.full((2, 64, 9), 5000.0)
    bending = np.full((2, 64, 9), 20000.0)
    lengths = np.array([30.0, 90.0])
    sized = size_from_arrays(
        axial, bending, lengths, fixity=np.array([2.0, 1.0]), names=("crank", "rod")
    )
    assert set(sized) == {"crank", "rod"}
    assert all(item.diameter > 0.0 for item in sized.values())


def test_the_slider_crank_budget_has_no_gears(solved: object) -> None:
    budget = mass_budget(solved)
    assert budget.items["gears"] == 0.0
    assert budget.gears is None
    assert budget.total_kg > 0.0


# -- the comparison -----------------------------------------------------------


def test_extended_expansion_buys_indicated_efficiency(mechanism: SliderCrank) -> None:
    """The Atkinson linkage's actual thermodynamic advantage, measured.

    Same compression ratio, larger expansion ratio, so a few points of
    indicated efficiency -- and only a few.  Quantifying how few is the point:
    it is much less than the popular account of extended expansion suggests.
    """
    otto = evaluate_slidercrank(mechanism, 1500.0)
    atkinson = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
    gain = atkinson.indicated_efficiency - otto.indicated_efficiency
    assert 0.0 < gain < 0.06


def test_the_range_advantage_is_firing_frequency_not_expansion() -> None:
    """The comparison's most load-bearing assumption, tested rather than assumed.

    In this model the EX-link fires every crankshaft revolution while a
    conventional four-stroke fires every other one, so per unit of work it
    accumulates half the journal rotation and half the piston sliding.  Re-run
    with that advantage removed, the range gain over a slider-crank very nearly
    disappears -- which means the headline result is about firing frequency,
    and must be reported as such.
    """
    atkinson = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
    sensitivity = firing_frequency_sensitivity(atkinson)
    otto = max(
        (
            evaluate_slidercrank(SliderCrank.for_compression_ratio(16.0), rpm)
            for rpm in (1500.0, 2000.0, 2500.0)
        ),
        key=lambda item: item.km_per_litre,
    )
    as_modelled = sensitivity["km_per_litre"] / otto.km_per_litre - 1.0
    four_stroke = sensitivity["km_per_litre_four_stroke"] / otto.km_per_litre - 1.0
    assert as_modelled > 0.15
    assert four_stroke < 0.10
    assert as_modelled > 3.0 * four_stroke


def test_friction_is_reported_per_cycle_not_per_revolution(solved: object) -> None:
    loss = friction_work(solved)
    assert 0.0 < loss < solved.indicated_work


def test_a_slider_crank_reaches_the_right_order_of_range(mechanism: SliderCrank) -> None:
    outcome = evaluate_slidercrank(mechanism, 2000.0)
    assert outcome.feasible
    assert 1000.0 < outcome.km_per_litre < 5000.0
    assert outcome.joints == 3
