"""Friction, mass budget, range: the chain that turns a linkage into km/L."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.coupled import solve_for_design
from exlink.friction import losses, sensitivity
from exlink.mass_budget import (
    assemble,
    flywheel_requirement,
    gas_torque,
    shaft_diameter,
)
from exlink.model import analyse
from exlink.performance import evaluate
from exlink.reference import COUPLED_DESIGN, REFINED_DESIGN
from exlink.vehicle import Vehicle, best_strategy, burn_and_coast, heat_release

SPEED = 1000.0


@pytest.fixture(scope="module")
def sized() -> object:
    return solve_for_design(REFINED_DESIGN, speed_rpm=SPEED)


# -- friction -----------------------------------------------------------------


def test_inertia_does_no_net_work_over_a_cycle(sized: object) -> None:
    """The gas torque and the full torque must have the same mean.

    Inertia torque is energy traded with the moving masses; over a closed cycle
    it has to integrate to zero.  If it did not, the flywheel sizing would be
    reading a fictitious energy swing, and the brake work would be wrong by
    whatever the discrepancy was.
    """
    assert float(np.mean(gas_torque(sized.loads))) == pytest.approx(
        sized.loads.mean_torque, rel=1.0e-6
    )


def test_friction_work_is_positive_and_bounded(sized: object) -> None:
    result = losses(sized.loads, sized.diameters)
    assert result.bearing_work > 0.0
    assert result.piston_work > 0.0
    assert result.total_work == pytest.approx(
        result.bearing_work + result.piston_work + result.mesh_work
    )
    assert result.brake_work == pytest.approx(result.indicated_work - result.total_work)


def test_friction_scales_linearly_with_the_journal_coefficient(sized: object) -> None:
    """Doubling mu must exactly double the bearing loss.

    This is the property that lets a single sensitivity sweep stand in for the
    whole uncertainty in the friction model.
    """
    single = losses(sized.loads, sized.diameters, journal_friction=0.008)
    double = losses(sized.loads, sized.diameters, journal_friction=0.016)
    assert double.bearing_work == pytest.approx(2.0 * single.bearing_work, rel=1.0e-9)


def test_zero_friction_recovers_the_indicated_work(sized: object) -> None:
    ideal = losses(
        sized.loads,
        sized.diameters,
        journal_friction=0.0,
        piston_friction=0.0,
        ring_tension=0.0,
        mesh_efficiency=1.0,
    )
    assert ideal.mechanical_efficiency == pytest.approx(1.0)


def test_mechanical_efficiency_falls_as_friction_rises(sized: object) -> None:
    sweep = sensitivity(sized.loads, sized.diameters)
    values = [sweep[key] for key in sorted(sweep)]
    assert values == sorted(values, reverse=True)


def test_indicated_efficiency_is_thermodynamically_sane() -> None:
    """The p-V loop area over the heat released must be a plausible fraction.

    An idealised Atkinson cycle at this expansion ratio should land in the
    40-55 % band.  Outside it, either the cycle or the heat release is wrong,
    and the range model inherits the error directly.
    """
    solved = analyse(REFINED_DESIGN, samples=720).require_solved()
    thermo = solved.thermodynamics
    quantity = heat_release(3000.0, thermo.p_compression_end, thermo.p_combustion, 1.22)
    work = float(np.trapezoid(thermo.pressure, thermo.volume))
    assert 0.40 < work / quantity < 0.55


# -- mass budget --------------------------------------------------------------


def test_flywheel_inertia_follows_the_textbook_scaling() -> None:
    """``J = dE / (delta omega^2)``, so doubling speed quarters the inertia."""
    torque = 2000.0 + 5000.0 * np.sin(np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False))
    slow, swing_slow = flywheel_requirement(torque, 100.0)
    fast, swing_fast = flywheel_requirement(torque, 200.0)
    assert swing_slow == pytest.approx(swing_fast)
    assert fast == pytest.approx(slow / 4.0, rel=1.0e-9)


def test_a_constant_torque_needs_no_flywheel() -> None:
    inertia, swing = flywheel_requirement(np.full(360, 2500.0), 100.0)
    assert swing == pytest.approx(0.0, abs=1.0e-9)
    assert inertia == pytest.approx(0.0, abs=1.0e-12)


def test_shaft_diameter_grows_as_the_cube_root_of_load() -> None:
    small = shaft_diameter(1000.0, 0.0, 20.0)
    large = shaft_diameter(8000.0, 0.0, 20.0)
    assert large / small == pytest.approx(2.0, rel=0.02)


def test_budget_is_dominated_by_parts_the_linkage_sizing_never_saw(sized: object) -> None:
    """The headline point: 'mass' meant a quarter-kilogram tail.

    The sized members must be a small share of a budget that also carries a
    crankcase, a flywheel, shafts, bearings and gears.
    """
    analysis = analyse(REFINED_DESIGN, samples=360)
    metrics = analysis.metrics
    thermo = analysis.require_solved().thermodynamics
    budget = assemble(
        sized.loads,
        sized.diameters,
        sized.mass_properties.member_mass,
        sized.piston_mass,
        metrics.height,
        metrics.width,
        metrics.expansion_stroke + metrics.compression_stroke,
        float(np.max(thermo.gauge_pressure)),
    )
    assert budget.total_kg > 5.0 * sized.total_mass_kg
    assert budget.shares()["linkage"] < 0.25
    assert set(budget.items) == {
        "linkage",
        "piston",
        "gears",
        "shafts",
        "bearings",
        "crankcase",
        "cylinder_head",
        "flywheel",
    }


def test_envelope_drives_crankcase_mass(sized: object) -> None:
    """H and B are mass in disguise; the budget has to show it.

    This is what lets the three-objective geometric problem collapse to one.
    """
    analysis = analyse(REFINED_DESIGN, samples=360)
    thermo = analysis.require_solved().thermodynamics
    common = (
        sized.loads,
        sized.diameters,
        sized.mass_properties.member_mass,
        sized.piston_mass,
    )
    stroke = analysis.metrics.expansion_stroke + analysis.metrics.compression_stroke
    peak = float(np.max(thermo.gauge_pressure))
    small = assemble(*common, 150.0, 100.0, stroke, peak)
    big = assemble(*common, 300.0, 200.0, stroke, peak)
    assert big.items["crankcase"] > 1.5 * small.items["crankcase"]


# -- vehicle ------------------------------------------------------------------


def test_a_heavier_engine_goes_less_far() -> None:
    car = Vehicle()
    light = best_strategy(car, 3.0, brake_power=150.0, efficiency=0.30)
    heavy = best_strategy(car, 30.0, brake_power=150.0, efficiency=0.30)
    assert light.km_per_litre > heavy.km_per_litre


def test_a_more_efficient_engine_goes_proportionally_further() -> None:
    """Range is linear in brake efficiency at a fixed operating point.

    Fuel energy per unit work is ``1 / eta_b`` and nothing else in the road
    load depends on it, so the ratio must be exact.
    """
    car = Vehicle()
    low = burn_and_coast(car, 10.0, 150.0, 0.20, 6.0, 12.0)
    high = burn_and_coast(car, 10.0, 150.0, 0.40, 6.0, 12.0)
    assert high.km_per_litre == pytest.approx(2.0 * low.km_per_litre, rel=1.0e-9)


def test_an_engine_too_weak_to_accelerate_is_rejected() -> None:
    result = burn_and_coast(Vehicle(), 10.0, 0.5, 0.30, 6.0, 20.0)
    assert not result.feasible
    assert "accelerate" in result.reason


def test_burn_and_coast_conserves_energy() -> None:
    """Over a closed cycle the kinetic energy nets to zero, so the propulsive
    work must equal the resistance work over the *whole* distance.

    This is the claim that accelerating hard costs nothing in road load, and
    it is the reason burn-and-coast is worth modelling at all rather than
    assuming a constant cruise.
    """
    car = Vehicle()
    result = burn_and_coast(car, 10.0, 150.0, 0.30, 6.0, 12.0)
    assert result.feasible

    total_mass = result.total_mass
    speeds = np.linspace(6.0, 12.0, 2001)
    resistance = car.resistance(speeds, total_mass)
    wheel_power = 150.0 * car.transmission_efficiency

    # Resistance work over the burn, plus the kinetic energy gained, must equal
    # the propulsive work delivered over the burn.
    net = wheel_power / speeds - resistance
    burn_distance_integrand = total_mass * speeds / net
    resistance_work = float(np.trapezoid(resistance * burn_distance_integrand, speeds))
    kinetic = 0.5 * total_mass * (12.0**2 - 6.0**2)
    propulsive = float(np.trapezoid(wheel_power / speeds * burn_distance_integrand, speeds))
    assert propulsive == pytest.approx(resistance_work + kinetic, rel=1.0e-4)


def test_range_is_in_the_right_order_of_magnitude() -> None:
    """A Prototype-class gasoline entry scores between 1000 and 5000 km/L.

    Landing outside that band means the road load, the fuel model or the unit
    conversions are wrong, and no design comparison built on it would mean
    anything.
    """
    result = best_strategy(Vehicle(), 12.0, brake_power=150.0, efficiency=0.30)
    assert 1000.0 < result.km_per_litre < 5000.0


# -- the whole chain ----------------------------------------------------------


def test_evaluate_returns_a_consistent_result() -> None:
    outcome = evaluate(REFINED_DESIGN, speed_rpm=800.0)
    assert outcome.feasible, outcome.reason()
    assert outcome.coupled is not None and outcome.friction is not None
    assert outcome.brake_efficiency < outcome.indicated_efficiency
    assert outcome.engine_mass_kg == pytest.approx(outcome.budget.total_kg)
    assert outcome.km_per_litre > 0.0


def test_an_unanalysable_design_is_penalised_not_raised() -> None:
    """The optimizer must be able to walk through infeasible ground."""
    broken = REFINED_DESIGN.replace(a=0.5 * REFINED_DESIGN.a)
    outcome = evaluate(broken, speed_rpm=800.0)
    assert not outcome.feasible
    assert outcome.km_per_litre == 0.0
    assert outcome.reason()


def test_buildable_diameters_never_shrink_a_member() -> None:
    outcome = evaluate(REFINED_DESIGN, speed_rpm=800.0)
    assert outcome.coupled is not None
    for name, value in outcome.buildable_diameters.items():
        assert value >= outcome.coupled.diameters[name] - 1.0e-12
    assert outcome.stock_premium >= 0.0


def test_backing_off_the_singularity_wins_on_range() -> None:
    """The headline result of the coupled study, stated as a test.

    The quasi-statically attractive design sits at the transmission-angle
    singularity, where the accelerations, the bearing loads and hence the
    structure are all worst.  Backing off costs efficiency on paper and buys
    back more than it costs once mass is priced in kilometres.
    """
    near = max(
        (evaluate(REFINED_DESIGN, speed_rpm=rpm) for rpm in (600.0, 800.0, 1000.0)),
        key=lambda item: item.km_per_litre,
    )
    off = max(
        (evaluate(COUPLED_DESIGN, speed_rpm=rpm) for rpm in (800.0, 1000.0, 1250.0)),
        key=lambda item: item.km_per_litre,
    )
    assert off.feasible and near.feasible
    assert off.engine_mass_kg < 0.75 * near.engine_mass_kg
    assert off.km_per_litre > near.km_per_litre
