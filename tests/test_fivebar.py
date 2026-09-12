"""The geared five-bar: does the closed form describe the mechanism found?"""

from __future__ import annotations

import numpy as np
import pytest

from exlink import fivebar, topology
from exlink.fivebar import SYNTHESISED, AssemblyError, FiveBar, in_the_domain, solve


def test_the_closed_form_matches_the_spring_model() -> None:
    """Two independent solvers, one mechanism, agreeing to microns.

    The closed form here and the spring model of :mod:`exlink.topology` share no
    code: one intersects circles, the other minimises an elastic energy by
    Newton's method on a penalised ground structure.  Agreement to a few microns
    is the check that the five-bar written down here really is the mechanism the
    synthesis found, and what is left is the spring model's own compliance.
    """
    _structure, layout, x = in_the_domain(SYNTHESISED)
    motion = topology.sweep(layout, x, samples=720, penalty=topology.PENALTY_SCHEDULE[-1])
    closed = solve(SYNTHESISED, samples=720)

    assert float(np.max(np.abs(closed.lam - motion.lam))) < 0.01, "the same motion"
    assert motion.strain < 1.0e-4, "it runs as a mechanism, not a strained structure"
    assert motion.output_slack < topology.SLACK_TOLERANCE, "the input determines the piston"


def test_writing_it_into_the_domain_and_reading_it_back_is_exact() -> None:
    """A round trip through the ground structure changes no dimension."""
    _structure, _layout, x = in_the_domain(SYNTHESISED)
    recovered, offset = fivebar.from_topology(x)

    assert offset == pytest.approx(0.0, abs=1.0e-12), "the datum is the crank pin"
    assert np.allclose(recovered.to_array(), SYNTHESISED.to_array(), atol=1.0e-9)


def test_it_meets_the_specification_it_was_synthesised_for() -> None:
    """The strokes and the compression ratio §5.9 reports, from the closed form."""
    closed = solve(SYNTHESISED, samples=1440)
    assert topology.requirement_error(closed.lam) < 0.6, "the four conditions"
    assert float(np.ptp(closed.lam)) == pytest.approx(74.12, abs=0.05), "the stroke"


def test_and_is_mechanically_indefensible_as_synthesised() -> None:
    """The finding that motivates pricing it: the rod is nowhere near the axis.

    Not a defect of the topology -- a defect of judging a mechanism by its
    motion.  A conventional engine keeps this ratio below about 0.3.
    """
    closed = solve(SYNTHESISED, samples=1440)
    assert closed.side_load_ratio > 1.5, "the piston is pressed into the liner"
    assert float(np.degrees(np.max(np.abs(closed.rod_angle)))) > 55.0


def test_a_linkage_that_cannot_close_says_so() -> None:
    """Not a worse design: not a design.  So it raises rather than returning nan."""
    with pytest.raises(AssemblyError, match="do not reach each other"):
        solve(SYNTHESISED.replace(L_1=5.0), samples=90)
    with pytest.raises(AssemblyError, match="cylinder axis"):
        solve(SYNTHESISED.replace(e=10.0), samples=90)


def test_the_design_vector_round_trips_through_its_representations() -> None:
    """Array, mapping and dataclass all carry the same nine numbers."""
    values = SYNTHESISED.to_array()
    assert values.size == len(fivebar.VARIABLE_NAMES)
    assert np.allclose(FiveBar.from_array(values).to_array(), values)
    assert np.allclose(FiveBar.from_mapping(SYNTHESISED.to_mapping()).to_array(), values)
    with pytest.raises(ValueError, match="expected 9 design variables"):
        FiveBar.from_array(np.zeros(3))
    with pytest.raises(ValueError, match="unknown design variables"):
        SYNTHESISED.replace(nonsense=1.0)


def test_every_variable_is_described() -> None:
    """A reader of the results table should not have to read the source."""
    assert set(fivebar.VARIABLE_DESCRIPTIONS) == set(fivebar.VARIABLE_NAMES)
    assert set(fivebar.VARIABLE_NAMES) >= fivebar.ANGULAR_VARIABLES


def _loaded(speed_rpm: float, samples: int = 720) -> fivebar.FiveBarLoads:
    """The synthesised five-bar, sized crudely and loaded at one speed."""
    from exlink import cycle
    from exlink.constants import DEFAULT_SPEC
    from exlink.materials import DEFAULT_MATERIAL

    motion = solve(SYNTHESISED, samples=samples)
    thermo = cycle.solve(motion.lam, DEFAULT_SPEC)
    diameters = dict.fromkeys([member[0] for member in fivebar.MEMBERS], 12.0)
    properties = fivebar.mass_properties(
        motion,
        SYNTHESISED,
        diameters,
        DEFAULT_MATERIAL.density,
        piston_mass=2.0e-4,
        piston_length=DEFAULT_SPEC.piston_length,
    )
    return fivebar.solve_loads(
        motion,
        SYNTHESISED,
        thermo.piston_force,
        properties,
        speed_rpm * 2.0 * np.pi / 60.0,
        pressure_angle=DEFAULT_SPEC.pressure_angle,
        piston_length=DEFAULT_SPEC.piston_length,
    )


@pytest.mark.parametrize("rpm", [0.0, 3000.0])
def test_the_equilibrium_conserves_power(rpm: float) -> None:
    """Shaft work must equal the p-V loop area, and that is what fixes the signs.

    The two are computed by different routes: one integrates the unknown torque
    the equilibrium returns, the other integrates the applied gas force against
    the piston velocity.  Nothing forces them to agree unless every sign and
    every moment arm in the eighteen equations is right.  It caught the torque
    sign, which came out reading minus one.
    """
    loads = _loaded(rpm)
    assert loads.shaft_work / loads.indicated_work == pytest.approx(1.0, abs=1.0e-3)


def test_inertia_does_no_net_work_over_a_cycle() -> None:
    """A second check, and an independent one: the balance cannot depend on speed.

    Inertia forces are internal to a closed cycle, so spinning the mechanism up
    changes every reaction but not the work delivered.  If the inertia terms were
    entered with a wrong arm, the two speeds would disagree.
    """
    still, running = _loaded(0.0), _loaded(3000.0)
    assert running.shaft_work == pytest.approx(still.shaft_work, rel=1.0e-9)
    assert running.peak_bearing_load != pytest.approx(still.peak_bearing_load, rel=1.0e-3)


def test_the_inertia_relief_is_the_right_way_round() -> None:
    """Spinning it up unloads the journals and loads the gear teeth.

    Both follow from where the mass is: the reciprocating inertia opposes the gas
    force near top dead centre, which is where the bearings see their worst load,
    while the gear pair has to accelerate the crank throws and so sees more.
    """
    still, running = _loaded(0.0), _loaded(3000.0)
    assert running.peak_bearing_load < still.peak_bearing_load
    assert np.max(np.abs(running.gear_force)) > np.max(np.abs(still.gear_force))


def test_the_liner_load_is_as_bad_as_the_geometry_promised() -> None:
    """A side-load ratio of 1.77 is not an abstraction; it is kilonewtons."""
    loads = _loaded(3000.0)
    assert float(np.max(np.abs(loads.liner_force))) > 3000.0, "thousands of newtons"
    assert loads.conditioning < 1.0e8, "the system is still solvable"


def test_the_five_bar_saves_members_but_not_journals() -> None:
    """Seven journals, exactly the EX-link's, because the floating pin carries two.

    The naive count says a five-member mechanism must rub less than a
    seven-member one.  It does not: three links meet at the floating pin, which
    in metal is two bearings side by side, each turning at its own relative
    speed.  Counting that pin once would hand the five-bar a friction advantage
    it has not got.
    """
    from exlink.dynamics import MEMBERS as EXLINK_MEMBERS

    assert len(fivebar.MEMBERS) == 5
    assert len(EXLINK_MEMBERS) == 7
    assert len(fivebar.JOURNALS) == 7
    assert sum(1 for key, _a, _b in fivebar.JOURNALS if key.startswith("F")) == 2


def test_the_chain_to_range_is_internally_consistent() -> None:
    """Brake efficiency is the product of the two it is built from."""
    performance = fivebar.evaluate(SYNTHESISED, 1000.0, samples=360)
    assert performance.result.converged
    assert performance.brake_efficiency == pytest.approx(
        performance.indicated_efficiency * performance.mechanical_efficiency, rel=1.0e-9
    )
    assert performance.range_km_per_litre > 0.0
    assert performance.engine_mass > 0.0


def test_the_side_load_is_where_the_five_bar_loses() -> None:
    """Friction at the liner, not at the bearings, is what eats this design.

    The reason to carry a synthesised mechanism through to physics rather than
    stopping at its motion: the piston rubs away most of the indicated work,
    because the connecting rod the fit chose runs 60 degrees off the cylinder.
    """
    result = fivebar.solve_sized(SYNTHESISED, 1000.0, samples=360)
    losses = fivebar.friction_work(result)
    assert losses["piston"] > losses["bearings"], "the liner dominates"
    assert losses["piston"] > 0.5 * result.indicated_work, "over half the output"


def test_the_mass_spiral_diverges_above_about_a_thousand_rpm() -> None:
    """A physical result, reported rather than hidden behind a converged flag.

    Heavier members carry more inertia, which needs heavier members.  For this
    geometry the loop closes below roughly 1000 rpm and runs away above it, so the
    mechanism as synthesised cannot be built for the study's own 2000 rpm
    operating point.  The flag says so instead of returning a fixed point that is
    not one.
    """
    assert fivebar.solve_sized(SYNTHESISED, 1000.0, samples=360).converged
    fast = fivebar.solve_sized(SYNTHESISED, 2000.0, samples=360)
    assert not fast.converged
    assert fast.diameters["con_rod"] > 2.0 * 19.85, "it is running away, not merely rough"
