"""The sizing / dynamics fixed point, and its GEMSEO counterpart."""

from __future__ import annotations

import numpy as np
import pytest

from exlink import analyse
from exlink.coupled import solve_coupled, solve_for_design
from exlink.disciplines import COUPLED_SAMPLES, DynamicsDiscipline, StructureDiscipline
from exlink.dynamics import MEMBER_NAMES
from exlink.reference import REFINED_DESIGN

VIABLE_RPM = 1000.0
"""A speed at which the reference design can still be built."""


@pytest.fixture(scope="module")
def at_rest():
    return solve_for_design(REFINED_DESIGN, speed_rpm=0.0, samples=180)


@pytest.fixture(scope="module")
def at_speed():
    return solve_for_design(
        REFINED_DESIGN, speed_rpm=VIABLE_RPM, samples=180, max_iterations=400
    )


def test_there_is_no_coupling_at_rest(at_rest) -> None:
    """Without inertia the masses do not feed back, so one sweep settles it.

    This is the geometric problem: sizing is a post-processing step, not a
    loop. It takes a second sweep only to confirm nothing moved.
    """
    assert at_rest.converged
    assert at_rest.iterations <= 3
    assert at_rest.feasible


def test_the_loop_is_real_once_the_engine_turns(at_speed, at_rest) -> None:
    """At speed the sections drive the loads that drive the sections."""
    assert at_speed.converged
    assert at_speed.iterations > 10
    assert at_speed.total_mass_kg > 2.0 * at_rest.total_mass_kg


def test_mass_grows_steeply_with_speed() -> None:
    """Scaling says ``m ~ (C a)^3``; the numbers should be at least super-linear."""
    slow = solve_for_design(REFINED_DESIGN, speed_rpm=500.0, samples=180, max_iterations=400)
    fast = solve_for_design(REFINED_DESIGN, speed_rpm=1000.0, samples=180, max_iterations=400)
    assert fast.total_mass_kg > 2.0 * slow.total_mass_kg


def test_the_loop_runs_away_at_high_speed() -> None:
    """Past a point no section is thick enough to carry its own inertia.

    That has to be reported, not silently returned as a design: a saturated
    result means the mechanism cannot be built at that speed.
    """
    result = solve_for_design(REFINED_DESIGN, speed_rpm=3000.0, samples=180, max_iterations=200)
    assert result.saturated
    assert not result.feasible


def test_relaxation_reaches_the_same_fixed_point(at_speed) -> None:
    """Under-relaxation changes the path to the answer, not the answer."""
    damped = solve_for_design(
        REFINED_DESIGN,
        speed_rpm=VIABLE_RPM,
        samples=180,
        relaxation=0.5,
        max_iterations=800,
    )
    assert damped.converged
    for name in MEMBER_NAMES:
        assert damped.diameters[name] == pytest.approx(at_speed.diameters[name], rel=1e-4)


def test_relaxation_is_validated() -> None:
    solved = analyse(REFINED_DESIGN, samples=180).require_solved()
    with pytest.raises(ValueError, match=r"must lie in \(0, 1\]"):
        solve_coupled(solved.kinematics, solved.thermodynamics, relaxation=0.0)


def test_solving_an_unanalysable_design_raises() -> None:
    with pytest.raises(ValueError, match="unanalysable"):
        solve_for_design(REFINED_DESIGN.replace(a=25.0, c=25.0), samples=180)


def test_a_warm_start_costs_fewer_sweeps(at_speed) -> None:
    warm = solve_for_design(
        REFINED_DESIGN,
        speed_rpm=VIABLE_RPM,
        samples=180,
        initial_diameters=at_speed.diameters,
        max_iterations=400,
    )
    assert warm.iterations < at_speed.iterations


def test_every_member_is_fully_utilised(at_speed) -> None:
    name, utilisation = at_speed.worst_utilisation()
    assert utilisation == pytest.approx(1.0, abs=1e-3), name


@pytest.mark.slow
def test_the_gemseo_mda_matches_the_reference_solver() -> None:
    """The two implementations of the same fixed point must agree.

    :func:`~exlink.coupled.solve_coupled` is a hand-written Gauss-Seidel sweep;
    the MDA is GEMSEO driving the same two disciplines. Agreeing to a micron
    means the discipline wrappers carry the physics faithfully.
    """
    from gemseo import create_mda

    mda = create_mda(
        "MDAGaussSeidel",
        [DynamicsDiscipline(speed_rpm=VIABLE_RPM), StructureDiscipline()],
        tolerance=1e-10,
        max_mda_iter=400,
    )
    output = mda.execute(REFINED_DESIGN.to_mapping())
    reference = solve_for_design(
        REFINED_DESIGN, speed_rpm=VIABLE_RPM, samples=COUPLED_SAMPLES, max_iterations=400
    )
    expected = np.array([reference.diameters[n] for n in MEMBER_NAMES])
    assert output["diameters"] == pytest.approx(expected, abs=1e-4)


@pytest.mark.slow
def test_the_mda_couples_in_both_directions() -> None:
    """Both variables must be strong couplings, or it is not an MDA at all."""
    from gemseo import create_mda

    mda = create_mda(
        "MDAGaussSeidel",
        [DynamicsDiscipline(speed_rpm=VIABLE_RPM), StructureDiscipline()],
    )
    couplings = set(mda.coupling_structure.strong_couplings)
    assert "diameters" in couplings
    assert {"member_axial", "member_bending"} & couplings


def test_the_coupled_sample_count_resolves_the_tdc_gap() -> None:
    """``g`` must be measured accurately enough to optimize against.

    It is the difference of two nearly equal maxima, so its absolute error is
    what counts against a 0.01 mm bound -- and on a coarse grid that error is
    not small. This pins the default the coupled disciplines run at against a
    converged reference.
    """
    converged = analyse(REFINED_DESIGN, samples=2880).metrics.tdc_gap
    coupled = analyse(REFINED_DESIGN, samples=COUPLED_SAMPLES).metrics.tdc_gap
    coarse = analyse(REFINED_DESIGN, samples=180).metrics.tdc_gap

    assert abs(coupled - converged) < 1e-3
    # The grid that was not good enough, kept as the reason for the one that is.
    assert abs(coarse - converged) > 1e-3
