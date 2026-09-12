"""The five-bar's optimization: is the problem posed so that a solve means something?"""

from __future__ import annotations

import numpy as np
import pytest

from exlink import fivebar, topology
from exlink import fivebar_scenarios as scenarios
from exlink.fivebar import SYNTHESISED


def test_an_unanalysable_design_is_never_reported_feasible() -> None:
    """The bug that quietly ruined the first solve, pinned.

    Returning zero for the constraints of a design that cannot be analysed makes
    that design satisfy every one of them, and GEMSEO duly declared a point
    "feasible" whose side-load ratio was 1.27 and whose range did not exist.
    Every failure branch must report its constraints violated.
    """
    unassemblable = scenarios.analyse(SYNTHESISED.replace(L_1=5.0), 2000.0, samples=90)
    assert unassemblable["assembly_margin"] > 0.0
    assert unassemblable["side_load_margin"] > 0.0
    assert unassemblable["neg_range"] == pytest.approx(scenarios.CYCLE_PENALTY)

    # And a design that sizes but whose mass spiral never closed is not a
    # structure either, however good its geometry looks.
    stalled = scenarios.analyse(SYNTHESISED, 3000.0, samples=120)
    assert stalled["converged"] == 0.0
    assert stalled["assembly_margin"] > 0.0


def test_the_ladder_never_goes_flat() -> None:
    """Every rung carries a gradient, or the search stops standing on it.

    Three rungs below a computable range -- sizes but makes no net work, does not
    size at all, is not a cycle -- and each is ordered below the one above it.
    """
    running = scenarios.analyse(SYNTHESISED, 1000.0, samples=180)
    stalled = scenarios.analyse(SYNTHESISED, 3000.0, samples=120)
    broken = scenarios.analyse(SYNTHESISED.replace(L_1=5.0), 1000.0, samples=90)

    assert running["neg_range"] < 0.0, "a real range is the best rung"
    assert running["neg_range"] < stalled["neg_range"] < broken["neg_range"]


def test_restoring_feasibility_holds_the_whole_cycle() -> None:
    """What the restoration is for, and what it cannot do.

    It can hold the four-stroke cycle -- all five precision-point conditions, not
    just the expansion stroke and the compression ratio -- while keeping the
    mechanism clear of its transmission-angle singularity.  It cannot bring the
    side load inside the limit a conventional engine lives by, and that is a
    property of the topology rather than of the search: with the cycle held the
    best of thirty-two starts reached 1.46, while dropping the cycle constraint
    reaches 0.045.

    Both halves are asserted, because a test that only demanded the cycle would
    pass on the singular design this constraint set exists to exclude.
    """
    restored = scenarios.restore_feasibility(samples=360, band=0.5)
    motion = fivebar.solve(restored, samples=720)
    residual = topology.requirement_residuals(motion.lam)

    assert topology.requirement_error(motion.lam) <= 0.7, "the cycle is held"
    assert abs(residual[2]) * 74.0 < 1.0, "both top dead centres, level to a millimetre"
    assert motion.transmission >= scenarios.MIN_TRANSMISSION - 1.0e-3, "still a mechanism"
    assert motion.closure > 0.0 and motion.reach > 0.0, "and still assembled"
    assert motion.side_load_ratio > scenarios.SIDE_LOAD_LIMIT, (
        "and the side load is *not* fixed: holding the cycle costs it"
    )


def test_two_constraints_do_not_specify_a_four_stroke_cycle() -> None:
    """Why the whole requirement is constrained and not the study's usual two.

    The EX-link's formulation holds the expansion stroke and the compression
    ratio, which suffices there because its topology makes the two top dead
    centres level as a matter of geometry.  A five-bar's does not, and asked for
    only those two the restoration dropped one top dead centre 47 mm below the
    other -- a motion ``find_phases`` still accepts, and whose peak pressure then
    came out at 0.098 MPa against the 5.5 it should be.
    """
    from exlink import cycle
    from exlink.constants import DEFAULT_SPEC

    # What the restoration actually returned when it was held to those two
    # conditions alone, kept verbatim.
    broken = fivebar.FiveBar(
        q_1=31.898987, q_2=19.511763, I=160.0, theta_r=5.867512,
        theta_f=-103.647164, L_1=87.049493, L_2=121.849183, e=340.0,
        x_1=66.951712, branch=-1, piston_branch=1,
    )  # fmt: skip
    motion = fivebar.solve(broken, samples=720)
    residual = topology.requirement_residuals(motion.lam)
    phases = cycle.find_phases(motion.lam)
    thermo = cycle.solve(motion.lam, DEFAULT_SPEC, phases)

    # It passes the two constraints the EX-link is solved under...
    assert abs(phases.expansion_stroke - 74.0) < 0.5
    assert abs(thermo.compression_ratio - 16.0) < 0.5
    # ...and is not a four-stroke engine at all.
    assert abs(residual[2]) * 74.0 > 40.0, "its top dead centres differ by tens of mm"
    assert float(np.max(thermo.gauge_pressure)) < 0.2, "and it barely compresses anything"


SYNTHESISED_SIDE_LOAD = 1.7738


def test_the_synthesised_design_was_already_past_alts_limit() -> None:
    """Its transmission angle was below 30 degrees before anything was optimised.

    Worth recording: the precision-point fit did not merely pick a bad rod angle,
    it picked a geometry whose driving dyad was already marginal.
    """
    motion = fivebar.solve(SYNTHESISED, samples=720)
    assert motion.transmission < scenarios.MIN_TRANSMISSION
    assert motion.side_load_ratio == pytest.approx(SYNTHESISED_SIDE_LOAD, abs=1.0e-3)


def test_every_bound_names_a_real_variable() -> None:
    """A typo in the design box would silently leave a variable unbounded."""
    assert set(scenarios.BOUNDS) == set(fivebar.VARIABLE_NAMES)
    for name, (low, high) in scenarios.BOUNDS.items():
        assert low < high, name
        assert low <= getattr(SYNTHESISED, name) <= high, f"{name} starts inside its box"


def test_the_outcome_calls_an_infeasible_answer_infeasible() -> None:
    """A solve that ends on a design that does not size is not a solution."""
    outcome = scenarios.FiveBarOutcome(
        design=SYNTHESISED, performance=None, evaluations=0, message=""
    )
    assert not outcome.feasible
    assert not np.isnan(outcome.design.q_1)
