"""The five-bar's optimization: is the problem posed so that a solve means something?"""

from __future__ import annotations

import numpy as np
import pytest

from exlink import fivebar
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


def test_restoring_feasibility_fixes_the_geometry_the_fit_chose() -> None:
    """The side load comes down by more than an order of magnitude, in a second.

    And it stays a mechanism while doing it, which is the whole point of the
    transmission-angle constraint: the first attempt without one returned a design
    with a beautiful 0.045 side load sitting on its own singularity.
    """
    restored = scenarios.restore_feasibility(samples=240)
    motion = fivebar.solve(restored, samples=720)

    assert motion.side_load_ratio < scenarios.SIDE_LOAD_LIMIT, "inside the limit"
    assert motion.side_load_ratio < 0.2 * SYNTHESISED_SIDE_LOAD, "and far better than the start"
    assert motion.transmission >= scenarios.MIN_TRANSMISSION - 1.0e-3, "still a mechanism"
    assert motion.closure > 0.0 and motion.reach > 0.0, "and still assembled"


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
