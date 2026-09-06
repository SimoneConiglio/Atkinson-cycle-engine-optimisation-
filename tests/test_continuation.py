"""Climbing to the reliable region, one rung at a time."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.continuation import (
    DEFAULT_SCHEDULE,
    Continuation,
    ContinuationStep,
    format_continuation,
    reliability_continuation,
)
from exlink.design import Bounds
from exlink.reference import COUPLED_DESIGN, RANGE_DESIGN
from exlink.synthesis import target_from_design


def _step(
    target: float,
    beta: float,
    feasible: bool = True,
    moved: bool = True,
    strict: bool = True,
):
    return ContinuationStep(
        beta_target=target,
        design=COUPLED_DESIGN,
        km_per_litre=3300.0,
        steered_beta=beta,
        system_beta=beta - 0.1,
        system_pf=1.0e-3,
        worst_constraint=-1.0e-8,
        feasible=feasible,
        strictly_feasible=strict,
        moved=moved,
        evaluations=40,
        converged=True,
    )


def test_the_schedule_ascends_and_starts_below_the_target():
    assert list(DEFAULT_SCHEDULE) == sorted(DEFAULT_SCHEDULE)
    assert DEFAULT_SCHEDULE[0] < DEFAULT_SCHEDULE[-1]
    assert DEFAULT_SCHEDULE[-1] == pytest.approx(3.0)


def test_reaching_means_the_last_rung_was_met():
    climbed = Continuation(COUPLED_DESIGN, (_step(1.0, 1.5), _step(3.0, 3.02)))
    assert climbed.reached
    assert climbed.best.beta_target == pytest.approx(3.0)

    stalled = Continuation(COUPLED_DESIGN, (_step(1.0, 1.5), _step(3.0, 2.1)))
    assert not stalled.reached
    assert stalled.best.beta_target == pytest.approx(1.0)


def test_an_infeasible_rung_is_never_counted_as_met():
    """Holding beta while violating a geometric bound is not success."""
    outcome = Continuation(COUPLED_DESIGN, (_step(3.0, 3.5, feasible=False),))
    assert not outcome.reached
    assert outcome.best is None
    assert "infeasible" in format_continuation(outcome)


def test_feasibility_is_judged_against_the_bands_the_solve_used():
    """A rung solved at a widened band is not a failure for missing 0.05."""
    outcome = Continuation(COUPLED_DESIGN, (_step(3.0, 3.1, strict=False),))
    assert outcome.reached
    text = format_continuation(outcome)
    assert "met" in text
    assert "strict" not in text


def test_a_ladder_that_never_moved_is_the_failure_this_module_diagnoses():
    outcome = Continuation(COUPLED_DESIGN, (_step(0.0, -0.4, moved=False),))
    assert not any(step.moved for step in outcome.steps)
    assert "NO" in format_continuation(outcome)


def test_the_evaluations_are_summed_across_rungs():
    outcome = Continuation(COUPLED_DESIGN, (_step(1.0, 1.5), _step(2.0, 2.5)))
    assert outcome.evaluations == 80


def test_the_table_names_every_rung():
    outcome = Continuation(COUPLED_DESIGN, (_step(1.0, 1.5), _step(3.0, 3.1)))
    text = format_continuation(outcome)
    assert "1.00" in text
    assert "3.00" in text
    assert "highest rung met" in text


@pytest.mark.slow
def test_one_rung_runs_end_to_end():
    """A single short rung, to check the wiring rather than the result."""
    target = target_from_design(RANGE_DESIGN)
    outcome = reliability_continuation(
        target,
        RANGE_DESIGN,
        schedule=(0.0,),
        bounds=Bounds.around(RANGE_DESIGN, relative=0.30),
        max_iterations=4,
        module=0.8,
        teeth=48,
    )
    assert len(outcome.steps) == 1
    step = outcome.steps[0]
    assert np.isfinite(step.km_per_litre)
    assert step.evaluations > 0
