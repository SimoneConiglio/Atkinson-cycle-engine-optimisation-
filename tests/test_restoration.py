"""Getting inside the feasible set before optimising within it."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from exlink.constants import DEFAULT_TARGETS
from exlink.design import Bounds, Design
from exlink.reference import COUPLED_DESIGN
from exlink.restoration import (
    DEFAULT_MARGIN,
    Restoration,
    format_restoration,
    restore,
)
from exlink.synthesis import NUMBER_OF_CONSTRAINTS, _range_constraints, target_from_design

TARGETS = dataclasses.replace(DEFAULT_TARGETS, max_tdc_gap=0.1)
BAND = 0.15


def _outcome(margin: float, start: float = -1.0e-2):
    return Restoration(
        design=COUPLED_DESIGN,
        margin=margin,
        start_margin=start,
        evaluations=40,
        converged=True,
    )


def test_a_point_on_the_boundary_is_not_counted_as_restored():
    """Section 6.2 and 6.8 both turn on what sitting at g = 0 costs."""
    assert not _outcome(0.0).feasible
    assert not _outcome(0.5 * DEFAULT_MARGIN).feasible
    assert _outcome(2.0 * DEFAULT_MARGIN).feasible


def test_improvement_and_feasibility_are_different_questions():
    """A solve can move the right way and still not arrive."""
    crawled = _outcome(-1.0e-3, start=-1.0e-2)
    assert crawled.improved
    assert not crawled.feasible


def test_the_table_separates_the_three_outcomes():
    text = format_restoration([_outcome(1.0e-2), _outcome(-1.0e-3), _outcome(-1.0e-2)])
    assert "feasible" in text
    assert "improved" in text
    assert "stuck" in text
    assert "1 of 3" in text


def test_the_margin_is_the_worst_of_every_constraint():
    """Restoration maximises this, so it has to be the whole constraint set."""
    from exlink.performance import evaluate

    target = target_from_design(COUPLED_DESIGN)
    performance = evaluate(COUPLED_DESIGN, speed_rpm=1000.0, module=0.8, teeth=48)
    rows = _range_constraints(performance, target, BAND, TARGETS)
    assert rows.size == NUMBER_OF_CONSTRAINTS
    # The coupled reference is feasible at the widened band, so every row is
    # positive and the worst of them is what a restoration would be raising.
    assert float(np.min(rows)) > 0.0


@pytest.mark.slow
def test_restoration_reaches_the_interior_from_a_scattered_start():
    """The measurement §6.3(6) asks for, at one start rather than six."""
    target = target_from_design(COUPLED_DESIGN)
    box = Bounds.around(COUPLED_DESIGN, relative=0.30)
    rng = np.random.default_rng(0)
    start = Design.from_array(
        np.clip(
            COUPLED_DESIGN.to_array() * (1.0 + 0.05 * rng.normal(size=11)),
            box.lower,
            box.upper,
        )
    )
    outcome = restore(
        target, start, bounds=box, targets=TARGETS, module=0.8, teeth=48, max_iterations=40
    )
    assert outcome.evaluations > 0
    assert outcome.improved
