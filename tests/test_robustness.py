"""Manufacturing tolerance, propagated two ways, on a design that sits near a
singularity."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.design import ANGULAR_VARIABLES, VARIABLE_NAMES
from exlink.reference import COUPLED_DESIGN, REFINED_DESIGN
from exlink.robustness import (
    CONSTRAINT_NAMES,
    DEFAULT_GRADE,
    IT_FACTORS,
    covariance,
    required_grade,
    tolerance_half_widths,
    tolerance_report,
    tolerance_unit,
)

SAMPLES = 250
CRANK = 360


@pytest.fixture(scope="module")
def report() -> object:
    return tolerance_report(REFINED_DESIGN, samples=SAMPLES, crank_samples=CRANK)


# -- the tolerance model itself -----------------------------------------------


def test_tolerance_unit_matches_iso_286() -> None:
    """``i = 0.45 D^(1/3) + 0.001 D`` micrometres, spot-checked.

    At a 50 mm nominal size the standard gives i = 1.56 um and an IT8 band of
    39 um; reproducing the published value is what makes the rest of this
    module a tolerance study rather than a guess.
    """
    assert tolerance_unit(50.0) * 1000.0 == pytest.approx(1.71, abs=0.05)
    band = IT_FACTORS[8] * tolerance_unit(50.0)
    assert band == pytest.approx(0.039, abs=0.006)


def test_half_widths_grow_with_the_grade_number() -> None:
    tight = tolerance_half_widths(REFINED_DESIGN, grade=6)
    loose = tolerance_half_widths(REFINED_DESIGN, grade=10)
    dimensional = [i for i, n in enumerate(VARIABLE_NAMES) if n not in ANGULAR_VARIABLES]
    assert np.all(loose[dimensional] > tight[dimensional])


def test_angular_variables_are_not_given_a_length_tolerance() -> None:
    """The clocking angles are an assembly quantity, not a machined size."""
    widths = tolerance_half_widths(REFINED_DESIGN, angular=0.05)
    for index, name in enumerate(VARIABLE_NAMES):
        if name in ANGULAR_VARIABLES:
            assert widths[index] == pytest.approx(0.05)


def test_covariance_is_diagonal_and_scales_as_the_square() -> None:
    matrix = covariance(REFINED_DESIGN)
    assert np.allclose(matrix, np.diag(np.diag(matrix)))
    widths = tolerance_half_widths(REFINED_DESIGN)
    assert np.allclose(np.sqrt(np.diag(matrix)), widths / 3.0)


# -- propagation ---------------------------------------------------------------


def test_first_order_and_monte_carlo_agree_within_a_factor_of_two(report: object) -> None:
    """Linearisation is what should be distrusted near a singularity.

    It need not be accurate, but if it were wrong by an order of magnitude the
    exact Jacobians would be giving no useful robustness information at all.
    """
    for name in CONSTRAINT_NAMES:
        sampled = report.monte_carlo_sigma[name]
        if sampled <= 0.0:
            continue
        assert 0.4 < report.linear_sigma[name] / sampled < 2.5, name


def test_first_order_is_conservative_on_this_problem(report: object) -> None:
    """It overestimates sigma rather than under, which is the safe direction.

    A first-order robust formulation built on these gradients would therefore
    be pessimistic, not optimistic -- worth stating, because the opposite would
    make the whole approach unusable here.
    """
    ratios = [
        report.linear_sigma[name] / report.monte_carlo_sigma[name]
        for name in CONSTRAINT_NAMES
        if report.monte_carlo_sigma[name] > 0.0
    ]
    assert np.mean(ratios) > 0.95


def test_the_tdc_gap_constraint_cannot_be_held(report: object) -> None:
    """The headline robustness finding.

    ``g <= 0.01 mm`` is tighter than the variation of the dimensions that
    produce ``g``.  Its capability is far below the industrial target of 1.33,
    and it is violated by a large fraction of nominally-conforming builds.
    This is a defect in the specification, not in any design that meets it.
    """
    assert report.capability()["tdc_gap"] < 0.5
    assert report.violation_rate["tdc_gap"] > 0.25


def test_every_other_constraint_is_comfortably_holdable(report: object) -> None:
    for name in CONSTRAINT_NAMES:
        if name == "tdc_gap":
            continue
        assert report.violation_rate[name] < 0.25, name


def test_the_gap_constraint_is_off_the_bottom_of_the_iso_ladder() -> None:
    """No machining grade fixes it, which is what makes it a specification bug.

    If a tighter grade would work, the answer is to tighten the drawing.  Here
    the required tolerance unit multiple is below the tightest grade in the
    table, so the gap has to be taken up by adjustment at assembly instead.
    """
    grade, factor = required_grade(REFINED_DESIGN, "tdc_gap", crank_samples=CRANK)
    assert grade is None
    assert factor < min(IT_FACTORS.values())


def test_a_slack_constraint_needs_no_precision_at_all() -> None:
    grade, factor = required_grade(REFINED_DESIGN, "compatibility", crank_samples=CRANK)
    assert grade == max(IT_FACTORS)
    assert factor > IT_FACTORS[DEFAULT_GRADE]


def test_report_is_reproducible() -> None:
    first = tolerance_report(COUPLED_DESIGN, samples=60, crank_samples=CRANK, seed=7)
    second = tolerance_report(COUPLED_DESIGN, samples=60, crank_samples=CRANK, seed=7)
    assert first.violation_rate == second.violation_rate


def test_an_unanalysable_design_is_refused_rather_than_scored() -> None:
    """Unlike the optimizer path, this one raises.

    A tolerance study on a design that does not exist would return numbers
    with no meaning, and nothing downstream needs it to keep walking.
    """
    broken = REFINED_DESIGN.replace(a=0.4 * REFINED_DESIGN.a)
    with pytest.raises(ValueError, match="does not close"):
        tolerance_report(broken, samples=5, crank_samples=CRANK)
