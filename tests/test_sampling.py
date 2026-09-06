"""Checking FORM against sampled builds."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from exlink.constants import DEFAULT_TARGETS
from exlink.reference import RELIABLE_DESIGN
from exlink.robustness import RELIABILITY_NAMES, constraint_moments, failure_probability
from exlink.sampling import (
    SampledReliability,
    _design_point,
    _exact_constraints,
    format_sampled,
    sampled_reliability,
)

BAND = {"expansion_stroke": 0.15, "compression_ratio": 0.15}
TARGETS = dataclasses.replace(DEFAULT_TARGETS, max_tdc_gap=0.1)


def _result(system: float, error: float = 0.0, draws: int = 1000):
    return SampledReliability(
        names=RELIABILITY_NAMES,
        per_constraint=dict.fromkeys(RELIABILITY_NAMES, 0.0),
        system=system,
        standard_error=error,
        draws=draws,
        unbuildable=0,
        importance=False,
        effective_draws=float(draws),
    )


def test_the_sampled_constraints_are_the_ones_form_linearises():
    """Any disagreement here would make the whole check meaningless."""
    from exlink.constants import DEFAULT_SPEC

    exact = _exact_constraints(RELIABLE_DESIGN, 360, TARGETS, DEFAULT_SPEC, BAND)
    moments = constraint_moments(RELIABLE_DESIGN, targets=TARGETS, band=BAND)
    assert exact == pytest.approx(moments.value, rel=1.0e-12, abs=1.0e-15)


def test_a_design_that_does_not_close_yields_nothing():
    from exlink.constants import DEFAULT_SPEC

    broken = dataclasses.replace(RELIABLE_DESIGN, a=1.0, c=1.0)
    assert _exact_constraints(broken, 180, TARGETS, DEFAULT_SPEC, BAND) is None


def test_the_design_point_lies_at_beta_standard_deviations():
    beta = np.array([3.0, 8.0])
    jacobian = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    sigma = np.array([2.0, 2.0, 2.0])
    shift = _design_point(beta, jacobian, sigma)
    assert np.linalg.norm(shift) == pytest.approx(3.0)
    # Along the binding constraint's own direction, not the other's.
    assert shift[0] == pytest.approx(3.0)
    assert shift[1] == pytest.approx(0.0)


def test_the_design_point_is_the_origin_when_nothing_is_finite():
    shift = _design_point(np.array([np.inf]), np.array([[1.0, 1.0]]), np.ones(2))
    assert shift == pytest.approx(np.zeros(2))


def test_the_interval_is_clipped_to_a_probability():
    low, high = _result(0.001, error=0.002).interval
    assert low == 0.0
    assert high < 1.0
    assert _result(0.999, error=0.01).interval[1] == 1.0


def test_agreement_is_judged_by_ratio_not_by_difference():
    sampled = _result(8.6e-3)
    assert not sampled.agrees_with(1.3e-3)
    assert sampled.agrees_with(1.0e-2)
    assert _result(0.0, draws=1000).agrees_with(1.0e-4)


def test_the_table_reports_the_comparison():
    text = format_sampled(_result(8.6e-3, error=1.0e-3), form=1.3e-3)
    assert "crude Monte Carlo" in text
    assert "FORM / sampled" in text
    for name in RELIABILITY_NAMES:
        assert name in text


def test_a_short_run_is_wired_end_to_end():
    outcome = sampled_reliability(
        RELIABLE_DESIGN, draws=200, targets=TARGETS, band=BAND, importance=False
    )
    assert outcome is not None
    assert outcome.draws == 200
    assert outcome.effective_draws == pytest.approx(200.0)
    assert 0.0 <= outcome.system <= 1.0


def test_importance_sampling_is_unbiased_where_crude_sampling_can_check_it():
    """Same probability, two densities: the weights have to undo the shift."""
    crude = sampled_reliability(
        RELIABLE_DESIGN, draws=3000, targets=TARGETS, band=BAND, importance=False, seed=1
    )
    shifted = sampled_reliability(
        RELIABLE_DESIGN, draws=3000, targets=TARGETS, band=BAND, importance=True, seed=1
    )
    assert crude is not None
    assert shifted is not None
    assert shifted.importance
    assert shifted.effective_draws < crude.effective_draws
    spread = 3.0 * (crude.standard_error + shifted.standard_error)
    assert abs(shifted.system - crude.system) <= spread


@pytest.mark.slow
def test_form_is_optimistic_at_the_study_result():
    """Section 6.8's finding, at a sample count that can see it."""
    form = failure_probability(RELIABLE_DESIGN, targets=TARGETS, band=BAND)
    outcome = sampled_reliability(
        RELIABLE_DESIGN, draws=20_000, targets=TARGETS, band=BAND, importance=False
    )
    assert outcome is not None
    assert outcome.system > 3.0 * form.system
    assert not outcome.agrees_with(form.system)
