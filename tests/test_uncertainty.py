"""The widened uncertainty model of section 6.7."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from exlink.constants import DEFAULT_TARGETS
from exlink.performance import evaluate
from exlink.reference import RELIABLE_DESIGN
from exlink.robustness import RELIABILITY_NAMES, constraint_moments, reliability_from_moments
from exlink.uncertainty import (
    COUPLED_NAMES,
    DEFAULT_SCATTER,
    PARAMETER_NAMES,
    UNCERTAIN_NAMES,
    WIDENED_NAMES,
    ParameterScatter,
    constraint_vector,
    format_shares,
    widened_moments,
    widened_reliability,
)

BAND = {"expansion_stroke": 0.15, "compression_ratio": 0.15}
TARGETS = dataclasses.replace(DEFAULT_TARGETS, max_tdc_gap=0.1)
RUN = {
    "band": BAND,
    "targets": TARGETS,
    "speed_rpm": 1000.0,
    "module": 0.8,
    "teeth": 48,
}


@pytest.fixture(scope="module")
def widened():
    return widened_moments(RELIABLE_DESIGN, **RUN)


def test_the_uncertain_vector_is_the_eleven_dimensions_plus_six():
    assert len(UNCERTAIN_NAMES) == 17
    assert len(PARAMETER_NAMES) == 6
    assert DEFAULT_SCATTER.as_array().size == 6
    assert UNCERTAIN_NAMES[-6:] == PARAMETER_NAMES


def test_the_widened_model_covers_five_constraints_the_narrow_one_cannot():
    assert WIDENED_NAMES[: len(RELIABILITY_NAMES)] == RELIABILITY_NAMES
    assert set(COUPLED_NAMES).isdisjoint(RELIABILITY_NAMES)
    assert len(COUPLED_NAMES) == 5
    assert len(WIDENED_NAMES) == len(RELIABILITY_NAMES) + 5


def test_the_constraint_vector_agrees_with_the_analytic_one_where_they_overlap():
    """Same constraints, same sign convention, reached two different ways."""
    narrow = constraint_moments(RELIABLE_DESIGN, targets=TARGETS, band=BAND)
    performance = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, module=0.8, teeth=48)
    wide = constraint_vector(performance, TARGETS, BAND)
    assert wide is not None
    assert wide[: narrow.value.size] == pytest.approx(narrow.value, rel=1.0e-6, abs=1.0e-9)


def test_a_design_that_cannot_be_sized_returns_nothing():
    from exlink.reference import PUBLISHED_DESIGN

    performance = evaluate(PUBLISHED_DESIGN, speed_rpm=1000.0)
    if performance.coupled is None:
        assert constraint_vector(performance, TARGETS, BAND) is None


def test_the_purely_geometric_constraints_are_all_dimensions(widened):
    """Nothing a bar is made of moves where its holes are."""
    _moments, shares = widened
    dominant = shares.dominant()
    for name in RELIABILITY_NAMES:
        assert dominant[name] == "dimensions"


def test_each_of_the_five_new_constraints_is_driven_by_a_parameter(widened):
    """If the dimensions drove them, the narrow model would already have them."""
    _moments, shares = widened
    dominant = shares.dominant()
    for name in COUPLED_NAMES:
        assert dominant[name] != "dimensions"


def test_the_variance_split_is_exact_rather_than_attributed(widened):
    _moments, shares = widened
    assert shares.shares.shape == (len(WIDENED_NAMES), 7)
    rows = shares.shares.sum(axis=1)
    assert rows == pytest.approx(np.ones_like(rows))


def test_widening_finds_the_gear_pair_the_narrow_model_could_not_see(widened):
    """Section 6.7's headline: the pinned pair sits on its face-width limit."""
    moments, _shares = widened
    reliability = reliability_from_moments(moments)
    assert reliability.system > 0.4
    worst = moments.names[int(np.argmin(moments.beta))]
    assert worst == "gear"


def test_a_gear_pair_off_its_limit_recovers_the_geometric_answer():
    """With the pair the exhaustive search chose, the geometry binds again."""
    reliability = widened_reliability(
        RELIABLE_DESIGN, band=BAND, targets=TARGETS, speed_rpm=1000.0, module=1.0, teeth=39
    )
    assert reliability is not None
    assert reliability.system < 0.05
    worst = reliability.moments.names[int(np.argmin(reliability.moments.beta))]
    assert worst.startswith("stroke")


def test_the_widened_model_contains_the_narrow_one_exactly(widened):
    """Section 6.7's first finding, and the one that had to be checked first.

    The eight geometric constraints turn out to be untouched by the six new
    parameters, so §5.2's and §5.5's reliability figures stand.  That is a
    result only because it is *exact*: the widened model takes the same
    analytic gradients for these rows rather than differencing them, so any
    disagreement here would be physics rather than a step size.
    """
    moments, _shares = widened
    narrow = constraint_moments(RELIABLE_DESIGN, targets=TARGETS, band=BAND)
    for index, name in enumerate(narrow.names):
        wide = list(moments.names).index(name)
        assert moments.value[wide] == pytest.approx(narrow.value[index], rel=1.0e-9)
        assert moments.sigma[wide] == pytest.approx(narrow.sigma[index], rel=1.0e-9)


def test_scatter_of_zero_leaves_only_the_dimensions(widened):
    quiet = ParameterScatter(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    _moments, shares = widened_moments(RELIABLE_DESIGN, scatter=quiet, **RUN)
    for row, sigma in zip(shares.shares, _moments.sigma, strict=True):
        if sigma > 0.0:
            assert row[0] == pytest.approx(1.0)


def test_the_table_names_every_constraint_and_source(widened):
    _moments, shares = widened
    text = format_shares(shares)
    for name in WIDENED_NAMES:
        assert name in text
    assert "dimension" in text
