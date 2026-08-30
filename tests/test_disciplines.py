"""The GEMSEO discipline wrappers."""

from __future__ import annotations

import numpy as np
import pytest

from exlink import PUBLISHED_DESIGN, VARIABLE_NAMES, analyse
from exlink.disciplines import OUTPUT_NAMES, ExlinkDiscipline, PenalisedExlinkDiscipline


@pytest.fixture(scope="module")
def discipline() -> ExlinkDiscipline:
    return ExlinkDiscipline(samples=360)


def test_the_grammars_declare_the_expected_names(discipline) -> None:
    assert set(discipline.input_grammar.names) == set(VARIABLE_NAMES)
    assert set(OUTPUT_NAMES) <= set(discipline.output_grammar.names)


def test_defaults_are_the_published_design(discipline) -> None:
    for name in VARIABLE_NAMES:
        assert discipline.default_input_data[name][0] == pytest.approx(
            getattr(PUBLISHED_DESIGN, name)
        )


def test_every_output_is_a_finite_scalar(discipline) -> None:
    output = discipline.execute()
    for name in OUTPUT_NAMES:
        value = np.ravel(output[name])
        assert value.shape == (1,), name
        assert np.isfinite(value[0]), name


def test_outputs_agree_with_a_direct_analysis(discipline) -> None:
    output = discipline.execute()
    metrics = analyse(PUBLISHED_DESIGN, samples=360).metrics
    assert output["efficiency"][0] == pytest.approx(metrics.efficiency)
    assert output["neg_efficiency"][0] == pytest.approx(-metrics.efficiency)
    assert output["height"][0] == pytest.approx(metrics.height)
    assert output["width"][0] == pytest.approx(metrics.width)


def test_an_unanalysable_design_returns_penalised_outputs(discipline) -> None:
    broken = PUBLISHED_DESIGN.replace(a=25.0, c=25.0)
    output = discipline.execute(broken.to_mapping())
    assert output["valid"][0] == 0.0
    assert output["height"][0] == pytest.approx(1000.0)
    assert np.isfinite(output["side_load_ratio"][0])


def test_the_penalised_objective_rewards_feasibility() -> None:
    """``F(X)`` must be worse for a design that violates the constraints."""
    from exlink.reference import REFINED_DESIGN

    discipline = PenalisedExlinkDiscipline(samples=360, penalty_parameter=0.5)
    good = discipline.execute(REFINED_DESIGN.to_mapping())["penalised_objective"][0]
    bad = discipline.execute(PUBLISHED_DESIGN.to_mapping())["penalised_objective"][0]
    assert good < bad


def test_the_penalty_parameter_is_validated() -> None:
    with pytest.raises(ValueError, match=r"must lie in \(0, 1\]"):
        PenalisedExlinkDiscipline(penalty_parameter=0.0)


def test_a_smaller_penalty_parameter_bites_harder() -> None:
    infeasible = PUBLISHED_DESIGN.to_mapping()
    loose = PenalisedExlinkDiscipline(samples=360, penalty_parameter=0.9)
    tight = PenalisedExlinkDiscipline(samples=360, penalty_parameter=0.2)
    assert (
        tight.execute(infeasible)["penalised_objective"][0]
        > loose.execute(infeasible)["penalised_objective"][0]
    )


def test_moving_limits_enter_the_penalty() -> None:
    unlimited = PenalisedExlinkDiscipline(samples=360, penalty_parameter=0.5)
    limited = PenalisedExlinkDiscipline(samples=360, penalty_parameter=0.5, max_height=150.0)
    data = PUBLISHED_DESIGN.to_mapping()
    assert (
        limited.execute(data)["penalised_objective"][0]
        > unlimited.execute(data)["penalised_objective"][0]
    )
