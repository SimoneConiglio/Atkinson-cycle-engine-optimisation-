"""The GEMSEO problem formulations."""

from __future__ import annotations

import pytest

from exlink import GLOBAL_BOUNDS, PUBLISHED_DESIGN, VARIABLE_NAMES, Bounds, analyse
from exlink.reference import REFINED_DESIGN
from exlink.scenarios import (
    EQUALITY_OUTPUTS,
    INEQUALITY_OUTPUTS,
    build_design_space,
    build_scenario,
    format_analysis,
    is_feasible,
    maximise_efficiency,
    refine,
)


def test_the_design_space_has_one_variable_per_design_variable() -> None:
    space = build_design_space()
    assert list(space.variable_names) == list(VARIABLE_NAMES)
    assert space.dimension == 11


def test_the_design_space_starts_at_the_given_design() -> None:
    space = build_design_space(GLOBAL_BOUNDS, PUBLISHED_DESIGN)
    assert space.get_current_value() == pytest.approx(PUBLISHED_DESIGN.to_array())


def test_the_design_space_respects_the_bounds() -> None:
    box = Bounds.around(PUBLISHED_DESIGN, relative=0.1)
    space = build_design_space(box, PUBLISHED_DESIGN)
    assert space.get_lower_bounds() == pytest.approx(box.lower)
    assert space.get_upper_bounds() == pytest.approx(box.upper)


def test_a_scenario_carries_every_constraint() -> None:
    problem = build_scenario(samples=180).formulation.optimization_problem
    assert set(INEQUALITY_OUTPUTS) == {
        c.name for c in problem.constraints.get_inequality_constraints()
    }
    assert set(EQUALITY_OUTPUTS) == {
        c.name for c in problem.constraints.get_equality_constraints()
    }


def test_relaxing_equalities_replaces_them_with_bracketing_inequalities() -> None:
    """NSGA-II takes no equality constraints, so they become |residual| <= tol."""
    problem = build_scenario(
        samples=180, relax_equalities=True
    ).formulation.optimization_problem
    names = set(problem.scalar_constraint_names)
    for base in EQUALITY_OUTPUTS:
        assert f"{base}_upper" in names
        assert f"{base}_lower" in names
    assert not list(problem.constraints.get_equality_constraints())


def test_moving_limits_become_constraints() -> None:
    problem = build_scenario(
        samples=180, max_height=250.0, max_width=160.0
    ).formulation.optimization_problem
    assert {"height_limit", "width_limit"} <= set(problem.scalar_constraint_names)


def test_a_multi_objective_scenario_carries_three_objectives() -> None:
    problem = build_scenario(
        ["neg_efficiency", "height", "width"], samples=180
    ).formulation.optimization_problem
    assert problem.objective.output_names == ["neg_efficiency", "height", "width"]
    assert not problem.is_mono_objective


def test_is_feasible_accepts_the_refined_reference() -> None:
    assert is_feasible(analyse(REFINED_DESIGN, samples=1440))


def test_is_feasible_rejects_the_published_table() -> None:
    """Re-analysed as printed, the published design misses several constraints."""
    assert not is_feasible(analyse(PUBLISHED_DESIGN, samples=1440))


def test_is_feasible_rejects_a_penalised_design() -> None:
    assert not is_feasible(analyse(PUBLISHED_DESIGN.replace(a=25.0, c=25.0), samples=180))


def test_format_analysis_reports_every_constraint(refined_analysis) -> None:
    text = format_analysis(refined_analysis)
    for symbol in ("eta", "STE", "eps", "mra", "W", "g", "d", "gamma", "H", "B"):
        assert symbol in text
    assert "feasible: True" in text


def test_format_analysis_explains_a_penalised_design() -> None:
    analysis = analyse(PUBLISHED_DESIGN.replace(a=25.0, c=25.0), samples=180)
    assert "PENALISED" in format_analysis(analysis)


@pytest.mark.slow
def test_the_augmented_lagrangian_makes_the_published_design_feasible() -> None:
    """The report's own last step, reproduced.

    Started from the published table -- which does not satisfy its own
    constraints when re-analysed -- the augmented Lagrangian must land on a
    genuinely feasible design without giving up much efficiency.
    """
    outcome = refine(
        PUBLISHED_DESIGN,
        samples=360,
        max_iter=60,
        relative=0.25,
        sub_algorithm_settings={"max_iter": 200},
    )
    assert outcome.feasible
    assert outcome.analysis.metrics.efficiency > 0.20
    assert outcome.n_evaluations > 0


@pytest.mark.slow
def test_a_local_search_improves_on_its_starting_point() -> None:
    """A short local run must not return something worse than it started from."""
    outcome = maximise_efficiency(
        algorithm="NLOPT_COBYLA",
        bounds=Bounds.around(REFINED_DESIGN, relative=0.05),
        initial=REFINED_DESIGN,
        max_iter=120,
        samples=360,
    )
    start = analyse(REFINED_DESIGN, samples=360).metrics.efficiency
    assert outcome.analysis.metrics.efficiency >= start - 1e-3


def test_moea_targets_relax_only_the_tdc_gap() -> None:
    """``g`` is the constraint that makes a population-based run hopeless."""
    from exlink import DEFAULT_TARGETS
    from exlink.scenarios import MOEA_TDC_GAP, moea_targets

    relaxed = moea_targets()
    assert relaxed.max_tdc_gap == MOEA_TDC_GAP
    assert relaxed.max_tdc_gap > DEFAULT_TARGETS.max_tdc_gap
    for field in (
        "expansion_stroke",
        "compression_ratio",
        "max_rod_angle",
        "max_transmission",
        "min_clearance",
        "max_side_load",
    ):
        assert getattr(relaxed, field) == getattr(DEFAULT_TARGETS, field)


def test_the_tdc_gap_is_effectively_an_equality_constraint() -> None:
    """Document the measurement behind :func:`moea_targets`.

    Sampled over a box around a good design, ``g <= 0.01 mm`` is satisfied by a
    fraction of a percent of designs -- orders of magnitude rarer than any other
    inequality. That is why the multi-objective stage runs relaxed.
    """
    import numpy as np

    from exlink import Design

    rng = np.random.default_rng(0)
    box = Bounds.around(REFINED_DESIGN, relative=0.15, absolute_angle=15.0)
    samples = rng.uniform(box.lower, box.upper, size=(250, 11))

    analysed = [analyse(Design.from_array(x), samples=180) for x in samples]
    valid = [a.metrics for a in analysed if a.valid]
    assert valid, "expected some analysable designs in the box"

    tight = sum(m.tdc_gap <= 0.01 for m in valid) / len(valid)
    loose = sum(m.clearance >= 10.0 for m in valid) / len(valid)
    assert tight < 0.05
    assert loose > tight


@pytest.mark.slow
@pytest.mark.moea
def test_nsga2_returns_a_usable_front() -> None:
    """The relaxed multi-objective run must return a front, not a single point."""
    pymoo = pytest.importorskip("gemseo_pymoo")
    del pymoo

    from exlink.scenarios import local_pareto

    outcome = local_pareto(
        REFINED_DESIGN,
        relative=0.2,
        absolute_angle=20.0,
        pop_size=60,
        max_gen=20,
        samples=180,
        seed=3,
    )
    assert len(outcome.front) > 1
    metrics = [analyse(d, samples=180).metrics for d in outcome.front]
    assert any(m.valid for m in metrics)
    # A front, not a cluster: the objectives must actually trade off.
    heights = [m.height for m in metrics if m.valid]
    assert max(heights) - min(heights) > 1.0
