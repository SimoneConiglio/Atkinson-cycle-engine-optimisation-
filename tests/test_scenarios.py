"""The GEMSEO problem formulations."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from exlink import (
    DEFAULT_TARGETS,
    GLOBAL_BOUNDS,
    PUBLISHED_DESIGN,
    VARIABLE_NAMES,
    Bounds,
    analyse,
)
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

    # Generations matter more than population size on this problem: at 20
    # generations the run still returns a single point whatever the population,
    # because it has not yet worked its way into the thin feasible region.
    # 35 generations gives fronts of 6-33 points across seeds.
    outcome = local_pareto(
        REFINED_DESIGN,
        relative=0.2,
        absolute_angle=20.0,
        pop_size=80,
        max_gen=35,
        samples=180,
        seed=3,
    )
    assert len(outcome.front) >= 3
    metrics = [analyse(d, samples=180).metrics for d in outcome.front]
    assert all(m.valid for m in metrics)
    # A front, not a cluster: the objectives must actually trade off.
    heights = [m.height for m in metrics]
    assert max(heights) - min(heights) > 5.0


def test_feasibility_tolerates_a_constraint_sitting_on_its_bound() -> None:
    """An augmented Lagrangian converges *onto* the bounds, not inside them.

    A design a few parts in 1e5 outside ``g`` is converged as far as every
    solver in the chain is concerned, so :func:`is_feasible` allows GEMSEO's own
    ``ineq_tolerance``. Ten times that slack is still a violation.
    """
    from exlink.scenarios import INEQUALITY_TOLERANCE

    analysis = analyse(REFINED_DESIGN, samples=1440)
    margin = analysis.metrics.tdc_gap - DEFAULT_TARGETS.max_tdc_gap
    assert margin < 0.0, "the reference should sit strictly inside its bounds"

    # Tighten the bound until the design lands just outside it, by exactly the
    # tolerance, and check both sides of the decision.
    on_bound = replace(
        DEFAULT_TARGETS,
        max_tdc_gap=analysis.metrics.tdc_gap - 0.5 * INEQUALITY_TOLERANCE,
    )
    outside = replace(
        DEFAULT_TARGETS,
        max_tdc_gap=analysis.metrics.tdc_gap - 10.0 * INEQUALITY_TOLERANCE,
    )
    assert is_feasible(analysis, targets=on_bound)
    assert not is_feasible(analysis, targets=outside)


def test_the_coupled_scenario_carries_the_dynamic_constraints() -> None:
    """The MDF problem keeps the report's constraints and adds its own."""
    from exlink.scenarios import COUPLED_INEQUALITY_OUTPUTS, build_coupled_scenario

    problem = build_coupled_scenario(speed_rpm=1000.0).formulation.optimization_problem
    names = set(problem.scalar_constraint_names)
    assert set(INEQUALITY_OUTPUTS) <= names
    assert set(COUPLED_INEQUALITY_OUTPUTS) <= names
    assert problem.objective.name == "total_mass"


def test_the_coupled_scenario_uses_an_mda() -> None:
    """Without one, the coupling would simply be evaluated in the wrong order."""
    from exlink.scenarios import build_coupled_scenario

    scenario = build_coupled_scenario(speed_rpm=1000.0)
    mda = scenario.formulation.mda
    couplings = set(mda.coupling_structure.strong_couplings)
    assert "diameters" in couplings
    assert {"member_axial", "member_bending"} & couplings


def test_the_bearing_margin_is_signed_against_its_limit() -> None:
    from exlink.scenarios import BearingMarginDiscipline

    discipline = BearingMarginDiscipline(limit=10_000.0)
    safe = discipline.execute({"peak_bearing_load": np.array([5_000.0])})
    unsafe = discipline.execute({"peak_bearing_load": np.array([20_000.0])})
    assert float(safe["bearing_margin"][0]) < 0.0
    assert float(unsafe["bearing_margin"][0]) > 0.0


@pytest.mark.slow
def test_minimising_mass_beats_the_quasi_static_design() -> None:
    """The payoff: with inertia in the loop the optimizer finds a far lighter design.

    The report's geometry is near-singular because that maximises the
    quasi-static lever arm. Once the parts have to survive their own inertia
    that becomes the expensive choice, and a short local search should already
    find a much lighter design at a comparable efficiency.
    """
    from exlink.coupled import solve_for_design
    from exlink.scenarios import minimise_mass

    start = solve_for_design(REFINED_DESIGN, speed_rpm=1000.0, samples=360, max_iterations=400)
    outcome = minimise_mass(speed_rpm=1000.0, max_iter=60, relative=0.30, min_efficiency=0.24)
    found = solve_for_design(outcome.design, speed_rpm=1000.0, samples=360, max_iterations=400)
    assert found.feasible
    assert found.total_mass_kg < start.total_mass_kg
