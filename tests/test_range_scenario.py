"""The range problem: its grammar, its mixed-integer structure, and its gradients."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.design import VARIABLE_NAMES, Bounds
from exlink.disciplines import COUPLING_DIAMETERS, RANGE_OUTPUTS, RangeDiscipline
from exlink.dynamics import MEMBER_NAMES
from exlink.gears import lattice_inter_axle, size_pair, tooth_count
from exlink.reference import COUPLED_DESIGN
from exlink.scenarios import (
    RANGE_INEQUALITY_OUTPUTS,
    _best_design,
    build_design_space,
    build_range_scenario,
)

SPEED = 1000.0


@pytest.fixture(scope="module")
def diameters() -> np.ndarray:
    from exlink.coupled import solve_for_design

    result = solve_for_design(COUPLED_DESIGN, speed_rpm=SPEED)
    return np.array([result.diameters[name] for name in MEMBER_NAMES])


# -- the discipline -----------------------------------------------------------


def test_range_discipline_reports_every_declared_output(diameters: np.ndarray) -> None:
    discipline = RangeDiscipline(speed_rpm=SPEED)
    output = discipline.execute({**COUPLED_DESIGN.to_mapping(), COUPLING_DIAMETERS: diameters})
    for name in RANGE_OUTPUTS:
        assert name in output
        assert np.isfinite(float(output[name][0])), name


def test_an_unanalysable_design_returns_the_penalty_row() -> None:
    """The optimizer has to be able to walk through infeasible ground."""
    discipline = RangeDiscipline(speed_rpm=SPEED)
    broken = COUPLED_DESIGN.replace(a=0.4 * COUPLED_DESIGN.a)
    output = discipline.execute(
        {**broken.to_mapping(), COUPLING_DIAMETERS: np.full(len(MEMBER_NAMES), 8.0)}
    )
    assert float(output["km_per_litre"][0]) == 0.0
    assert float(output["runs_margin"][0]) < 0.0


def test_pinning_the_module_makes_the_gradient_meaningful(diameters: np.ndarray) -> None:
    """The bug that stopped the optimizer dead, as a regression test.

    Left to choose its own module, the discipline's objective is a step
    function of ``I``: the lightest workable module changes at a threshold and
    the range jumps across it.  A central difference straddling that threshold
    returns a gradient two orders of magnitude too large, and SLSQP rejects the
    resulting subproblem as "inequality constraints incompatible" without
    evaluating anything.

    With the pair pinned the same derivative is finite and sane.  The test is
    the ratio, not the absolute value: the floating-module gradient is
    quantisation noise and its magnitude depends on the step size, but it must
    be enormously larger than the pinned one.
    """
    from exlink.coupled import solve_for_design
    from exlink.reference import REFINED_DESIGN

    # The pathology needs a design sitting near a module-selection threshold,
    # which the near-singular reference does and the coupled one does not.
    sized = solve_for_design(REFINED_DESIGN, speed_rpm=SPEED)
    data = {
        **REFINED_DESIGN.to_mapping(),
        COUPLING_DIAMETERS: np.array([sized.diameters[n] for n in MEMBER_NAMES]),
    }

    floating = RangeDiscipline(speed_rpm=SPEED)
    floating.linearize(data, compute_all_jacobians=True)
    noisy = abs(float(floating.jac["neg_range"]["I"][0, 0]))

    module = size_pair(REFINED_DESIGN.I, 1000.0).module
    pinned_discipline = RangeDiscipline(
        speed_rpm=SPEED, module=module, teeth=tooth_count(REFINED_DESIGN.I, module)
    )
    pinned_discipline.linearize(data, compute_all_jacobians=True)
    clean = abs(float(pinned_discipline.jac["neg_range"]["I"][0, 0]))

    assert clean < 1.0e4
    assert noisy > 100.0 * clean


def test_pinned_gradients_agree_with_a_coarse_difference(diameters: np.ndarray) -> None:
    """The finite-difference Jacobian has to be the derivative of the run.

    Cheap to get wrong -- a mismatched perturbation index or a stale cached
    state would go unnoticed -- and everything downstream depends on it.
    """
    module = size_pair(COUPLED_DESIGN.I, 1000.0).module
    teeth = tooth_count(COUPLED_DESIGN.I, module)
    discipline = RangeDiscipline(speed_rpm=SPEED, module=module, teeth=teeth)
    data = {**COUPLED_DESIGN.to_mapping(), COUPLING_DIAMETERS: diameters}
    discipline.linearize(data, compute_all_jacobians=True)

    # A step inside the converged plateau.  The derivative is flat from about
    # 1e-4 to 1e-6 relative; at 1e-3 the chain is genuinely nonlinear over the
    # interval and a difference there disagrees by a factor of two -- which
    # says nothing about the Jacobian and everything about the step.
    step = 1.0e-4 * abs(COUPLED_DESIGN.a)
    forward = discipline._evaluate(
        COUPLED_DESIGN.replace(a=COUPLED_DESIGN.a + step), diameters
    )["km_per_litre"]
    backward = discipline._evaluate(
        COUPLED_DESIGN.replace(a=COUPLED_DESIGN.a - step), diameters
    )["km_per_litre"]
    coarse = -(forward - backward) / (2.0 * step)
    analytic = float(discipline.jac["neg_range"]["a"][0, 0])
    assert coarse == pytest.approx(analytic, rel=0.02)


# -- the scenario -------------------------------------------------------------


def test_choosing_the_gear_pair_removes_i_from_the_design_space() -> None:
    """``I = 1.5 m z``, so once the integers are chosen ``I`` is not free.

    Ten continuous variables, not eleven.  That is the mixed-integer structure
    of the problem made concrete.
    """
    scenario = build_range_scenario(
        bounds=Bounds.around(COUPLED_DESIGN, relative=0.3),
        initial=COUPLED_DESIGN,
        speed_rpm=SPEED,
    )
    names = list(scenario.design_space.variable_names)
    assert len(names) == len(VARIABLE_NAMES) - 1
    assert "I" not in names


def test_the_scenario_pins_i_to_the_lattice() -> None:
    scenario = build_range_scenario(
        bounds=Bounds.around(COUPLED_DESIGN, relative=0.3),
        initial=COUPLED_DESIGN,
        speed_rpm=SPEED,
        module=1.5,
        teeth=25,
    )
    for discipline in scenario.disciplines:
        if "I" in discipline.input_grammar:
            value = float(np.ravel(discipline.default_input_data["I"])[0])
            assert value == pytest.approx(lattice_inter_axle(1.5, 25))


def test_best_design_recovers_the_pinned_variable() -> None:
    """A ten-variable solution has to come back as an eleven-variable design.

    And the missing entry must be the *pinned* value the disciplines actually
    ran with, not the starting design's original one -- otherwise every
    reported result would describe a mechanism that was never evaluated.
    """
    scenario = build_range_scenario(
        bounds=Bounds.around(COUPLED_DESIGN, relative=0.3),
        initial=COUPLED_DESIGN,
        speed_rpm=SPEED,
        module=1.5,
        teeth=25,
    )
    scenario.execute(algo_name="SLSQP", max_iter=2)
    design = _best_design(scenario, objective="neg_range")
    assert pytest.approx(lattice_inter_axle(1.5, 25)) == design.I
    assert pytest.approx(COUPLED_DESIGN.I) != design.I


def test_design_space_can_fix_any_variable() -> None:
    space = build_design_space(
        Bounds.around(COUPLED_DESIGN, relative=0.2), COUPLED_DESIGN, fixed=("I", "e")
    )
    names = list(space.variable_names)
    assert "I" not in names and "e" not in names
    assert len(names) == len(VARIABLE_NAMES) - 2


def test_the_range_margins_are_attached_the_right_way_round() -> None:
    """Two opposite sign conventions meet in this scenario.

    The coupled margins are violations bounded above by zero; the range margins
    are margins bounded below by zero.  Attaching them the same way makes the
    starting point look infeasible and the optimizer stops immediately, which
    is exactly what happened.
    """
    scenario = build_range_scenario(
        bounds=Bounds.around(COUPLED_DESIGN, relative=0.3),
        initial=COUPLED_DESIGN,
        speed_rpm=SPEED,
    )
    problem = scenario.formulation.optimization_problem

    # GEMSEO negates a constraint declared positive, so the names record which
    # convention each one was attached under.  Exactly the range margins should
    # be negated, and none of the coupled ones.
    negated = {c.name for c in problem.constraints if c.name.startswith("-")}
    assert negated == {f"-{name}" for name in RANGE_INEQUALITY_OUTPUTS}

    # And the range margins must be satisfied at the start, which is what fails
    # when they are attached under the coupled convention.
    start = problem.design_space.get_current_value()
    for constraint in problem.constraints:
        if constraint.name not in negated:
            continue
        value = np.atleast_1d(constraint.evaluate(start))
        assert np.all(value <= 1.0e-6), f"{constraint.name} = {value}"


def test_snapping_i_to_the_lattice_breaks_the_other_constraints() -> None:
    """The mixed-integer structure, made visible.

    ``I`` is one of the variables the equalities ``STE = 74`` and
    ``epsilon = 16`` are satisfied with, so moving it onto the gear lattice --
    0.3 mm here -- pushes the design off constraints it previously met.  The
    remaining ten variables have to repair them, which is why each candidate
    gear pair defines its own continuous problem rather than a mere re-scoring
    of one solution.
    """
    from exlink.performance import evaluate

    module = size_pair(COUPLED_DESIGN.I, 1000.0).module
    pinned = COUPLED_DESIGN.replace(
        I=lattice_inter_axle(module, tooth_count(COUPLED_DESIGN.I, module))
    )
    assert pytest.approx(COUPLED_DESIGN.I) != pinned.I

    before = evaluate(COUPLED_DESIGN, speed_rpm=SPEED)
    after = evaluate(pinned, speed_rpm=SPEED, module=module)
    assert before.analysis.valid and after.analysis.valid
    # Some constraint moved measurably: the snap is not a rounding detail.
    assert abs(after.metrics.rod_angle - before.metrics.rod_angle) > 1.0e-3


@pytest.mark.slow
def test_the_optimizer_moves_away_from_the_singularity() -> None:
    """Range-driven optimization has to *discover* the retreat, not be told it.

    Started at the near-singular geometry the quasi-static problem prefers, the
    optimizer should reduce ``W`` and gain range -- which is the whole claim of
    the coupled study, arrived at without any constraint pointing that way.
    """
    from exlink.performance import evaluate
    from exlink.reference import REFINED_DESIGN

    before = evaluate(REFINED_DESIGN, speed_rpm=SPEED)
    scenario = build_range_scenario(
        bounds=Bounds.around(REFINED_DESIGN, relative=0.45),
        initial=REFINED_DESIGN,
        speed_rpm=SPEED,
    )
    scenario.execute(algo_name="SLSQP", max_iter=12)
    after = evaluate(_best_design(scenario, objective="neg_range"), speed_rpm=SPEED)

    assert after.feasible
    assert after.km_per_litre > 1.2 * before.km_per_litre
    assert after.metrics.compatibility < before.metrics.compatibility
