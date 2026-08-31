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


def test_the_automatic_module_choice_is_a_step_function() -> None:
    """The reason the gear pair must be pinned during a gradient solve.

    Left free, the module is chosen by a discrete search, so the range is a
    step function of ``I``: it jumps at every threshold where the lightest
    workable module changes.  A finite difference straddling such a threshold
    reports the height of the step divided by the step size, which is
    quantisation noise, not a derivative.  SLSQP handed a subproblem built from
    two such gradients rejects it as "inequality constraints incompatible" and
    terminates without evaluating anything.

    The threshold is located deterministically rather than hoping a reference
    design happens to sit within a micrometre of one -- which is what an
    earlier version of this test did, and it passed or failed depending on the
    machine.
    """
    from exlink.gears import size_pair

    load = 3000.0
    grid = np.linspace(50.0, 70.0, 4001)
    modules = np.array([size_pair(float(value), load).module for value in grid])
    switches = np.flatnonzero(np.diff(modules) != 0.0)
    assert switches.size > 0, "no module-selection threshold in the scanned range"

    index = int(switches[0])
    below, above = float(grid[index]), float(grid[index + 1])
    assert size_pair(below, load).module != size_pair(above, load).module
    # The two sit a fraction of a millimetre apart yet select different hobs.
    assert above - below < 0.02


def test_pinning_the_module_removes_the_step(diameters: np.ndarray) -> None:
    """With the pair pinned, the same neighbourhood is smooth.

    The discipline's own output either side of a module-selection threshold
    must be continuous once the module is fixed, because nothing discrete is
    left in the chain.
    """
    from exlink.gears import size_pair, tooth_count

    load = 3000.0
    grid = np.linspace(50.0, 70.0, 4001)
    modules = np.array([size_pair(float(value), load).module for value in grid])
    index = int(np.flatnonzero(np.diff(modules) != 0.0)[0])
    below, above = float(grid[index]), float(grid[index + 1])

    module = size_pair(below, load).module
    discipline = RangeDiscipline(
        speed_rpm=SPEED, module=module, teeth=tooth_count(below, module)
    )
    left = discipline._evaluate(COUPLED_DESIGN.replace(I=below), diameters)
    right = discipline._evaluate(COUPLED_DESIGN.replace(I=above), diameters)

    # A sub-hundredth-millimetre change in I must not move the range
    # appreciably once the discrete choice is held.
    assert abs(right["km_per_litre"] - left["km_per_litre"]) < 1.0
    assert abs(right["gear_margin"] - left["gear_margin"]) < 0.1


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
    optimizer should reduce ``W`` and gain range -- the whole claim of the
    coupled study, reached without any constraint pointing that way.

    The claim tested is directional, not that the run finishes feasible.  From
    that starting point it does not, within a budget CI can afford: snapping
    ``I`` onto the gear lattice throws the design off the equalities and off
    the top-dead-centre gap, and repairing those while also improving the range
    takes far more than a dozen iterations.  Asserting feasibility here would
    be asserting a convergence rate, which is not the claim and is not
    reproducible across machines.
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

    assert after.analysis.valid
    assert after.km_per_litre > 1.2 * before.km_per_litre
    assert after.metrics.compatibility < before.metrics.compatibility


def test_projection_restores_the_equalities_exactly() -> None:
    """A solver stops within its own tolerance of a constraint, not on it.

    SLSQP finishes a couple of parts in ten thousand outside the relaxed
    equality band here.  The minimum-norm Newton step, built from the exact
    Jacobian rows, puts the design back on the manifold to machine precision
    while changing everything else as little as possible.
    """
    from exlink.model import analyse, equality_constraints
    from exlink.scenarios import project_onto_equalities

    off = COUPLED_DESIGN.replace(e=COUPLED_DESIGN.e * 1.004)
    before = equality_constraints(analyse(off, samples=720))
    assert np.max(np.abs(before)) > 1.0e-3

    projected = project_onto_equalities(off, samples=720)
    after = equality_constraints(analyse(projected, samples=720))
    assert np.max(np.abs(after)) < 1.0e-9


def test_projection_holds_the_gear_lattice() -> None:
    """A projection free to move ``I`` hands back ungearable geometry."""
    from exlink.scenarios import project_onto_equalities

    off = COUPLED_DESIGN.replace(I=57.6, e=COUPLED_DESIGN.e * 1.004)
    projected = project_onto_equalities(off, samples=720)
    assert pytest.approx(57.6, abs=1.0e-12) == projected.I


def test_restoring_the_equalities_can_break_the_gap() -> None:
    """The hypersensitivity of ``g``, seen from a fourth direction.

    A minimum-norm step of a few hundredths of a millimetre -- the smallest
    change that restores ``STE`` and ``epsilon`` -- is enough to move the
    top-dead-centre gap past its own bound.  Together with the tolerance study
    and the gear-lattice snap, this says the same thing three ways: ``g``
    responds to any perturbation of the geometry far faster than its 0.01 mm
    band allows, and no geometric choice can hold it.

    The practical consequence is that a projected design must have its
    inequalities re-checked, never assumed.
    """
    from exlink.model import analyse
    from exlink.scenarios import project_onto_equalities

    # A design sitting on the relaxed equality band with the gap satisfied.
    on_band = COUPLED_DESIGN.replace(e=COUPLED_DESIGN.e * 1.0007)
    start = analyse(on_band, samples=720)
    assert start.valid

    projected = project_onto_equalities(on_band, samples=720)
    landed = analyse(projected, samples=720)
    assert landed.valid
    step = float(np.linalg.norm(projected.to_array() - on_band.to_array()))
    gap_change = abs(landed.metrics.tdc_gap - start.metrics.tdc_gap)
    # A sub-tenth-millimetre step moves the gap by a comparable amount, on a
    # quantity whose entire budget is 0.01 mm.
    assert gap_change / max(step, 1.0e-12) > 0.05
