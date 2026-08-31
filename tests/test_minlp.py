"""The gear choice as a mixed-integer problem, via bi-level outer approximation."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.design import VARIABLE_NAMES
from exlink.reference import COUPLED_DESIGN

minlp = pytest.importorskip(
    "exlink.minlp", reason="needs the gemseo-bilevel-outer-approximation plugin"
)


@pytest.fixture(scope="module")
def candidates() -> list:
    return minlp.candidates_from_design(COUPLED_DESIGN, speed_rpm=1000.0, limit=3)


@pytest.fixture(scope="module")
def scenario(candidates: list) -> object:
    return minlp.build_minlp_scenario(
        candidates, COUPLED_DESIGN, speed_rpm=1000.0, sub_max_iter=4
    )


# -- the lattice ---------------------------------------------------------------


def test_a_lattice_point_reproduces_its_own_centre_distance(candidates: list) -> None:
    for point in candidates:
        assert point.inter_axle == pytest.approx(1.5 * point.module * point.teeth)


def test_candidates_are_buildable(candidates: list) -> None:
    """The enumeration must offer pairs that can carry the load, not just near ones."""
    from exlink.gears import size_pair

    assert candidates
    for point in candidates:
        pair = size_pair(point.inter_axle, 1500.0, module=point.module, teeth=point.teeth)
        assert pair.feasible


# -- the bi-level split --------------------------------------------------------


def test_the_formulation_splits_discrete_from_continuous(scenario: object) -> None:
    """The whole point of the Benders formulation: it does the split itself.

    Only the categorical selection is left at the main level; every continuous
    linkage variable belongs to the sub-problem.  ``I`` appears in neither,
    because it is an *output* of the catalogue interpolation rather than a
    design variable at all.
    """
    main = list(scenario.design_space.variable_names)
    assert main == [minlp.GEAR_CHOICE]

    sub = list(scenario.formulation.sub_problem_design_space.variable_names)
    assert set(sub) == {name for name in VARIABLE_NAMES if name != "I"}
    assert "I" not in main and "I" not in sub


def test_the_main_problem_carries_the_feasibility_condition(scenario: object) -> None:
    """How outer approximation handles an infeasible sub-problem.

    Pinning ``I`` to a lattice point throws the design off the equalities that
    ``I`` was one of the variables used to satisfy, so some lattice points have
    no feasible continuous solution.  Attaching the constraints with
    ``main_level=True`` puts an ``is_feasible`` condition on the main problem,
    so such a point is excluded on evidence rather than silently scored.
    """
    names = [c.name for c in scenario.formulation.optimization_problem.constraints]
    assert any("is_feasible" in name for name in names)


def test_the_sub_problem_keeps_the_real_constraints(scenario: object) -> None:
    sub = scenario.formulation.sub_problem_scenario_adapter.scenario
    names = [c.name for c in sub.formulation.optimization_problem.constraints]
    for expected in ("tdc_gap_margin", "side_load_margin", "bearing_margin"):
        assert expected in names


def test_the_catalogue_interpolation_is_exact_at_a_vertex(candidates: list) -> None:
    """At unit penalty the interpolation is ``I = sum_j y_j I_j``.

    A one-hot selection must therefore return exactly that lattice point's
    centre distance, module and tooth count -- otherwise the main problem is
    choosing between mechanisms that do not exist.
    """
    space = minlp.build_minlp_scenario(
        candidates, COUPLED_DESIGN, speed_rpm=1000.0, sub_max_iter=2
    ).formulation.design_space
    del space

    from gemseo_bilevel_outer_approximation.algos.design_space.catalogue_design_space import (
        CatalogueDesignSpace,
    )

    catalogue = CatalogueDesignSpace()
    catalogue.add_categorical_variable(
        name=minlp.GEAR_CHOICE, value=[0], catalogue=list(range(len(candidates)))
    )
    discipline = catalogue.get_catalogue_interpolation_discipline(
        penalty=1.0,
        variable=minlp.GEAR_CHOICE,
        output="I",
        catalogue=np.array([item.inter_axle for item in candidates]),
    )
    for index, point in enumerate(candidates):
        one_hot = np.zeros(len(candidates))
        one_hot[index] = 1.0
        output = discipline.execute({minlp.GEAR_CHOICE: one_hot})
        assert float(np.ravel(output["I"])[0]) == pytest.approx(point.inter_axle)


# -- the discipline change the formulation needs -------------------------------


def test_the_gear_pair_is_an_input_not_construction_data() -> None:
    """The main problem sets the gear pair, so it has to arrive as an input.

    With it fixed at construction time there is nothing for a mixed-integer
    master to choose.
    """
    from exlink.coupled import solve_for_design
    from exlink.disciplines import (
        COUPLING_DIAMETERS,
        GEAR_MODULE,
        GEAR_TEETH,
        RangeDiscipline,
    )
    from exlink.dynamics import MEMBER_NAMES

    sized = solve_for_design(COUPLED_DESIGN, speed_rpm=1000.0)
    diameters = np.array([sized.diameters[name] for name in MEMBER_NAMES])
    discipline = RangeDiscipline(speed_rpm=1000.0)

    def margin(module: float, teeth: int) -> float:
        output = discipline.execute(
            {
                **COUPLED_DESIGN.to_mapping(),
                COUPLING_DIAMETERS: diameters,
                GEAR_MODULE: np.array([module]),
                GEAR_TEETH: np.array([float(teeth)]),
            }
        )
        return float(output["gear_margin"][0])

    # A fine module needs a wide face, so its margin is much smaller.
    assert margin(0.8, 48) < margin(2.0, 19)


def test_both_sign_conventions_are_emitted() -> None:
    """``positive=True`` renames a constraint, which the adapter cannot address.

    So the range margins are published in both conventions and the bi-level
    formulation attaches the violation form.
    """
    from exlink.disciplines import RANGE_OUTPUTS

    for name in ("runs_margin", "gear_margin", "runs_violation", "gear_violation"):
        assert name in RANGE_OUTPUTS


@pytest.mark.slow
def test_the_decomposition_chooses_a_lattice_point(candidates: list) -> None:
    """End to end, on a small budget: it must return a real lattice point."""
    result = minlp.solve(
        COUPLED_DESIGN,
        candidates=candidates,
        speed_rpm=1000.0,
        max_iter=3,
        sub_max_iter=4,
    )
    assert result.point in candidates
    assert result.design is not None
    assert pytest.approx(result.point.inter_axle) == result.design.I
    assert result.iterations >= 1
