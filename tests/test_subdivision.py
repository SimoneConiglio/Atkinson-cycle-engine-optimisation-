"""Applying the box subdivision here: the box start, and the wiring it needs."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.design import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design
from exlink.reference import REFINED_DESIGN
from exlink.scenarios import (
    EQUALITY_OUTPUTS,
    INEQUALITY_OUTPUTS,
    analyse,
    is_feasible,
)
from exlink.subdivision import (
    CONSTRAINT_SCALES,
    MARGIN_FLOOR,
    BoxStart,
    _decode_one_hot,
    _relaxed_equality_disciplines,
    find_box_start,
    format_subdivision,
    scaled_margins,
)

pytest.importorskip("gemseo_box_subdivision")

GRID = {"q_1": 2, "q_2": 2, "theta_f": 3}
SAMPLES = 180
INDEX = {name: index for index, name in enumerate(VARIABLE_NAMES)}


def box_of(combination: tuple[int, ...], grid: dict[str, int] = GRID) -> Bounds:
    """The corner pair of one box of the grid, over all eleven variables."""
    lower = np.array(GLOBAL_BOUNDS.lower, dtype=float)
    upper = np.array(GLOBAL_BOUNDS.upper, dtype=float)
    for name, position in zip(grid, combination, strict=True):
        index = INDEX[name]
        width = (GLOBAL_BOUNDS.upper[index] - GLOBAL_BOUNDS.lower[index]) / grid[name]
        edge = float(GLOBAL_BOUNDS.lower[index])
        lower[index] = edge + position * width
        upper[index] = edge + (position + 1) * width
    return Bounds(lower=lower, upper=upper)


def test_the_margins_agree_with_the_feasibility_test() -> None:
    """A positive worst margin must mean exactly what ``is_feasible`` means."""
    margins = scaled_margins(REFINED_DESIGN, samples=SAMPLES)
    assert margins.size == len(INEQUALITY_OUTPUTS) + 2 * len(EQUALITY_OUTPUTS)
    assert margins.min() > 0.0
    assert is_feasible(analyse(REFINED_DESIGN, samples=SAMPLES))


def test_the_margins_are_floored_on_the_plateau() -> None:
    """An unanalysable design must not report a margin that swamps every other.

    Its ``tdc_gap`` is the penalty, 1000 against a bound of 0.01, so the
    unclamped margin is -1e5 and the worst-margin ranking becomes a ranking of
    one constraint.
    """
    centre = Design.from_array(0.5 * (GLOBAL_BOUNDS.lower + GLOBAL_BOUNDS.upper))
    assert not analyse(centre, samples=SAMPLES).valid
    margins = scaled_margins(centre, samples=SAMPLES)
    assert margins.min() == MARGIN_FLOOR
    assert set(CONSTRAINT_SCALES) == set(INEQUALITY_OUTPUTS)


def test_the_margins_are_scaled_by_their_own_bound() -> None:
    """Each row must be the violation as a fraction of what that bound is."""
    margins = scaled_margins(REFINED_DESIGN, samples=SAMPLES)
    analysis = analyse(REFINED_DESIGN, samples=SAMPLES)
    from exlink.model import inequality_constraints

    for position, (name, value) in enumerate(
        zip(INEQUALITY_OUTPUTS, inequality_constraints(analysis), strict=True)
    ):
        assert margins[position] == pytest.approx(-value / CONSTRAINT_SCALES[name])


def test_a_box_start_is_analysable_where_the_centre_is_not() -> None:
    """The point of the policy: the centre of a box is usually not evaluable."""
    box = box_of((1, 0, 1))
    centre = Design.from_array(0.5 * (box.lower + box.upper))
    assert not analyse(centre, samples=SAMPLES).valid

    start = find_box_start(
        box.lower, box.upper, incumbent=REFINED_DESIGN, seed=0, probes=128,
        restarts=1, samples=SAMPLES,
    )
    assert isinstance(start, BoxStart)
    assert analyse(start.design, samples=SAMPLES).valid
    assert start.stage in {"probe", "descent", "restoration"}


def test_a_box_start_stays_inside_its_box() -> None:
    """A start outside the box is a start for a different sub-problem."""
    box = box_of((1, 1, 0))
    start = find_box_start(
        box.lower, box.upper, incumbent=REFINED_DESIGN, seed=1, probes=64,
        restarts=1, samples=SAMPLES,
    )
    vector = start.design.to_array()
    assert np.all(vector >= box.lower - 1e-9)
    assert np.all(vector <= box.upper + 1e-9)


def test_a_box_start_is_reproducible() -> None:
    """A box the master re-proposes must be the same sub-problem as before."""
    box = box_of((0, 0, 2))
    kwargs = {
        "incumbent": REFINED_DESIGN, "seed": 3, "probes": 64,
        "restarts": 1, "samples": SAMPLES,
    }
    first = find_box_start(box.lower, box.upper, **kwargs)
    second = find_box_start(box.lower, box.upper, **kwargs)
    assert first.design.to_array() == pytest.approx(second.design.to_array())
    assert first.margin == pytest.approx(second.margin)


def test_the_incumbent_is_used_where_it_lies_in_the_box() -> None:
    """The box holding a feasible design must cost one analysis, not a probe."""
    box = box_of((0, 0, 2))
    assert box.contains(REFINED_DESIGN)
    start = find_box_start(
        box.lower, box.upper, incumbent=REFINED_DESIGN, seed=0, probes=512,
        restarts=1, samples=SAMPLES,
    )
    assert start.stage == "probe"
    assert start.feasible
    assert start.evaluations == 1


def test_a_restoration_beats_the_probe_it_started_from() -> None:
    """Restoration must improve the margin, or the stage is not worth its cost."""
    box = box_of((1, 0, 2))
    probed = find_box_start(
        box.lower, box.upper, incumbent=REFINED_DESIGN, seed=0, probes=64,
        restarts=0, samples=SAMPLES,
    )
    restored = find_box_start(
        box.lower, box.upper, incumbent=REFINED_DESIGN, seed=0, probes=64,
        restarts=2, samples=SAMPLES,
    )
    assert restored.margin >= probed.margin


def test_the_relaxed_equalities_become_their_own_outputs() -> None:
    """They cannot be renamed constraints; the Benders adapter rejects those.

    ``MDOScenarioAdapterBenders`` differentiates a sub-problem constraint by
    looking its name up among the discipline outputs, so a constraint named
    anything but an output raises ``KeyError`` on the first linearisation.
    """
    disciplines, names = _relaxed_equality_disciplines()
    assert names == (
        "stroke_error_upper",
        "stroke_error_lower",
        "compression_ratio_error_upper",
        "compression_ratio_error_lower",
    )
    assert len(disciplines) == len(names)
    for discipline, name in zip(disciplines, names, strict=True):
        assert name in discipline.output_grammar


def test_the_two_sides_bracket_the_band() -> None:
    """Both sides must be non-positive exactly on the band, and only there."""
    from gemseo.core.chains.chain import MDOChain

    from exlink.disciplines import ExlinkDiscipline
    from exlink.scenarios import DEFAULT_EQUALITY_TOLERANCE

    disciplines, names = _relaxed_equality_disciplines()
    chain = MDOChain([ExlinkDiscipline(samples=SAMPLES), *disciplines])
    output = chain.execute(REFINED_DESIGN.to_mapping())
    residual = float(np.ravel(output["stroke_error"])[0])
    half_width = DEFAULT_EQUALITY_TOLERANCE["stroke_error"]
    assert float(np.ravel(output["stroke_error_upper"])[0]) == pytest.approx(
        residual - half_width
    )
    assert float(np.ravel(output["stroke_error_lower"])[0]) == pytest.approx(
        -residual - half_width
    )
    assert abs(residual) <= half_width
    assert all(float(np.ravel(output[name])[0]) <= 0.0 for name in names)


def test_the_sides_carry_the_gradient_of_the_residual() -> None:
    """A side whose Jacobian were zero would let the solver ignore the band."""
    from gemseo.core.chains.chain import MDOChain

    from exlink.disciplines import ExlinkDiscipline

    disciplines, _ = _relaxed_equality_disciplines()
    chain = MDOChain([ExlinkDiscipline(samples=SAMPLES), *disciplines])
    data = REFINED_DESIGN.to_mapping()
    chain.execute(data)
    jacobian = chain.linearize(data, compute_all_jacobians=True)

    def row(name: str) -> np.ndarray:
        return np.array(
            [float(np.ravel(jacobian[name][v])[0]) for v in VARIABLE_NAMES]
        )

    residual = row("stroke_error")
    assert row("stroke_error_upper") == pytest.approx(residual)
    assert row("stroke_error_lower") == pytest.approx(-residual)


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        ((1, 0, 1, 0, 0, 0, 1), (0, 0, 2)),
        ((0, 1, 1, 0, 1, 0, 0), (1, 0, 0)),
        ((0, 1, 0, 1, 0, 1, 0), (1, 1, 1)),
    ],
)
def test_the_one_hot_decodes_to_a_box_index(key, expected) -> None:
    """The report is read box by box, so the labels have to be the right ones."""
    assert _decode_one_hot(key, GRID) == expected


def test_the_scenario_puts_every_constraint_on_the_master() -> None:
    """A box holding no feasible design must produce a feasibility cut.

    Without it the master waits on a sub-problem that cannot converge, and
    most boxes of this design space hold nothing.
    """
    from gemseo_box_subdivision import BoxSubdivisionSettings

    from exlink.subdivision import build_subdivided_scenario

    scenario = build_subdivided_scenario(
        GRID,
        samples=SAMPLES,
        settings=BoxSubdivisionSettings(max_iter=1, sub_problem_max_iter=1),
        probes=4,
        restarts=0,
    )
    assert scenario.formulation.optimization_problem.constraints.get_names() == [
        "[is_feasible-1.0]"
    ]
    sub_problem = (
        scenario.formulation.sub_problem_scenario_adapter.scenario.formulation
    ).optimization_problem
    names = sub_problem.constraints.get_names()
    for name in INEQUALITY_OUTPUTS:
        assert name in names
    assert "stroke_error_upper" in names
    assert "compression_ratio_error_lower" in names


def test_the_subdivided_space_normalizes_only_what_it_subdivides() -> None:
    """A variable left out must stay itself, in millimetres, with its own bounds."""
    from exlink.subdivision import build_subdivided_scenario

    scenario = build_subdivided_scenario(GRID, samples=SAMPLES, restoring=False)
    space = scenario.formulation.sub_problem_scenario_adapter.scenario.design_space
    for name in GRID:
        assert f"{name}_normalized" in space
        assert name not in space
    for name in ("a", "c", "I", "theta_r"):
        assert name in space
    assert scenario.subdivision.variable_names == tuple(GRID)


def test_the_adapter_starts_the_sub_problem_where_it_restored_to() -> None:
    """The whole point: the sub-problem must begin at an analysable point.

    The normalized formulation solves for ``u`` in [0, 1] within the box, so
    the check undoes that mapping and compares in millimetres.
    """
    from gemseo_box_subdivision import BoxSubdivisionSettings

    from exlink.subdivision import build_subdivided_scenario

    scenario = build_subdivided_scenario(
        GRID,
        samples=SAMPLES,
        settings=BoxSubdivisionSettings(max_iter=1, sub_problem_max_iter=2),
        probes=64,
        restarts=1,
    )
    adapter = scenario.formulation.sub_problem_scenario_adapter
    scenario.execute()

    starts = type(adapter).starts
    assert starts, "the adapter must have been asked for at least one box"
    for start in starts.values():
        vector = start.design.to_array()
        assert np.all(vector >= GLOBAL_BOUNDS.lower - 1e-9)
        assert np.all(vector <= GLOBAL_BOUNDS.upper + 1e-9)


def test_the_default_policy_finds_nothing_here() -> None:
    """The measurement the restoration exists for.

    Left on the center of the box, every sub-problem starts on the penalty
    plateau, where the objective and five of the seven constraints are flat
    with a zero gradient, and SLSQP returns its starting point.
    """
    from gemseo_box_subdivision import BoxSubdivisionSettings

    from exlink.subdivision import maximise_efficiency_by_subdivision

    outcome = maximise_efficiency_by_subdivision(
        GRID,
        samples=SAMPLES,
        restoring=False,
        settings=BoxSubdivisionSettings(max_iter=2, sub_problem_max_iter=10),
    )
    assert outcome.design is None
    assert outcome.feasible_boxes == 0


def test_the_report_says_so_when_nothing_was_found() -> None:
    """A run that found nothing must not print a design it does not have."""
    from exlink.subdivision import SubdivisionOutcome

    empty = SubdivisionOutcome(
        design=None,
        efficiency=float("nan"),
        evaluations=65,
        restoration_evaluations=0,
        boxes=0,
        feasible_boxes=0,
    )
    text = format_subdivision(empty)
    assert "no feasible design found" in text
    assert "efficiency" not in text.split("no feasible")[0].split("analyses")[-1]


def test_the_range_problem_refuses_to_subdivide_the_pinned_variable() -> None:
    """``I`` is not a design variable of the range problem; the gear pair owns it."""
    from exlink.subdivision import build_subdivided_range_scenario

    with pytest.raises(ValueError, match="pinned by the gear pair"):
        build_subdivided_range_scenario({"I": 2})


def test_the_range_problem_keeps_its_coupling() -> None:
    """Chaining the five disciplines instead of converging them is not MDF.

    The diameters that leave the structure discipline have to reach the
    dynamics one again, or the sub-problem optimises an engine whose inertia
    and whose sections disagree.
    """
    from gemseo.mda.base_mda import BaseMDA

    from exlink.subdivision import build_subdivided_range_scenario

    scenario = build_subdivided_range_scenario(
        {"q_1": 2, "theta_f": 2}, samples=SAMPLES, restoring=False
    )
    chain = scenario.disciplines[0]
    assert any(isinstance(item, BaseMDA) for item in chain.disciplines)

    space = scenario.formulation.sub_problem_scenario_adapter.scenario.design_space
    assert "I" not in space
    assert "I_normalized" not in space


@pytest.mark.parametrize("problem", ["efficiency", "range"])
def test_no_constraint_is_named_anything_but_a_discipline_output(problem) -> None:
    """The invariant the Benders adapter imposes, on both problems.

    ``MDOScenarioAdapterBenders._compute_jacobian`` differentiates every
    sub-problem constraint by looking its **name** up among the discipline
    outputs. A constraint GEMSEO renamed -- ``constraint_name=`` for a band,
    or the ``-`` prefix ``positive=True`` produces -- is not an output of
    anything, so it raises ``KeyError``.

    The failure does not appear until the master linearises the adapter, which
    a one-iteration run never reaches, so it has to be pinned structurally.
    """
    from exlink.subdivision import (
        build_subdivided_range_scenario,
        build_subdivided_scenario,
    )

    if problem == "efficiency":
        scenario = build_subdivided_scenario(GRID, samples=SAMPLES, restoring=False)
    else:
        scenario = build_subdivided_range_scenario(
            {"q_1": 2, "theta_f": 2}, samples=SAMPLES, restoring=False
        )

    sub_scenario = scenario.formulation.sub_problem_scenario_adapter.scenario
    outputs: set[str] = set()
    stack = list(sub_scenario.disciplines)
    while stack:
        discipline = stack.pop()
        outputs.update(discipline.io.output_grammar)
        stack.extend(getattr(discipline, "disciplines", ()))

    names = sub_scenario.formulation.optimization_problem.constraints.get_names()
    assert names
    unknown = [name for name in names if name not in outputs]
    assert not unknown, f"constraints the adapter cannot differentiate: {unknown}"


def test_the_range_margins_arrive_already_negated() -> None:
    """``runs_margin >= 0`` has to reach the problem as ``runs_violation <= 0``."""
    from gemseo.core.chains.chain import MDOChain

    from exlink.disciplines import ExlinkDiscipline
    from exlink.reference import PUBLISHED_DESIGN
    from exlink.scenarios import RANGE_INEQUALITY_OUTPUTS
    from exlink.subdivision import _negated_disciplines

    _, names = _negated_disciplines(RANGE_INEQUALITY_OUTPUTS)
    assert names == ("runs_violation", "gear_violation")

    # The sign flip itself, on a discipline cheap enough to check directly.
    chain = MDOChain([
        ExlinkDiscipline(samples=SAMPLES),
        *_negated_disciplines(("tdc_gap_margin",))[0],
    ])
    output = chain.execute(PUBLISHED_DESIGN.to_mapping())
    assert float(np.ravel(output["tdc_gap_violation"])[0]) == pytest.approx(
        -float(np.ravel(output["tdc_gap_margin"])[0])
    )


def test_the_master_can_linearise_the_range_adapter() -> None:
    """The run that actually exercised the ``KeyError``, at its smallest.

    Two master iterations over two boxes is enough to reach the adapter's
    Jacobian, which is where a renamed constraint fails; a single iteration is
    not, which is how this went unnoticed.
    """
    from gemseo_box_subdivision import BoxSubdivisionSettings

    from exlink.subdivision import build_subdivided_range_scenario

    scenario = build_subdivided_range_scenario(
        {"q_1": 2, "theta_f": 2},
        samples=SAMPLES,
        seed=0,
        probes=32,
        restarts=1,
        settings=BoxSubdivisionSettings(
            convexity_margin=200.0, max_iter=3, sub_problem_max_iter=2
        ),
    )
    scenario.execute()
    assert type(scenario.formulation.sub_problem_scenario_adapter).starts
