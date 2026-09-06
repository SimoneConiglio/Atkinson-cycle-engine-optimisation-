"""GEMSEO's N2 chart and XDSM of the coupled problem."""

from __future__ import annotations

import shutil

import pytest

from exlink.diagrams import build_published_scenario, write_n2, write_xdsm


@pytest.fixture(scope="module")
def scenario():
    return build_published_scenario()


def test_the_scenario_the_diagrams_are_drawn_from_is_the_one_solved(scenario) -> None:
    """A hand-drawn graph drifts; one generated from the scenario cannot.

    Both pictures come from the object the optimizer is handed, so this pins
    the disciplines they will show and the order section 4.2 reads off them.
    """
    names = [type(discipline).__name__ for discipline in scenario.disciplines]
    assert names == [
        "ExlinkDiscipline",
        "DynamicsDiscipline",
        "StructureDiscipline",
        "BearingMarginDiscipline",
        "RangeDiscipline",
    ]


def test_the_coupling_that_needs_an_mda_is_the_one_the_n2_will_show(scenario) -> None:
    """One block below the diagonal, and it is the sizing/inertia fixed point.

    The N2 chart is an image, so what it *means* has to be checked on the graph
    behind it: the dynamics and structure disciplines must exchange variables
    in both directions, and nothing else may, or the figure's caption -- and
    section 3.6's argument for MDF -- would be wrong.
    """
    by_name = {type(d).__name__: d for d in scenario.disciplines}
    dynamics = by_name["DynamicsDiscipline"]
    structure = by_name["StructureDiscipline"]

    forward = set(dynamics.io.output_grammar) & set(structure.io.input_grammar)
    backward = set(structure.io.output_grammar) & set(dynamics.io.input_grammar)
    assert forward, "dynamics must feed structure"
    assert backward, "structure must feed dynamics -- this is the feedback"
    assert "diameters" in backward

    # No other pair closes a loop: every remaining discipline is feed-forward.
    others = ["ExlinkDiscipline", "BearingMarginDiscipline", "RangeDiscipline"]
    for name in others:
        downstream = by_name[name]
        loops = [
            other
            for other in scenario.disciplines
            if other is not downstream
            and set(downstream.io.output_grammar) & set(other.io.input_grammar)
            and set(other.io.output_grammar) & set(downstream.io.input_grammar)
        ]
        assert not loops, f"{name} unexpectedly sits in a loop with {loops}"

    # And the geometry discipline shares nothing at all: section 4.2's caption
    # calls its row empty, which is a claim about the graph, not about the PNG.
    geometry = by_name["ExlinkDiscipline"]
    for other in scenario.disciplines:
        if other is geometry:
            continue
        assert not set(geometry.io.output_grammar) & set(other.io.input_grammar)
        assert not set(other.io.output_grammar) & set(geometry.io.input_grammar)


@pytest.mark.slow
def test_the_n2_chart_is_written(tmp_path, scenario) -> None:
    path = write_n2(tmp_path / "n2.png", scenario=scenario)
    assert path.is_file()
    assert path.stat().st_size > 5_000


@pytest.mark.slow
def test_the_xdsm_is_written_in_the_best_format_the_tools_allow(tmp_path, scenario) -> None:
    """pyXDSM where LaTeX exists, HTML where it does not -- but always a file.

    The published figure is the rasterised pyXDSM one; the HTML fallback keeps
    the command usable on a machine without a TeX toolchain, which is most of
    them, so both paths have to work.
    """
    path = write_xdsm(tmp_path, scenario=scenario)
    assert path.is_file()
    if shutil.which("pdflatex") is None:
        assert path.suffix == ".html"
    elif shutil.which("pdftoppm") is None:
        assert path.suffix == ".pdf"
    else:
        assert path.suffix == ".png"
        # Trimmed to its ink, so the figure is not mostly margin.
        from PIL import Image

        with Image.open(path) as image:
            assert image.width > image.height
            assert image.width > 800
