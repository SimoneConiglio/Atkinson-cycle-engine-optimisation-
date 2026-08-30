"""The command-line interface."""

from __future__ import annotations

import json

import pytest

from exlink import PUBLISHED_DESIGN, VARIABLE_NAMES
from exlink.cli import build_parser, load_design, main, save_design
from exlink.reference import REFINED_DESIGN


def test_the_parser_requires_a_command() -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args([])


def test_named_designs_resolve() -> None:
    assert load_design("published") == PUBLISHED_DESIGN
    assert load_design("refined") == REFINED_DESIGN
    assert load_design(None) == REFINED_DESIGN


def test_an_unknown_design_is_reported_clearly() -> None:
    with pytest.raises(SystemExit, match="no such design"):
        load_design("does-not-exist.json")


def test_a_design_round_trips_through_json(tmp_path) -> None:
    path = save_design(PUBLISHED_DESIGN, tmp_path / "d.json")
    assert set(json.loads(path.read_text())) == set(VARIABLE_NAMES)
    assert load_design(str(path)) == PUBLISHED_DESIGN


def test_a_design_can_be_read_from_a_bare_list(tmp_path) -> None:
    path = tmp_path / "d.json"
    path.write_text(json.dumps(list(PUBLISHED_DESIGN.to_array())))
    assert load_design(str(path)) == PUBLISHED_DESIGN


def test_analyse_prints_a_report(capsys) -> None:
    code = main(["analyse", "--design", "refined", "--samples", "360"])
    output = capsys.readouterr().out
    assert code == 0
    assert "efficiency" in output
    assert "feasible" in output


def test_analyse_reports_failure_for_a_penalised_design(tmp_path, capsys) -> None:
    path = save_design(PUBLISHED_DESIGN.replace(a=25.0, c=25.0), tmp_path / "bad.json")
    code = main(["analyse", "--design", str(path), "--samples", "180"])
    assert code == 1
    assert "PENALISED" in capsys.readouterr().out


def test_analyse_can_save_the_design(tmp_path) -> None:
    target = tmp_path / "saved.json"
    main(["analyse", "--design", "published", "--samples", "180", "--save", str(target)])
    assert load_design(str(target)) == PUBLISHED_DESIGN


def test_plot_writes_every_figure(tmp_path) -> None:
    outdir = tmp_path / "figures"
    code = main(
        ["plot", "--design", "refined", "--samples", "360", "-o", str(outdir), "--dpi", "60"]
    )
    assert code == 0
    written = {p.name for p in outdir.glob("*.png")}
    assert written == {"motion.png", "cycle.png", "torque.png", "mechanism.png", "overview.png"}


@pytest.mark.slow
def test_animate_writes_a_file(tmp_path) -> None:
    target = tmp_path / "a.gif"
    code = main(
        [
            "animate",
            "--design",
            "refined",
            "--frames",
            "10",
            "--fps",
            "8",
            "--dpi",
            "60",
            "-o",
            str(target),
        ]
    )
    assert code == 0
    assert target.is_file()
