"""Figures and animation. Smoke tests: they must build without raising."""

from __future__ import annotations

import matplotlib
import pytest
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

matplotlib.use("Agg")

from exlink import PUBLISHED_DESIGN
from exlink.animation import animate, animate_dashboard, animate_formulations, save
from exlink.plots import (
    plot_cycle,
    plot_mechanism,
    plot_motion,
    plot_overview,
    plot_pareto,
    plot_torque,
)
from exlink.reference import REFINED_DESIGN


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.mark.parametrize(
    "builder",
    [plot_motion, plot_cycle, plot_torque, plot_overview],
    ids=["motion", "cycle", "torque", "overview"],
)
def test_single_design_figures_build(builder, refined_analysis) -> None:
    figure = builder(REFINED_DESIGN, refined_analysis)
    assert isinstance(figure, Figure)
    assert figure.axes


def test_the_mechanism_figure_has_one_panel_per_angle(refined_analysis) -> None:
    figure = plot_mechanism(
        REFINED_DESIGN, crank_angles=(0.0, 120.0, 240.0), analysis=refined_analysis
    )
    assert len(figure.axes) == 3


def test_the_pareto_figure_builds(refined_analysis) -> None:
    designs = [
        REFINED_DESIGN,
        REFINED_DESIGN.replace(a=REFINED_DESIGN.a * 1.02),
        REFINED_DESIGN.replace(e=REFINED_DESIGN.e * 1.02),
    ]
    figure = plot_pareto(designs, highlight=REFINED_DESIGN, samples=360)
    assert isinstance(figure, Figure)


def test_plotting_an_unanalysable_design_raises() -> None:
    broken = PUBLISHED_DESIGN.replace(a=25.0, c=25.0)
    with pytest.raises(ValueError, match="unanalysable"):
        plot_motion(broken, samples=180)


def test_the_pareto_figure_needs_a_valid_design() -> None:
    with pytest.raises(ValueError, match="no analysable design"):
        plot_pareto([PUBLISHED_DESIGN.replace(a=25.0, c=25.0)], samples=180)


def test_the_animation_has_one_frame_per_crank_angle() -> None:
    animation = animate(REFINED_DESIGN, frames=24)
    assert len(list(animation.new_frame_seq())) == 24
    plt.close(animation._fig)


def test_the_dashboard_animation_builds() -> None:
    animation = animate_dashboard(REFINED_DESIGN, frames=24)
    assert len(animation._fig.axes) >= 4
    plt.close(animation._fig)


def test_the_formulation_panels_share_one_scale() -> None:
    """The comparison is only readable if the panels are not each rescaled.

    A design that fills less of its frame is a smaller mechanism; scaling each
    panel to its own contents would erase exactly the difference the figure
    exists to show.
    """
    animation = animate_formulations(frames=12)
    axes = animation._fig.axes
    assert len(axes) == 4
    limits = {(ax.get_xlim(), ax.get_ylim()) for ax in axes}
    assert len(limits) == 1
    plt.close(animation._fig)


def test_the_formulation_panels_show_different_designs() -> None:
    """Four panels of the same mechanism would be a bug, not a comparison."""
    from exlink.reference import COUPLED_DESIGN, REFINED_DESIGN

    animation = animate_formulations(
        designs={"a": REFINED_DESIGN, "b": COUPLED_DESIGN}, frames=12
    )
    axes = animation._fig.axes
    assert len(axes) == 2
    assert [ax.get_title() for ax in axes] == ["a", "b"]
    plt.close(animation._fig)


def test_animating_no_formulations_raises() -> None:
    with pytest.raises(ValueError, match="no designs"):
        animate_formulations(designs={}, frames=12)


def test_animating_an_unanalysable_design_raises() -> None:
    with pytest.raises(ValueError, match="unanalysable"):
        animate(PUBLISHED_DESIGN.replace(a=25.0, c=25.0), frames=12)


@pytest.mark.slow
def test_the_animation_writes_a_gif(tmp_path) -> None:
    animation = animate(REFINED_DESIGN, frames=12)
    path = save(animation, tmp_path / "out" / "exlink.gif", fps=8, dpi=60)
    assert path.is_file()
    assert path.stat().st_size > 0
