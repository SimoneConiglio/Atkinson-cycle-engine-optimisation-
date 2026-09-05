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


def test_every_angle_axis_spans_the_720_degrees_of_a_cycle(refined_analysis) -> None:
    """A cycle is 720 deg of crankshaft, and the figures have to say so.

    The analysis is parametrised on ``theta_1``, one turn of which is one cycle,
    but power is taken from the shaft that turns twice as fast.  A turning-moment
    diagram or a piston-motion curve drawn against ``theta_1`` runs to 360 and
    cannot be laid beside a conventional engine's without a conversion the
    reader has to make.  So the abscissa is the crankshaft's, and this pins it.
    """
    from exlink.plots import crankshaft_degrees

    figure = plot_overview(REFINED_DESIGN, refined_analysis)
    motion, cycle, torque = figure.axes
    for ax in (motion, torque):
        assert ax.get_xlim() == pytest.approx((0.0, 720.0))
        assert "crankshaft angle" in ax.get_xlabel()
    # The p-V loop has no angle axis and must not have acquired one.
    assert "volume" in cycle.get_xlabel()

    theta_1 = refined_analysis.require_solved().kinematics.theta_1
    assert crankshaft_degrees(theta_1).max() == pytest.approx(
        720.0 * (theta_1.size - 1) / theta_1.size
    )


def test_the_crankshaft_torque_is_half_the_torque_referred_to_theta_1(
    refined_analysis,
) -> None:
    """``M_r`` is referred to ``theta_1``; the crankshaft sees half of it.

    Plotting ``M_r`` itself against a crankshaft abscissa would mix the two
    shafts on one pair of axes, which is the same error as running the angle to
    360.  The power is what must be invariant, and it is: half the torque at
    twice the speed.
    """
    import numpy as np

    figure = plot_torque(REFINED_DESIGN, refined_analysis)
    curves = [np.asarray(line.get_ydata(), dtype=float) for line in figure.axes[0].get_lines()]
    (drawn,) = [values for values in curves if values.size > 2]
    loads = refined_analysis.require_solved().loads
    assert drawn.max() == pytest.approx(loads.torque.max() / 2.0 / 1000.0, rel=1e-9)


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
