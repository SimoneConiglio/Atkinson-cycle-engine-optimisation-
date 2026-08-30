"""Static matplotlib views of a design and of an optimization run.

Every figure here is also a diagnostic used while reading the report:

* :func:`plot_motion` shows the four monotone phases of ``lambda(theta_1)``.
  A design that is *not* Atkinson shows up immediately -- either as a single
  up-and-down (a plain Otto motion) or as two top dead centres at visibly
  different heights, which is the ``g`` constraint made visible.
* :func:`plot_cycle` shows the approximate Atkinson cycle in the ``p-V`` plane,
  where the expansion loop being longer than the compression loop is the whole
  point of the linkage.
* :func:`plot_pareto` shows the efficiency-versus-size trade-off that the
  report's multi-objective run produces.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .constants import DEFAULT_SPEC, EngineSpec
from .cycle import Phase
from .design import Design
from .kinematics import solve as solve_kinematics
from .model import Analysis, SolvedAnalysis, analyse

STYLE: dict[str, Any] = {
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": ":",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.facecolor": "white",
}
"""A light, uncluttered style applied by every plotting helper."""

PHASE_COLOURS: dict[Phase, str] = {
    Phase.EXPANSION: "#c0392b",
    Phase.EXHAUST: "#7f8c8d",
    Phase.INTAKE: "#2980b9",
    Phase.COMPRESSION: "#e67e22",
}

PHASE_LABELS: dict[Phase, str] = {
    Phase.EXPANSION: "expansion",
    Phase.EXHAUST: "exhaust",
    Phase.INTAKE: "intake",
    Phase.COMPRESSION: "compression",
}


def _solved(analysis: Analysis) -> SolvedAnalysis:
    """Narrow an analysis, refusing to plot penalised placeholder values."""
    if not analysis.valid:
        msg = f"cannot plot an unanalysable design: {analysis.metrics.reason}"
        raise ValueError(msg)
    return analysis.require_solved()


def _segments(labels: np.ndarray) -> list[tuple[int, int, Phase]]:
    """Split the revolution into contiguous runs of one phase each."""
    runs: list[tuple[int, int, Phase]] = []
    start = 0
    for index in range(1, labels.size + 1):
        if index == labels.size or labels[index] != labels[start]:
            runs.append((start, index, Phase(int(labels[start]))))
            start = index
    return runs


def _motion_axes(ax: Axes, analysis: SolvedAnalysis) -> Axes:
    """Draw ``lambda(theta_1)``, coloured by phase, onto existing axes."""
    kinematics = analysis.kinematics
    thermo = analysis.thermodynamics
    degrees = np.degrees(kinematics.theta_1)
    seen: set[Phase] = set()
    for start, stop, phase in _segments(thermo.phases.labels):
        stop = min(stop + 1, degrees.size)
        ax.plot(
            degrees[start:stop],
            kinematics.lam[start:stop],
            color=PHASE_COLOURS[phase],
            lw=2.0,
            label=PHASE_LABELS[phase] if phase not in seen else None,
        )
        seen.add(phase)
    ax.axhline(
        thermo.phases.lam_tdc,
        color="#2c3e50",
        ls="--",
        lw=0.9,
        alpha=0.7,
        label="top dead centre",
    )
    ax.set_xlabel(r"crank angle $\theta_1$ [deg]")
    ax.set_ylabel(r"piston height $\lambda$ [mm]")
    ax.set_title(
        f"piston motion   STE = {analysis.metrics.expansion_stroke:.2f} mm, "
        f"STC = {analysis.metrics.compression_stroke:.2f} mm, "
        f"g = {analysis.metrics.tdc_gap:.4f} mm"
    )
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    ax.legend(loc="lower right", ncol=2, framealpha=0.9)
    return ax


def _cycle_axes(ax: Axes, analysis: SolvedAnalysis) -> Axes:
    """Draw the ``p-V`` diagram onto existing axes."""
    thermo = analysis.thermodynamics
    volume = thermo.volume / 1000.0  # cc
    pressure = thermo.pressure / 0.1  # bar
    order = np.argsort(analysis.kinematics.theta_1)
    closed_v = np.append(volume[order], volume[order][0])
    closed_p = np.append(pressure[order], pressure[order][0])
    ax.plot(closed_v, closed_p, color="#2c3e50", lw=1.2, alpha=0.6)
    # A phase that straddles theta_1 = 0 comes back as two runs; label it once.
    seen: set[Phase] = set()
    for start, stop, phase in _segments(thermo.phases.labels):
        stop = min(stop + 1, volume.size)
        ax.plot(
            volume[start:stop],
            pressure[start:stop],
            color=PHASE_COLOURS[phase],
            lw=2.2,
            label=PHASE_LABELS[phase] if phase not in seen else None,
        )
        seen.add(phase)
    ax.set_xlabel("volume $V$ [cc]")
    ax.set_ylabel("pressure $p$ [bar]")
    ax.set_title(
        f"Atkinson cycle   " r"$\epsilon$ = " f"{analysis.metrics.compression_ratio:.2f}"
    )
    return ax


def _torque_axes(ax: Axes, analysis: SolvedAnalysis) -> Axes:
    """Draw the crankshaft torque onto existing axes."""
    loads = analysis.loads
    degrees = np.degrees(analysis.kinematics.theta_1)
    ax.plot(degrees, loads.torque / 1000.0, color="#8e44ad", lw=1.8)
    ax.axhline(
        loads.mean_torque / 1000.0,
        color="#c0392b",
        ls="--",
        lw=1.0,
        label=f"mean = {loads.mean_torque / 1000.0:.2f} N.m",
    )
    ax.axhline(0.0, color="#7f8c8d", lw=0.8)
    ax.set_xlabel(r"crank angle $\theta_1$ [deg]")
    ax.set_ylabel(r"torque $M_r$ [N.m]")
    ax.set_title(
        "crankshaft torque   "
        r"$\eta$ = "
        f"{100 * analysis.metrics.efficiency:.2f} %"
    )
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    ax.legend(loc="best", framealpha=0.9)
    return ax


def plot_motion(design: Design, analysis: Analysis | None = None, **kwargs: Any) -> Figure:
    """Plot the piston motion ``lambda(theta_1)`` over one revolution."""
    solved = _solved(analysis or analyse(design, **kwargs))
    with plt.style.context(STYLE):
        figure, ax = plt.subplots(figsize=(7.0, 4.0))
        _motion_axes(ax, solved)
        figure.tight_layout()
    return figure


def plot_cycle(design: Design, analysis: Analysis | None = None, **kwargs: Any) -> Figure:
    """Plot the approximate Atkinson cycle in the ``p-V`` plane."""
    solved = _solved(analysis or analyse(design, **kwargs))
    with plt.style.context(STYLE):
        figure, ax = plt.subplots(figsize=(6.0, 4.5))
        _cycle_axes(ax, solved)
        ax.legend(loc="best", framealpha=0.9)
        figure.tight_layout()
    return figure


def plot_torque(design: Design, analysis: Analysis | None = None, **kwargs: Any) -> Figure:
    """Plot the crankshaft torque over one revolution."""
    solved = _solved(analysis or analyse(design, **kwargs))
    with plt.style.context(STYLE):
        figure, ax = plt.subplots(figsize=(7.0, 4.0))
        _torque_axes(ax, solved)
        figure.tight_layout()
    return figure


def plot_mechanism(
    design: Design,
    crank_angles: Sequence[float] = (0.0, 90.0, 180.0, 270.0),
    spec: EngineSpec = DEFAULT_SPEC,
    analysis: Analysis | None = None,
    **kwargs: Any,
) -> Figure:
    """Draw the mechanism at several crank angles side by side.

    Args:
        design: The mechanism to draw.
        crank_angles: Crank angles to show [deg].
        spec: Fixed engine data.
        analysis: A precomputed analysis at full resolution.
        **kwargs: Forwarded to :func:`exlink.model.analyse`.

    Returns:
        The figure.
    """
    from .animation import LINK_STYLE, _draw_static, _mechanism_limits

    solved = _solved(analysis or analyse(design, spec=spec, **kwargs))

    # Re-solve at exactly the requested angles rather than snapping to the grid.
    angles = np.radians(np.asarray(crank_angles, dtype=float))
    frames = solve_kinematics(design, spec=spec, theta_1=angles)
    bodies = frames.bodies
    xlim, ylim = _mechanism_limits(solved, spec)

    with plt.style.context(STYLE):
        figure, axes = plt.subplots(
            1, len(angles), figsize=(3.1 * len(angles), 4.4), squeeze=False
        )
        for column, (ax, angle) in enumerate(zip(axes[0], crank_angles, strict=True)):
            ax.set_aspect("equal")
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_title(r"$\theta_1$ = " f"{angle:.0f}" r"$^\circ$")
            _draw_static(ax, solved, spec)
            for name, style in LINK_STYLE.items():
                polyline = bodies[name][column]
                ax.plot(polyline[:, 0], polyline[:, 1], solid_capstyle="round", **style)
            crown = frames.H[column, 1]
            half = 0.5 * spec.bore
            ax.plot(
                [design.x_1 - half, design.x_1 + half],
                [crown, crown],
                color="#2c3e50",
                lw=7.0,
                solid_capstyle="butt",
            )
            if column:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("y [mm]")
            ax.set_xlabel("x [mm]")
        figure.tight_layout()
    return figure


def plot_overview(design: Design, analysis: Analysis | None = None, **kwargs: Any) -> Figure:
    """Motion, cycle and torque on one page."""
    solved = _solved(analysis or analyse(design, **kwargs))
    with plt.style.context(STYLE):
        figure, axes = plt.subplots(3, 1, figsize=(7.5, 10.0))
        _motion_axes(axes[0], solved)
        _cycle_axes(axes[1], solved)
        axes[1].legend(loc="best", framealpha=0.9)
        _torque_axes(axes[2], solved)
        figure.tight_layout()
    return figure


def plot_pareto(
    designs: Sequence[Design],
    analyses: Sequence[Analysis] | None = None,
    highlight: Design | None = None,
    **kwargs: Any,
) -> Figure:
    """Plot a Pareto front of ``(eta, H, B)``.

    Two panels: efficiency against each envelope dimension, with the marker
    colour carrying the third objective, so the whole three-dimensional front
    is readable on a flat page.

    Args:
        designs: The Pareto-optimal designs.
        analyses: Their analyses; recomputed if omitted.
        highlight: A design to mark, typically the one finally selected.
        **kwargs: Forwarded to :func:`exlink.model.analyse`.

    Returns:
        The figure.

    Raises:
        ValueError: If no valid design is given.
    """
    analyses = list(analyses or [analyse(d, **kwargs) for d in designs])
    valid = [a for a in analyses if a.valid]
    if not valid:
        msg = "no analysable design in the front"
        raise ValueError(msg)

    eta = np.array([100 * a.metrics.efficiency for a in valid])
    height = np.array([a.metrics.height for a in valid])
    width = np.array([a.metrics.width for a in valid])

    with plt.style.context(STYLE):
        figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))
        for ax, size, other, xlabel, clabel in (
            (axes[0], height, width, "envelope along stroke $H$ [mm]", "$B$ [mm]"),
            (axes[1], width, height, "envelope across stroke $B$ [mm]", "$H$ [mm]"),
        ):
            scatter = ax.scatter(
                size, eta, c=other, cmap="viridis", s=34, edgecolor="white", linewidth=0.5
            )
            figure.colorbar(scatter, ax=ax, label=clabel)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r"efficiency $\eta$ [%]")
        if highlight is not None:
            chosen = analyse(highlight, **kwargs)
            if chosen.valid:
                axes[0].plot(
                    chosen.metrics.height,
                    100 * chosen.metrics.efficiency,
                    "*",
                    ms=18,
                    color="#c0392b",
                    mec="white",
                    label="selected",
                )
                axes[1].plot(
                    chosen.metrics.width,
                    100 * chosen.metrics.efficiency,
                    "*",
                    ms=18,
                    color="#c0392b",
                    mec="white",
                    label="selected",
                )
                axes[0].legend(loc="best")
        figure.suptitle(f"Pareto front, {len(valid)} designs", fontsize=11)
        figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    return figure


def plot_convergence(outcome: Any) -> Figure:
    """Plot the objective history of a solved scenario.

    Args:
        outcome: An :class:`exlink.scenarios.Outcome`.

    Returns:
        The figure.
    """
    problem = outcome.scenario.formulation.optimization_problem
    history = np.array(
        [
            np.ravel(value)[0]
            for value in problem.database.get_function_history(
                problem.objective.name, with_x_vect=False
            )
        ]
    )
    with plt.style.context(STYLE):
        figure, ax = plt.subplots(figsize=(7.0, 4.0))
        ax.plot(np.arange(1, history.size + 1), -history, lw=1.4, color="#2980b9")
        ax.plot(
            np.arange(1, history.size + 1),
            np.maximum.accumulate(-history),
            lw=2.0,
            color="#c0392b",
            label="best so far",
        )
        ax.set_xlabel("evaluation")
        ax.set_ylabel(r"efficiency $\eta$ [-]")
        ax.set_title(f"{outcome.algorithm}: {history.size} evaluations")
        ax.legend(loc="best")
        figure.tight_layout()
    return figure
