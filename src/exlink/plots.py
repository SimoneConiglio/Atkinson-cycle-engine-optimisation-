"""Static matplotlib views of a design and of an optimization run.

Every figure here doubles as a diagnostic:

* :func:`plot_motion` shows the four monotone phases of ``lambda(theta_1)``.
  A design that is *not* Atkinson shows up immediately -- either as a single
  up-and-down (a plain Otto motion) or as two top dead centres at visibly
  different heights, which is the ``g`` constraint made visible.
* :func:`plot_cycle` shows the approximate Atkinson cycle in the ``p-V`` plane,
  where the expansion loop being longer than the compression loop is the whole
  point of the linkage.
* :func:`plot_pareto` shows the efficiency-versus-size trade-off produced by a
  multi-objective run.
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


def _load_axes(ax: Axes, result: Any, joint: str = "R1") -> Axes:
    """Draw a joint's reaction over the revolution onto existing axes."""
    degrees = np.degrees(result.loads.kinematics.theta_1)
    magnitude = np.linalg.norm(result.loads.reaction[joint], axis=1)
    ax.plot(degrees, magnitude, lw=1.8, color="#c0392b")
    ax.set_xlabel(r"crank angle $\theta_1$ [deg]")
    ax.set_ylabel(f"$|F_{{{joint}}}|$ [N]")
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    return ax


def plot_bearing_loads(results: Sequence[Any], labels: Sequence[str]) -> Figure:
    """Compare a joint reaction across several engine speeds.

    Shows what adding inertia actually does: the mean torque is untouched, but
    the peak reaction climbs as the square of speed, and that peak is what sizes
    the parts.

    Args:
        results: Coupled results, one per speed.
        labels: A label for each.

    Returns:
        The figure.
    """
    with plt.style.context(STYLE):
        figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
        colours = plt.get_cmap("viridis")(np.linspace(0.15, 0.85, len(results)))
        for result, label, colour in zip(results, labels, colours, strict=True):
            degrees = np.degrees(result.loads.kinematics.theta_1)
            axes[0].plot(
                degrees,
                np.linalg.norm(result.loads.reaction["R1"], axis=1) / 1000.0,
                lw=1.8,
                color=colour,
                label=label,
            )
            axes[1].plot(
                degrees, result.loads.torque / 1000.0, lw=1.8, color=colour, label=label
            )
        axes[0].set_ylabel("crankshaft bearing load [kN]")
        axes[0].set_title("peak load grows as the square of speed")
        axes[1].set_ylabel(r"torque $M_r$ [N.m]")
        axes[1].set_title("mean torque is unchanged by it")
        for ax in axes:
            ax.set_xlabel(r"crank angle $\theta_1$ [deg]")
            ax.set_xlim(0.0, 360.0)
            ax.set_xticks(np.arange(0.0, 361.0, 90.0))
            ax.legend(loc="best", framealpha=0.9)
        figure.tight_layout()
    return figure


def plot_sizing(result: Any) -> Figure:
    """Show the sized sections, what drove each, and where the mass went.

    Args:
        result: A :class:`~exlink.coupled.CoupledResult`.

    Returns:
        The figure.
    """
    names = list(result.sizing)
    items = [result.sizing[n] for n in names]
    positions = np.arange(len(names))
    mode_colour = {"static": "#2980b9", "fatigue": "#c0392b", "buckling": "#e67e22"}

    with plt.style.context(STYLE):
        figure, axes = plt.subplots(1, 3, figsize=(13.0, 4.4))
        colours = [mode_colour[i.critical_mode] for i in items]

        axes[0].barh(positions, [i.diameter for i in items], color=colours)
        axes[0].set_xlabel("required diameter [mm]")
        axes[0].set_title("sections")

        axes[1].barh(positions, [1000.0 * i.mass_kg for i in items], color=colours)
        axes[1].set_xlabel("mass [g]")
        axes[1].set_title(f"total {result.total_mass_kg:.3f} kg, piston included")

        width = 0.26
        for offset, (attribute, label) in enumerate(
            (
                ("static_utilisation", "static"),
                ("fatigue_utilisation", "fatigue"),
                ("buckling_utilisation", "buckling"),
            )
        ):
            axes[2].barh(
                positions + (offset - 1) * width,
                [getattr(i, attribute) for i in items],
                height=width,
                color=mode_colour[label],
                label=label,
            )
        axes[2].axvline(1.0, color="#2c3e50", ls="--", lw=1.0)
        axes[2].set_xlabel("utilisation")
        axes[2].set_title("which mode binds")
        axes[2].legend(loc="lower right", framealpha=0.9)

        for ax in axes:
            ax.set_yticks(positions)
            ax.set_yticklabels(names)
            ax.invert_yaxis()
        figure.suptitle(
            f"sizing at {result.speed * 60.0 / (2.0 * np.pi):.0f} rpm   "
            f"({result.iterations} MDA sweeps)",
            fontsize=11,
        )
        figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    return figure


def plot_mass_vs_speed(speeds: Sequence[float], results: Sequence[Any]) -> Figure:
    """Structural mass against engine speed, on log axes.

    The slope is the story: a bending-critical member needs ``d ~ F^(1/3)``, so
    its mass goes as ``F^(2/3)`` while the inertia force it creates goes as its
    own mass -- composing to ``m ~ (C Omega^2)^3``. A cubic-in-acceleration,
    sixth-power-in-speed sensitivity is why the answer collapses so sharply.

    Args:
        speeds: Crankshaft speeds [rev/min].
        results: The coupled result at each, feasible or not.

    Returns:
        The figure.
    """
    with plt.style.context(STYLE):
        figure, ax = plt.subplots(figsize=(7.0, 4.4))
        usable = [(s, r) for s, r in zip(speeds, results, strict=True) if r.feasible]
        failed = [(s, r) for s, r in zip(speeds, results, strict=True) if not r.feasible]
        if usable:
            ax.plot(
                [s for s, _ in usable],
                [r.total_mass_kg for _, r in usable],
                "o-",
                color="#2980b9",
                lw=2.0,
                label="buildable",
            )
        if failed:
            ax.plot(
                [s for s, _ in failed],
                [r.total_mass_kg for _, r in failed],
                "x",
                ms=10,
                color="#c0392b",
                mew=2.0,
                label="loop runs away",
            )
        ax.set_xlabel("crankshaft speed [rpm]")
        ax.set_ylabel("total moving mass [kg]")
        ax.set_yscale("log")
        ax.set_title("structural mass against engine speed")
        ax.legend(loc="best", framealpha=0.9)
        figure.tight_layout()
    return figure


def plot_efficiency_mass(
    designs: Sequence[Design],
    results: Sequence[Any],
    labels: Sequence[str] | None = None,
    reference: tuple[Design, Any] | None = None,
    **kwargs: Any,
) -> Figure:
    """Efficiency against structural mass -- the trade the geometric problem
    cannot see.

    Nothing in the quasi-static formulation determines a cross-section, so mass
    is not an objective it can express. Once the parts are sized it becomes the
    one the dynamics most affects, and it pulls against efficiency through the
    transmission angle: the geometry with the longest lever arm sits nearest the
    singularity, where the accelerations -- and so the sections -- are worst.

    Args:
        designs: The designs on the trade-off curve.
        results: Their coupled results, in the same order.
        labels: Optional annotation for each point.
        reference: ``(design, coupled result)`` to mark as the starting point.
        **kwargs: Forwarded to :func:`exlink.model.analyse`.

    Returns:
        The figure.
    """
    analyses = [analyse(d, **kwargs) for d in designs]
    efficiency = np.array([100 * a.metrics.efficiency for a in analyses])
    mass = np.array([r.total_mass_kg for r in results])
    compatibility = np.array([a.metrics.compatibility for a in analyses])

    with plt.style.context(STYLE):
        figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))
        scatter = axes[0].scatter(
            mass,
            efficiency,
            c=compatibility,
            cmap="plasma",
            s=70,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
        order = np.argsort(mass)
        axes[0].plot(mass[order], efficiency[order], "-", lw=1.2, color="#7f8c8d", zorder=2)
        figure.colorbar(scatter, ax=axes[0], label=r"$W$ (1 = singular)")
        axes[0].set_xlabel("total moving mass [kg]")
        axes[0].set_ylabel(r"efficiency $\eta$ [%]")
        axes[0].set_title("the trade sizing creates")

        axes[1].scatter(
            compatibility,
            mass,
            c=efficiency,
            cmap="viridis",
            s=70,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
        axes[1].set_xlabel(r"transmission-angle margin $W$")
        axes[1].set_ylabel("total moving mass [kg]")
        axes[1].set_yscale("log")
        axes[1].set_title("and where it comes from")

        if reference is not None:
            ref_design, ref_result = reference
            ref_analysis = analyse(ref_design, **kwargs)
            axes[0].plot(
                ref_result.total_mass_kg,
                100 * ref_analysis.metrics.efficiency,
                "*",
                ms=18,
                color="#c0392b",
                mec="white",
                zorder=4,
                label="quasi-static optimum",
            )
            axes[1].plot(
                ref_analysis.metrics.compatibility,
                ref_result.total_mass_kg,
                "*",
                ms=18,
                color="#c0392b",
                mec="white",
                zorder=4,
                label="quasi-static optimum",
            )
            axes[0].legend(loc="best", framealpha=0.9)

        if labels is not None:
            for x, y, text in zip(mass, efficiency, labels, strict=True):
                axes[0].annotate(
                    text, (x, y), textcoords="offset points", xytext=(7, 5), fontsize=8
                )
        figure.tight_layout()
    return figure
