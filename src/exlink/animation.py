"""Matplotlib animation of the mechanism turning through a cycle.

An animation shows in a couple of seconds what static sketches take a chapter to
explain -- that the piston
reaches top dead centre *twice* per crankshaft revolution, and that the two
bottom dead centres differ, which is the whole point of the Atkinson linkage.

Three entry points:

:func:`animate`
    The mechanism alone.

:func:`animate_dashboard`
    The mechanism beside the piston motion, the p-V diagram and the torque
    curve, with a marker tracking the current crank angle on each -- the view
    that makes the link between geometry and thermodynamics obvious.

:func:`animate_formulations`
    Every formulation's final design, turning side by side on a common scale.
    Each panel is what a different objective converged to from the same
    starting point, so the panels differ only in what was asked of them --
    which makes the proportions readable as a consequence of the objective
    rather than as four unrelated mechanisms.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes
from matplotlib.patches import Circle

from .constants import DEFAULT_SPEC, EngineSpec
from .cycle import Phase
from .design import Design
from .model import Analysis, SolvedAnalysis, analyse
from .plots import (
    PHASE_COLOURS,
    PHASE_LABELS,
    STYLE,
    _cycle_axes,
    _motion_axes,
    _solved,
    _torque_axes,
)

LINK_STYLE: dict[str, dict[str, Any]] = {
    "crank_1": {"color": "#c0392b", "lw": 3.0, "zorder": 4},
    "crank_2": {"color": "#16a085", "lw": 3.0, "zorder": 4},
    "swing_rod": {"color": "#8e44ad", "lw": 2.6, "zorder": 5},
    "trigonal": {"color": "#2c3e50", "lw": 2.6, "zorder": 5},
    "piston_rod": {"color": "#d35400", "lw": 2.6, "zorder": 5},
}
"""Per-link line styles, shared by the static and animated views."""


def _draw_static(ax: Axes, analysis: SolvedAnalysis, spec: EngineSpec) -> None:
    """Draw everything that does not move: gears, shafts, cylinder."""
    design = analysis.design
    kinematics = analysis.kinematics
    centre_2 = kinematics.R2[0]
    ax.add_patch(
        Circle((0.0, 0.0), design.r_1, fill=False, ec="#c0392b", ls="--", lw=1.0, alpha=0.7)
    )
    ax.add_patch(
        Circle(
            tuple(centre_2),
            design.r_2,
            fill=False,
            ec="#16a085",
            ls="--",
            lw=1.0,
            alpha=0.7,
        )
    )
    ax.plot([0.0, centre_2[0]], [0.0, centre_2[1]], color="#7f8c8d", ls="-.", lw=0.9)
    ax.plot([0.0], [0.0], marker="+", ms=12, color="#c0392b", mew=1.6)
    ax.plot([centre_2[0]], [centre_2[1]], marker="+", ms=12, color="#16a085", mew=1.6)

    # Cylinder liner: the piston sweeps between these walls.
    half_bore = 0.5 * spec.bore
    lam = kinematics.lam
    bottom = float(lam.min()) - spec.piston_length
    top = float(lam.max()) + 0.25 * spec.piston_length
    for side in (-half_bore, half_bore):
        ax.plot(
            [design.x_1 + side, design.x_1 + side],
            [bottom, top],
            color="#34495e",
            lw=2.0,
            alpha=0.85,
        )
    ax.plot(
        [design.x_1 - half_bore, design.x_1 + half_bore],
        [top, top],
        color="#34495e",
        lw=2.0,
        alpha=0.85,
    )


def _mechanism_limits(analysis: SolvedAnalysis, spec: EngineSpec, margin: float = 0.08):
    """Axis limits that contain the whole mechanism over the whole revolution."""
    kinematics = analysis.kinematics
    design = analysis.design
    cloud = np.concatenate(
        [kinematics.Q, kinematics.A, kinematics.D, kinematics.E, kinematics.P, kinematics.H]
    )
    x_min = min(cloud[:, 0].min(), -design.r_1, design.x_1 - spec.bore)
    x_max = max(cloud[:, 0].max(), design.r_1, design.x_1 + spec.bore)
    y_min = min(cloud[:, 1].min(), -design.r_1)
    y_max = cloud[:, 1].max()
    span = max(x_max - x_min, y_max - y_min)
    pad = margin * span
    cx, cy = 0.5 * (x_min + x_max), 0.5 * (y_min + y_max)
    half = 0.5 * span + pad
    return (cx - half, cx + half), (cy - half, cy + half)


def _make_frame_updater(ax: Axes, analysis: SolvedAnalysis, spec: EngineSpec):
    """Create the artists and return a function that poses them at frame ``i``."""
    kinematics = analysis.kinematics
    thermo = analysis.thermodynamics
    bodies = kinematics.bodies

    lines = {
        name: ax.plot([], [], solid_capstyle="round", **LINK_STYLE[name])[0]
        for name in LINK_STYLE
    }
    (joints,) = ax.plot([], [], "o", ms=5.5, mfc="white", mec="#2c3e50", mew=1.4, zorder=6)
    (piston,) = ax.plot([], [], color="#2c3e50", lw=8.0, solid_capstyle="butt", zorder=6)
    label = ax.text(
        0.02,
        0.975,
        "",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        family="monospace",
        bbox={"fc": "white", "ec": "#bdc3c7", "alpha": 0.9, "boxstyle": "round,pad=0.4"},
    )

    half_bore = 0.5 * spec.bore

    def update(index: int):
        for name, line in lines.items():
            polyline = bodies[name][index]
            line.set_data(polyline[:, 0], polyline[:, 1])
        points = np.stack(
            [
                kinematics.Q[index],
                kinematics.A[index],
                kinematics.D[index],
                kinematics.E[index],
                kinematics.P[index],
            ]
        )
        joints.set_data(points[:, 0], points[:, 1])
        crown = kinematics.H[index, 1]
        piston.set_data(
            [analysis.design.x_1 - half_bore, analysis.design.x_1 + half_bore],
            [crown, crown],
        )
        phase = Phase(int(thermo.phases.labels[index]))
        label.set_text(
            f"theta1 = {np.degrees(kinematics.theta_1[index]):6.1f} deg\n"
            f"{PHASE_LABELS[phase]:<11}\n"
            f"p      = {thermo.gauge_pressure[index] / 0.1:6.2f} bar\n"
            f"Mr     = {analysis.loads.torque[index] / 1000.0:6.2f} N.m"
        )
        label.set_color(PHASE_COLOURS[phase])
        return (*lines.values(), joints, piston, label)

    return update


def animate(
    design: Design,
    frames: int = 180,
    interval: int = 40,
    spec: EngineSpec = DEFAULT_SPEC,
    figsize: tuple[float, float] = (6.0, 7.0),
    analysis: Analysis | None = None,
) -> FuncAnimation:
    """Animate the mechanism through one crankshaft revolution.

    Args:
        design: The mechanism to animate.
        frames: Number of crank angles in the animation.
        interval: Delay between frames [ms].
        spec: Fixed engine data.
        figsize: Figure size [in].
        analysis: A precomputed analysis; recomputed at ``frames`` resolution if
            omitted.

    Returns:
        The animation.  Keep a reference to it, or matplotlib will garbage
        collect it before it plays; :func:`save` handles that for you.

    Raises:
        ValueError: If the design cannot be analysed.
    """
    solved = _solved(analysis or analyse(design, samples=frames, spec=spec))

    with plt.style.context(STYLE):
        figure, ax = plt.subplots(figsize=figsize)
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title("EX-link mechanism")
        xlim, ylim = _mechanism_limits(solved, spec)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        _draw_static(ax, solved, spec)
        update = _make_frame_updater(ax, solved, spec)
        figure.tight_layout()

    n_frames = solved.kinematics.theta_1.size
    return FuncAnimation(
        figure, update, frames=n_frames, interval=interval, blit=True, repeat=True
    )


def animate_dashboard(
    design: Design,
    frames: int = 180,
    interval: int = 40,
    spec: EngineSpec = DEFAULT_SPEC,
    figsize: tuple[float, float] = (12.0, 7.5),
    analysis: Analysis | None = None,
) -> FuncAnimation:
    """Animate the mechanism alongside its motion, cycle and torque curves.

    A marker on each curve tracks the crank angle shown by the linkage, so the
    two top dead centres in the piston motion, the two loops of the p-V diagram
    and the torque peaks can be read off against the geometry that produces
    them.

    Args:
        design: The mechanism to animate.
        frames: Number of crank angles in the animation.
        interval: Delay between frames [ms].
        spec: Fixed engine data.
        figsize: Figure size [in].
        analysis: A precomputed analysis.

    Returns:
        The animation.

    Raises:
        ValueError: If the design cannot be analysed.
    """
    solved = _solved(analysis or analyse(design, samples=frames, spec=spec))
    kinematics = solved.kinematics
    thermo = solved.thermodynamics
    loads = solved.loads

    with plt.style.context(STYLE):
        figure = plt.figure(figsize=figsize)
        # Explicit margins rather than tight_layout: the mechanism panel spans all
        # three rows, which tight_layout cannot reconcile with the suptitle.
        grid = figure.add_gridspec(
            3,
            2,
            width_ratios=[1.0, 1.25],
            hspace=0.55,
            wspace=0.22,
            left=0.06,
            right=0.98,
            bottom=0.07,
            top=0.92,
        )
        ax_mech = figure.add_subplot(grid[:, 0])
        ax_motion = figure.add_subplot(grid[0, 1])
        ax_cycle = figure.add_subplot(grid[1, 1])
        ax_torque = figure.add_subplot(grid[2, 1])

        ax_mech.set_aspect("equal")
        ax_mech.set_xlabel("x [mm]")
        ax_mech.set_ylabel("y [mm]")
        ax_mech.set_title("EX-link mechanism")
        xlim, ylim = _mechanism_limits(solved, spec)
        ax_mech.set_xlim(*xlim)
        ax_mech.set_ylim(*ylim)
        _draw_static(ax_mech, solved, spec)
        update_mechanism = _make_frame_updater(ax_mech, solved, spec)

        _motion_axes(ax_motion, solved)
        _cycle_axes(ax_cycle, solved)
        _torque_axes(ax_torque, solved)

        markers = [
            ax_motion.plot([], [], "o", ms=7, color="#2c3e50", zorder=9)[0],
            ax_cycle.plot([], [], "o", ms=7, color="#2c3e50", zorder=9)[0],
            ax_torque.plot([], [], "o", ms=7, color="#2c3e50", zorder=9)[0],
        ]
        figure.suptitle(
            f"eta = {100 * solved.metrics.efficiency:.2f} %   "
            f"H = {solved.metrics.height:.1f} mm   "
            f"B = {solved.metrics.width:.1f} mm   "
            f"STE = {solved.metrics.expansion_stroke:.2f} mm   "
            f"eps = {solved.metrics.compression_ratio:.2f}",
            fontsize=10,
        )

    degrees = np.degrees(kinematics.theta_1)

    def update(index: int):
        artists = update_mechanism(index)
        markers[0].set_data([degrees[index]], [kinematics.lam[index]])
        markers[1].set_data([thermo.volume[index] / 1000.0], [thermo.pressure[index] / 0.1])
        markers[2].set_data([degrees[index]], [loads.torque[index] / 1000.0])
        return (*artists, *markers)

    n_frames = kinematics.theta_1.size
    return FuncAnimation(
        figure, update, frames=n_frames, interval=interval, blit=True, repeat=True
    )


def animate_formulations(
    designs: dict[str, Design] | None = None,
    frames: int = 120,
    interval: int = 40,
    spec: EngineSpec = DEFAULT_SPEC,
    panel_size: tuple[float, float] = (3.2, 4.3),
    common_scale: bool = True,
) -> FuncAnimation:
    """Animate each formulation's final design, side by side.

    The three static views of a converged design -- a table of metrics, a
    mechanism sketch, a mass budget -- all answer "what did it reach?".  None
    of them answers "what does the objective *do* to the mechanism?", which is
    the question this study is actually about.  Four linkages turning together
    answer it directly: the geometric optimum is visibly long-limbed and close
    to alignment, the coupled optimum visibly shorter and squatter, and the
    difference is the inertia loading that only the second formulation could
    see.

    All panels share one set of axis limits by default, because the point is
    that the designs are *different sizes*; scaling each to fit its own frame
    would hide exactly the effect being shown.

    Args:
        designs: Label to design, in the order the panels should appear.
            Defaults to the four results of :mod:`exlink.reference`, each
            labelled with its formulation and its headline number.
        frames: Crank angles in the animation.
        interval: Delay between frames [ms].
        spec: Fixed engine data.
        panel_size: Size of one panel [in].
        common_scale: Share one set of axis limits across the panels.  Set
            ``False`` only to inspect a single mechanism's detail.

    Returns:
        The animation.  Keep a reference to it, or matplotlib will garbage
        collect it; :func:`save` handles that.

    Raises:
        ValueError: If a design cannot be analysed, or none is given.
    """
    from .reference import (
        COUPLED_DESIGN,
        GRADIENT_DESIGN,
        REFINED_DESIGN,
        RELIABLE_DESIGN,
    )

    chosen = (
        designs
        if designs is not None
        else {
            "geometric, augmented Lagrangian\n$\\eta$ = 27.9 %": REFINED_DESIGN,
            "geometric, SLSQP + exact gradients\n$\\eta$ = 30.8 %": GRADIENT_DESIGN,
            "coupled, minimum mass\n0.234 kg of moving mass": COUPLED_DESIGN,
            "range + reliability\n3395 km/L at $P_f = 1.3\\times10^{-3}$": RELIABLE_DESIGN,
        }
    )
    if not chosen:
        msg = "no designs to animate"
        raise ValueError(msg)

    solved = {
        label: _solved(analyse(design, samples=frames, spec=spec))
        for label, design in chosen.items()
    }

    limits = [_mechanism_limits(analysis, spec) for analysis in solved.values()]
    if common_scale:
        # The union of the individual boxes: the smallest window that holds
        # every mechanism.  It is deliberately *not* padded out to a square --
        # each panel then wastes the same space, and a mechanism that fills
        # less of the frame is one that is genuinely smaller.
        shared = (
            (min(x[0] for x, _y in limits), max(x[1] for x, _y in limits)),
            (min(y[0] for _x, y in limits), max(y[1] for _x, y in limits)),
        )
        limits = [shared] * len(limits)

    count = len(solved)
    with plt.style.context(STYLE):
        figure, axes = plt.subplots(
            1,
            count,
            figsize=(panel_size[0] * count, panel_size[1]),
            squeeze=False,
        )
        updaters = []
        for ax, (label, analysis), (xlim, ylim) in zip(
            axes[0], solved.items(), limits, strict=True
        ):
            ax.set_aspect("equal")
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_title(label, fontsize=9)
            ax.set_xlabel("x [mm]")
            ax.tick_params(labelsize=8)
            _draw_static(ax, analysis, spec)
            updaters.append(_make_frame_updater(ax, analysis, spec))
        axes[0][0].set_ylabel("y [mm]")
        for ax in axes[0][1:]:
            ax.set_yticklabels([])
        figure.suptitle(
            "the same mechanism under four objectives, at the same crank angle",
            fontsize=11,
        )
        figure.tight_layout(rect=(0.0, 0.02, 1.0, 0.94))

    def update(index: int):
        artists: list[Any] = []
        for updater in updaters:
            artists.extend(updater(index))
        return tuple(artists)

    n_frames = min(a.kinematics.theta_1.size for a in solved.values())
    return FuncAnimation(
        figure, update, frames=n_frames, interval=interval, blit=True, repeat=True
    )


def save(
    animation: FuncAnimation,
    path: str | Path,
    fps: int = 25,
    dpi: int = 110,
    **kwargs: Any,
) -> Path:
    """Write an animation to disk, picking a writer from the file extension.

    ``.gif`` uses Pillow and always works; ``.mp4`` needs ffmpeg on the PATH.

    Args:
        animation: The animation to save.
        path: Destination; ``.gif`` or ``.mp4``.
        fps: Frames per second.
        dpi: Resolution.
        **kwargs: Forwarded to :meth:`~matplotlib.animation.Animation.save`.

    Returns:
        The path written.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = "pillow" if path.suffix.lower() == ".gif" else "ffmpeg"
    animation.save(str(path), writer=writer, fps=fps, dpi=dpi, **kwargs)
    plt.close(animation._fig)  # type: ignore[attr-defined]  # no public accessor
    return path
