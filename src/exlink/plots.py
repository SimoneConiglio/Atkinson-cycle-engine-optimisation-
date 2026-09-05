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

import math
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
from .materials import FloatArray
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

JOINT_NAMES: tuple[str, ...] = ("R1", "R2", "Q", "A", "D", "E", "P", "H")
"""Joints of the linkage, in the order :mod:`exlink.kinematics` defines them."""

JOINT_LABEL_OFFSET: dict[str, tuple[float, float]] = {
    "R1": (8.0, -12.0),
    "R2": (8.0, -4.0),
    "Q": (-16.0, 2.0),
    "A": (8.0, 4.0),
    "D": (8.0, -8.0),
    "E": (-16.0, -6.0),
    "P": (8.0, 0.0),
    "H": (8.0, 4.0),
}
"""Where each joint label sits relative to its joint, in points."""

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


def _dimension(
    ax: Axes,
    start: Sequence[float] | FloatArray,
    end: Sequence[float] | FloatArray,
    label: str,
    offset: float = 0.0,
    colour: str = "#2c3e50",
    fontsize: float = 11.0,
) -> None:
    """Draw a double-headed dimension between two points and label its middle.

    Args:
        ax: Axes to draw on.
        start: First point ``(x, y)`` [mm].
        end: Second point [mm].
        label: Text to place at the midpoint, normally a LaTeX symbol.
        offset: Perpendicular shift of the dimension line [mm], so the arrow
            can be pulled clear of the member it measures.
        colour: Line and text colour.
        fontsize: Label size in points.
    """
    first = np.asarray(start, dtype=float)
    second = np.asarray(end, dtype=float)
    span = second - first
    length = float(np.hypot(*span))
    if length <= 0.0:
        return
    normal = np.array([-span[1], span[0]]) / length
    shift = normal * offset
    ax.annotate(
        "",
        xy=tuple(second + shift),
        xytext=tuple(first + shift),
        arrowprops={
            "arrowstyle": "<|-|>",
            "color": colour,
            "lw": 1.1,
            "shrinkA": 0,
            "shrinkB": 0,
        },
        zorder=8,
    )
    middle = 0.5 * (first + second) + shift + normal * (3.5 if offset >= 0.0 else -3.5)
    ax.text(
        float(middle[0]),
        float(middle[1]),
        label,
        color=colour,
        fontsize=fontsize,
        ha="center",
        va="center",
        zorder=9,
        bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
    )


def _angle_arc(
    ax: Axes,
    centre: Sequence[float] | FloatArray,
    angle: float,
    label: str,
    radius: float = 26.0,
    colour: str = "#7f8c8d",
) -> None:
    """Draw an arc from the +x direction to ``angle`` and label it."""
    origin = np.asarray(centre, dtype=float)
    sweep = np.linspace(0.0, angle, 64)
    ax.plot(
        origin[0] + radius * np.cos(sweep),
        origin[1] + radius * np.sin(sweep),
        color=colour,
        lw=1.0,
        ls="-",
        zorder=7,
    )
    ax.plot(
        [origin[0], origin[0] + 1.35 * radius],
        [origin[1], origin[1]],
        color=colour,
        lw=0.8,
        ls=":",
        zorder=6,
    )
    tip = origin + 1.16 * radius * np.array([np.cos(0.5 * angle), np.sin(0.5 * angle)])
    ax.text(
        float(tip[0]),
        float(tip[1]),
        label,
        color=colour,
        fontsize=11.0,
        ha="center",
        va="center",
        zorder=9,
        bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.85},
    )


def _trigonal_inset(ax: Axes, design: Design, colour: str = "#34495e") -> None:
    """Inset showing ``E`` in the frame the trigonal link carries.

    ``x_b`` and ``y_b`` are the coordinates of ``E`` in the frame whose origin
    is ``A`` and whose first axis runs along ``AD``.  On the mechanism itself
    those two dimensions land on top of ``c`` and of the link, because the foot
    of ``E`` on ``AD`` falls close to ``D``; they are drawn here instead, in the
    frame that defines them.

    Args:
        ax: Axes to place the inset inside.
        design: The mechanism, for ``c``, ``x_b`` and ``y_b``.
        colour: Colour of the two dimensions.
    """
    inset = ax.inset_axes((0.62, 0.68, 0.36, 0.30))
    inset.set_aspect("equal")
    inset.set_facecolor("white")
    corner = {"A": (0.0, 0.0), "D": (design.c, 0.0), "E": (design.x_b, design.y_b)}
    order = ("A", "D", "E", "A")
    inset.plot(
        [corner[name][0] for name in order],
        [corner[name][1] for name in order],
        color="#2c3e50",
        lw=2.0,
        zorder=4,
    )
    inset.plot(
        [design.x_b, design.x_b, 0.0],
        [0.0, design.y_b, 0.0],
        color="#95a5a6",
        lw=0.9,
        ls=":",
        zorder=3,
    )
    for name, (x, y) in corner.items():
        inset.plot(x, y, marker="o", ms=4.0, color="#2c3e50", zorder=5)
        inset.annotate(
            name,
            (x, y),
            textcoords="offset points",
            xytext=(5.0, 4.0),
            fontsize=8.5,
            color="#2c3e50",
        )
    span = max(abs(design.c), abs(design.x_b), abs(design.y_b))
    _dimension(
        inset,
        (0.0, 0.0),
        (design.x_b, 0.0),
        r"$x_b$",
        offset=0.16 * span,
        colour=colour,
        fontsize=9.0,
    )
    _dimension(
        inset,
        (design.x_b, 0.0),
        (design.x_b, design.y_b),
        r"$y_b$",
        offset=0.16 * span,
        colour=colour,
        fontsize=9.0,
    )
    _dimension(
        inset,
        (0.0, 0.0),
        (design.c, 0.0),
        r"$c$",
        offset=-0.20 * span,
        colour="#2c3e50",
        fontsize=9.0,
    )
    inset.set_title("trigonal link, in its own frame", fontsize=8.5, pad=3.0)
    inset.set_xticks([])
    inset.set_yticks([])
    inset.grid(False)
    for spine in inset.spines.values():
        spine.set_visible(True)
        spine.set_color("#bdc3c7")
    inset.margins(0.30)


def plot_variables(
    design: Design,
    theta_1: float = 45.0,
    spec: EngineSpec = DEFAULT_SPEC,
    analysis: Analysis | None = None,
    **kwargs: Any,
) -> Figure:
    """Draw the mechanism with all eleven design variables dimensioned.

    The design vector is a list of symbols until someone sees where each one
    sits on the linkage, so this is the figure the use case opens with.  Every
    entry of ``X`` appears exactly once: the four member lengths and two crank
    throws as dimensions along the members they measure, the inter-axle
    distance and its orientation at the shafts, the local coordinates of ``E``
    in the frame the trigonal link carries, the cylinder offset as a horizontal
    dimension from the ``R1`` axis, and the crank dephasing through the
    kinematic relation that defines it.

    Args:
        design: The mechanism to draw.
        theta_1: Crank angle to freeze the linkage at [deg].  The default opens
            the triangle enough that no two dimensions overlap.
        spec: Fixed engine data.
        analysis: A precomputed analysis at full resolution.
        **kwargs: Forwarded to :func:`exlink.model.analyse`.

    Returns:
        The figure.
    """
    from .animation import LINK_STYLE, _draw_static, _mechanism_limits

    solved = _solved(analysis or analyse(design, spec=spec, **kwargs))
    angle = math.radians(float(theta_1))
    frame = solve_kinematics(design, spec=spec, theta_1=np.array([angle]))
    point = {name: np.asarray(getattr(frame, name)[0], dtype=float) for name in JOINT_NAMES}
    xlim, ylim = _mechanism_limits(solved, spec)

    with plt.style.context(STYLE):
        figure, ax = plt.subplots(figsize=(8.4, 8.6))
        ax.set_aspect("equal")
        ax.set_xlim(xlim[0] - 18.0, xlim[1] + 18.0)
        ax.set_ylim(ylim[0] - 8.0, ylim[1] + 4.0)
        _draw_static(ax, solved, spec)
        for name, style in LINK_STYLE.items():
            polyline = frame.bodies[name][0]
            ax.plot(polyline[:, 0], polyline[:, 1], solid_capstyle="round", **style)

        crown = point["H"][1]
        half_bore = 0.5 * spec.bore
        ax.plot(
            [design.x_1 - half_bore, design.x_1 + half_bore],
            [crown, crown],
            color="#2c3e50",
            lw=7.0,
            solid_capstyle="butt",
        )

        for name, offset in JOINT_LABEL_OFFSET.items():
            ax.plot(*point[name], marker="o", ms=5.0, color="#2c3e50", zorder=10)
            ax.annotate(
                name,
                (float(point[name][0]), float(point[name][1])),
                textcoords="offset points",
                xytext=offset,
                fontsize=9.5,
                color="#2c3e50",
                zorder=10,
            )

        # -- the four lengths and the two crank throws ---------------------------
        _dimension(ax, point["R1"], point["Q"], r"$q_1$", offset=-8.0, colour="#c0392b")
        _dimension(ax, point["R2"], point["D"], r"$q_2$", offset=9.0, colour="#16a085")
        _dimension(ax, point["Q"], point["A"], r"$a$", offset=11.0, colour="#8e44ad")
        _dimension(ax, point["A"], point["D"], r"$c$", offset=13.0, colour="#2c3e50")
        _dimension(ax, point["E"], point["P"], r"$e$", offset=-11.0, colour="#d35400")
        _dimension(ax, point["R1"], point["R2"], r"$I$", offset=-15.0, colour="#7f8c8d")

        # -- E in the frame the trigonal link carries ----------------------------
        # Drawn in that frame rather than on the mechanism: the foot of E on
        # AD falls within a millimetre of D here, so the two dimensions would
        # lie on top of ``c`` and on the link itself.
        _trigonal_inset(ax, design)

        # -- the cylinder offset --------------------------------------------------
        datum = float(ylim[0]) + 0.06 * (float(ylim[1]) - float(ylim[0]))
        ax.plot([0.0, 0.0], [0.0, datum], color="#95a5a6", lw=0.8, ls=":", zorder=6)
        ax.plot(
            [design.x_1, design.x_1],
            [datum, crown],
            color="#95a5a6",
            lw=0.8,
            ls=":",
            zorder=6,
        )
        _dimension(ax, (0.0, datum), (design.x_1, datum), r"$x_1$", colour="#34495e")

        # -- the two angles -------------------------------------------------------
        _angle_arc(ax, point["R1"], design.theta_r_rad, r"$\theta_r$", radius=34.0)
        _angle_arc(ax, point["R1"], angle, r"$\theta_1$", radius=19.0, colour="#c0392b")
        _angle_arc(
            ax,
            point["R2"],
            float(frame.theta_2[0]),
            r"$\theta_2$",
            radius=22.0,
            colour="#16a085",
        )

        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title(
            r"$X = (a,\, c,\, I,\, x_b,\, y_b,\, x_1,\, e,\, q_1,\, q_2,"
            r"\, \theta_f,\, \theta_r)^\top$"
            "\n"
            r"drawn at $\theta_1$ = "
            f"{theta_1:.0f}"
            r"$^\circ$;"
            r" the dephasing enters through $\theta_2 = -2\theta_1 + \theta_f$"
        )
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
