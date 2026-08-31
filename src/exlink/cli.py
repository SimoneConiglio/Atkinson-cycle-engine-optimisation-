"""Command-line interface: ``exlink <command>``.

Commands
--------
``analyse``    Print every objective and constraint of a design.
``animate``    Write a GIF or MP4 of the mechanism turning.
``plot``       Write the motion / cycle / torque / mechanism figures.
``optimize``   Maximise efficiency subject to all the constraints.
``refine``     Polish a design with the augmented Lagrangian.
``pareto``     Approximate the Pareto front of ``(-eta, H, B)``.

A design is given either by name (``published``, ``refined``) or as a JSON file
written by ``--save``, so the commands chain::

    exlink optimize --save best.json
    exlink refine --design best.json --save final.json
    exlink animate --design final.json -o final.gif
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from .design import GLOBAL_BOUNDS, VARIABLE_DESCRIPTIONS, VARIABLE_NAMES, Bounds, Design
from .dynamics import DEFAULT_SPEED_RPM
from .model import analyse
from .reference import PUBLISHED_DESIGN

NAMED_DESIGNS = ("published", "refined")


def load_design(reference: str | None) -> Design:
    """Resolve a design from a name or a JSON file path.

    Args:
        reference: ``"published"``, ``"refined"``, a path to a JSON file, or
            ``None`` for the refined reference.

    Returns:
        The design.

    Raises:
        SystemExit: If the reference cannot be resolved.
    """
    from .reference import REFINED_DESIGN

    if reference is None or reference == "refined":
        return REFINED_DESIGN
    if reference == "published":
        return PUBLISHED_DESIGN
    path = Path(reference)
    if not path.is_file():
        raise SystemExit(
            f"no such design: {reference!r} (expected one of {NAMED_DESIGNS} or a JSON file)"
        )
    data = json.loads(path.read_text())
    if isinstance(data, list):
        return Design.from_array(data)
    return Design(**{name: float(data[name]) for name in VARIABLE_NAMES})


def save_design(design: Design, path: str | Path) -> Path:
    """Write a design to a JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {name: float(getattr(design, name)) for name in VARIABLE_NAMES}
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def _bounds_from(args: argparse.Namespace, design: Design | None) -> Bounds:
    """Pick the design box: global, or a window around a design."""
    if getattr(args, "local", 0.0):
        if design is None:
            raise SystemExit("--local needs a --design to centre the box on")
        return Bounds.around(design, relative=args.local)
    return GLOBAL_BOUNDS


# -- commands --------------------------------------------------------------------


def _cmd_analyse(args: argparse.Namespace) -> int:
    from .scenarios import format_analysis

    design = load_design(args.design)
    analysis = analyse(design, samples=args.samples)
    print(format_analysis(analysis, title=f"design: {args.design or 'refined'}"))
    if args.save:
        print(f"\nwritten: {save_design(design, args.save)}")
    return 0 if analysis.valid else 1


def _cmd_animate(args: argparse.Namespace) -> int:
    from .animation import animate, animate_dashboard, save

    design = load_design(args.design)
    builder = animate if args.plain else animate_dashboard
    animation = builder(design, frames=args.frames, interval=int(1000 / args.fps))
    path = save(animation, args.output, fps=args.fps, dpi=args.dpi)
    print(f"written: {path}")
    return 0


def _cmd_plot(args: argparse.Namespace) -> int:
    from . import plots

    design = load_design(args.design)
    analysis = analyse(design, samples=args.samples)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figures = {
        "motion": plots.plot_motion(design, analysis),
        "cycle": plots.plot_cycle(design, analysis),
        "torque": plots.plot_torque(design, analysis),
        "mechanism": plots.plot_mechanism(design, analysis=analysis),
        "overview": plots.plot_overview(design, analysis),
    }
    for name, figure in figures.items():
        path = outdir / f"{name}.png"
        figure.savefig(path, dpi=args.dpi, bbox_inches="tight")
        print(f"written: {path}")
    return 0


def _cmd_optimize(args: argparse.Namespace) -> int:
    from .scenarios import maximise_efficiency

    start = load_design(args.design) if args.design else None
    outcome = maximise_efficiency(
        algorithm=args.algorithm,
        bounds=_bounds_from(args, start),
        initial=start,
        max_iter=args.max_iter,
        samples=args.samples,
        max_height=args.max_height,
        max_width=args.max_width,
        seed=args.seed,
    )
    print(outcome.summary())
    print(f"\nevaluations: {outcome.n_evaluations}")
    if args.save:
        print(f"written: {save_design(outcome.design, args.save)}")
    return 0 if outcome.feasible else 1


def _cmd_refine(args: argparse.Namespace) -> int:
    from .scenarios import refine

    outcome = refine(
        load_design(args.design),
        algorithm=args.algorithm,
        relative=args.relative,
        max_iter=args.max_iter,
        samples=args.samples,
    )
    print(outcome.summary())
    print(f"\nevaluations: {outcome.n_evaluations}")
    if args.save:
        print(f"written: {save_design(outcome.design, args.save)}")
    return 0 if outcome.feasible else 1


def _cmd_size(args: argparse.Namespace) -> int:
    from .coupled import solve_for_design
    from .scenarios import format_coupled

    design = load_design(args.design)
    result = solve_for_design(
        design,
        speed_rpm=args.rpm,
        samples=args.samples,
        max_iterations=args.max_iterations,
        relaxation=args.relaxation,
    )
    print(format_coupled(result, f"sizing at {args.rpm:.0f} rpm"))
    if not result.converged:
        print("\nthe sizing loop did not converge: try --relaxation 0.5")
    elif result.saturated:
        print(
            "\nsome member hit the diameter ceiling: at this speed the mechanism"
            "\ncannot be built to carry the inertia its own mass creates."
        )
    if args.plot:
        from .plots import plot_sizing

        path = Path(args.plot)
        path.parent.mkdir(parents=True, exist_ok=True)
        plot_sizing(result).savefig(path, dpi=args.dpi, bbox_inches="tight")
        print(f"\nwritten: {path}")
    return 0 if result.feasible else 1


def _cmd_pareto(args: argparse.Namespace) -> int:
    from .plots import plot_pareto
    from .scenarios import local_pareto, pareto_front

    start = load_design(args.design) if args.design else None
    if args.local and start is not None:
        outcome = local_pareto(
            start,
            relative=args.local,
            algorithm=args.algorithm,
            pop_size=args.pop_size,
            max_gen=args.max_gen,
            samples=args.samples,
            seed=args.seed,
        )
    else:
        outcome = pareto_front(
            algorithm=args.algorithm,
            initial=start,
            pop_size=args.pop_size,
            max_gen=args.max_gen,
            samples=args.samples,
            seed=args.seed,
        )
    print(f"front size: {len(outcome.front)}   evaluations: {outcome.n_evaluations}")
    print()
    print(outcome.summary())
    if outcome.front:
        figure = plot_pareto(outcome.front, highlight=outcome.design, samples=args.samples)
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=args.dpi, bbox_inches="tight")
        print(f"\nwritten: {path}")
    if args.save:
        print(f"written: {save_design(outcome.design, args.save)}")
    return 0 if outcome.front else 1


# -- parser ----------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        prog="exlink",
        description=(
            "Optimization of the EX-link Atkinson-cycle engine mechanism "
            "of an extended-expansion (Atkinson) engine linkage."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="design variables:\n"
        + "\n".join(f"  {n:<9} {VARIABLE_DESCRIPTIONS[n]}" for n in VARIABLE_NAMES),
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="log GEMSEO output")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_design(sub: argparse.ArgumentParser, default_help: str = "") -> None:
        sub.add_argument(
            "-d",
            "--design",
            default=None,
            help=f"'published', 'refined' (default), or a JSON file{default_help}",
        )

    analyse_parser = subparsers.add_parser("analyse", help="print objectives and constraints")
    add_design(analyse_parser)
    analyse_parser.add_argument("--samples", type=int, default=720)
    analyse_parser.add_argument("--save", default=None, help="write the design to JSON")
    analyse_parser.set_defaults(func=_cmd_analyse)

    animate_parser = subparsers.add_parser("animate", help="write a GIF or MP4")
    add_design(animate_parser)
    animate_parser.add_argument("-o", "--output", default="exlink.gif")
    animate_parser.add_argument("--frames", type=int, default=180)
    animate_parser.add_argument("--fps", type=int, default=25)
    animate_parser.add_argument("--dpi", type=int, default=110)
    animate_parser.add_argument(
        "--plain", action="store_true", help="mechanism only, without the side plots"
    )
    animate_parser.set_defaults(func=_cmd_animate)

    plot_parser = subparsers.add_parser("plot", help="write the static figures")
    add_design(plot_parser)
    plot_parser.add_argument("-o", "--outdir", default="figures")
    plot_parser.add_argument("--samples", type=int, default=720)
    plot_parser.add_argument("--dpi", type=int, default=140)
    plot_parser.set_defaults(func=_cmd_plot)

    optimize_parser = subparsers.add_parser("optimize", help="maximise efficiency")
    add_design(optimize_parser, " to start from")
    optimize_parser.add_argument("-a", "--algorithm", default="DIFFERENTIAL_EVOLUTION")
    optimize_parser.add_argument("--max-iter", type=int, default=4000, dest="max_iter")
    optimize_parser.add_argument("--samples", type=int, default=360)
    optimize_parser.add_argument(
        "--max-height", type=float, default=float("inf"), dest="max_height"
    )
    optimize_parser.add_argument(
        "--max-width", type=float, default=float("inf"), dest="max_width"
    )
    optimize_parser.add_argument(
        "--local",
        type=float,
        default=0.0,
        help="search a box of this relative half-width around --design",
    )
    optimize_parser.add_argument("--seed", type=int, default=1)
    optimize_parser.add_argument("--save", default=None)
    optimize_parser.set_defaults(func=_cmd_optimize)

    refine_parser = subparsers.add_parser("refine", help="augmented Lagrangian polish")
    add_design(refine_parser)
    refine_parser.add_argument("-a", "--algorithm", default="Augmented_Lagrangian_order_0")
    refine_parser.add_argument("--relative", type=float, default=0.25)
    refine_parser.add_argument("--max-iter", type=int, default=200, dest="max_iter")
    refine_parser.add_argument("--samples", type=int, default=720)
    refine_parser.add_argument("--save", default=None)
    refine_parser.set_defaults(func=_cmd_refine)

    pareto_parser = subparsers.add_parser("pareto", help="approximate the Pareto front")
    add_design(pareto_parser, " to seed a local front")
    pareto_parser.add_argument("-a", "--algorithm", default="PYMOO_NSGA2")
    pareto_parser.add_argument("--pop-size", type=int, default=200, dest="pop_size")
    pareto_parser.add_argument("--max-gen", type=int, default=60, dest="max_gen")
    pareto_parser.add_argument("--samples", type=int, default=360)
    pareto_parser.add_argument(
        "--local",
        type=float,
        default=0.0,
        help="run in a box of this relative half-width around --design",
    )
    pareto_parser.add_argument("--seed", type=int, default=1)
    pareto_parser.add_argument("-o", "--output", default="figures/pareto.png")
    pareto_parser.add_argument("--dpi", type=int, default=140)
    pareto_parser.add_argument("--save", default=None)
    pareto_parser.set_defaults(func=_cmd_pareto)

    size_parser = subparsers.add_parser(
        "size", help="size the parts with inertia in the load path"
    )
    add_design(size_parser)
    size_parser.add_argument(
        "--rpm",
        type=float,
        default=DEFAULT_SPEED_RPM,
        help="crankshaft speed; the single strongest driver of the answer",
    )
    size_parser.add_argument("--samples", type=int, default=180)
    size_parser.add_argument("--max-iterations", type=int, default=400, dest="max_iterations")
    size_parser.add_argument(
        "--relaxation",
        type=float,
        default=1.0,
        help="under-relaxation in (0, 1]; below 1 damps a stiff loop",
    )
    size_parser.add_argument("--plot", default=None, help="write a sizing figure here")
    size_parser.add_argument("--dpi", type=int, default=140)
    size_parser.set_defaults(func=_cmd_size)

    return parser


def configure_logging(verbose: bool) -> None:
    """Quieten GEMSEO unless the user asked for its output.

    GEMSEO installs its own stream handler and sets the ``gemseo`` logger to
    INFO **when it is imported**, which overwrites anything configured before
    that.  Importing it here first means the level set below is the one that
    sticks; otherwise every solver iteration is printed.

    Args:
        verbose: Whether to let GEMSEO log at INFO.
    """
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(levelname)s %(message)s",
    )
    import gemseo  # noqa: F401  - imported for its logging side effect

    logging.getLogger("gemseo").setLevel(logging.INFO if verbose else logging.WARNING)


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``exlink`` console script."""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.verbose)
    handler: Any = args.func
    return int(handler(args))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
