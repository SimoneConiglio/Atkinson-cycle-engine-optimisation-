"""Trace the efficiency-versus-size trade-off, two ways.

The report tried both: moving limits on the envelope with a single-objective
solver, and a multi-objective evolutionary algorithm. Moving limits are the
robust route on this problem -- see the note in ``exlink.scenarios``.

Takes about three minutes.

    python examples/04_pareto.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from exlink.cli import configure_logging
from exlink.design import Bounds
from exlink.plots import plot_pareto
from exlink.reference import REFINED_DESIGN
from exlink.scenarios import sweep_moving_limits

OUTDIR = Path("figures")


def main() -> None:
    configure_logging(verbose=False)
    OUTDIR.mkdir(exist_ok=True)

    # Walk an upper limit on H downwards, re-solving for maximum efficiency at
    # each step.  Each solve is an ordinary single-objective problem, warm-started
    # from the previous solution.
    limits = [240.0, 232.0, 224.0, 216.0, 208.0]
    outcomes = sweep_moving_limits(
        limits,
        bounds=Bounds.around(REFINED_DESIGN, relative=0.35, absolute_angle=40.0),
        initial=REFINED_DESIGN,
        samples=360,
        max_iter=40,
    )

    print(f"  {'H limit':>9}{'eta [%]':>10}{'H [mm]':>10}{'B [mm]':>10}  feasible")
    designs = []
    for limit, outcome in zip(limits, outcomes, strict=True):
        metrics = outcome.analysis.metrics
        print(
            f"  {limit:>9.0f}{100 * metrics.efficiency:>10.3f}"
            f"{metrics.height:>10.1f}{metrics.width:>10.1f}  {outcome.feasible}"
        )
        if outcome.feasible:
            designs.append(outcome.design)

    if designs:
        figure = plot_pareto(designs, samples=360)
        path = OUTDIR / "pareto.png"
        figure.savefig(path, dpi=140, bbox_inches="tight")
        print(f"\nwritten: {path}")


if __name__ == "__main__":
    main()
