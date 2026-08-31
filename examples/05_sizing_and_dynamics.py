"""Size the parts with inertia in the load path, and watch the coupling bite.

The geometric problem stops before this: the masses are not known until the
parts have a shape. Restoring inertia closes a loop -- sections set masses, masses set
inertia loads, loads set sections -- which has to be solved rather than
sequenced.

Two things this prints that the quasi-static study cannot show:

* the mean torque is untouched by engine speed, while the peak bearing load
  grows as its square and the structural mass as its sixth power;
* the near-singular geometry that is optimal without inertia becomes the worst
  possible choice with it.

Takes about a minute.

    python examples/05_sizing_and_dynamics.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from exlink import analyse
from exlink.cli import configure_logging
from exlink.coupled import solve_for_design
from exlink.plots import plot_bearing_loads, plot_mass_vs_speed, plot_sizing
from exlink.reference import REFINED_DESIGN
from exlink.scenarios import format_coupled

OUTDIR = Path("figures")
SPEEDS = (0.0, 500.0, 1000.0, 1500.0, 2000.0)


def main() -> None:
    configure_logging(verbose=False)
    OUTDIR.mkdir(exist_ok=True)

    print("Sizing the reference design across engine speeds.")
    print("Efficiency is a mean-torque quantity, so it cannot move; the parts can.\n")
    print(
        f"  {'rpm':>6}{'mass [kg]':>12}{'peak bearing [N]':>19}"
        f"{'mean Mr [N.mm]':>17}{'  sweeps':>9}  state"
    )

    results = []
    for rpm in SPEEDS:
        result = solve_for_design(
            REFINED_DESIGN, speed_rpm=rpm, samples=360, max_iterations=400
        )
        results.append(result)
        state = (
            "ok" if result.feasible else ("RUNAWAY" if result.saturated else "not converged")
        )
        print(
            f"  {rpm:>6.0f}{result.total_mass_kg:>12.3f}"
            f"{result.peak_bearing_load:>19.0f}{result.loads.mean_torque:>17.1f}"
            f"{result.iterations:>9}  {state}"
        )

    print()
    print(format_coupled(results[2], "sizing at 1000 rpm"))

    # The baseline sits at W = 0.981, a hair from the singularity, because
    # that is where the quasi-static lever arm is longest.  Proximity to the
    # singularity is also what amplifies the accelerations, so with inertia in
    # the load path the same choice becomes the expensive one.
    print("\n\nBacking away from the transmission-angle singularity, at 1000 rpm:")
    print(
        f"  {'swing rod':>10}{'W':>9}{'eta [%]':>10}{'H [mm]':>9}"
        f"{'mass [kg]':>12}{'bearing [N]':>14}"
    )
    for scale in (1.00, 0.94, 0.88, 0.82):
        design = REFINED_DESIGN.replace(a=REFINED_DESIGN.a * scale)
        analysis = analyse(design, samples=360)
        if not analysis.valid:
            continue
        sized = solve_for_design(design, speed_rpm=1000.0, samples=360, max_iterations=400)
        print(
            f"  {'x' + format(scale, '.2f'):>10}{analysis.metrics.compatibility:>9.4f}"
            f"{100 * analysis.metrics.efficiency:>10.2f}{analysis.metrics.height:>9.1f}"
            f"{sized.total_mass_kg:>12.3f}{sized.peak_bearing_load:>14.0f}"
        )

    for name, figure in {
        "sizing": plot_sizing(results[2]),
        "bearing_loads": plot_bearing_loads(results[:3], [f"{s:.0f} rpm" for s in SPEEDS[:3]]),
        "mass_vs_speed": plot_mass_vs_speed(list(SPEEDS), results),
    }.items():
        path = OUTDIR / f"{name}.png"
        figure.savefig(path, dpi=140, bbox_inches="tight")
        print(f"written: {path}")


if __name__ == "__main__":
    main()
