"""Augmented Lagrangian, started from the historical baseline.

Started from a design that does not satisfy its own constraints when
re-analysed, the augmented Lagrangian lands on a genuinely feasible one at
essentially the same efficiency.

Takes about two minutes.

    python examples/03_optimize.py
"""

from __future__ import annotations

from exlink import PUBLISHED_DESIGN, Design, analyse
from exlink.cli import configure_logging
from exlink.scenarios import format_analysis, refine


def main() -> None:
    configure_logging(verbose=False)
    start = analyse(PUBLISHED_DESIGN, samples=720)
    print(f"start:  eta = {100 * start.metrics.efficiency:6.3f} %   feasible = False")

    design = PUBLISHED_DESIGN
    best: Design | None = None
    for step, relative in enumerate((0.25, 0.10, 0.05), start=1):
        outcome = refine(
            design,
            samples=480,
            max_iter=60,
            relative=relative,
            sub_algorithm_settings={"max_iter": 250},
        )
        design = outcome.design
        metrics = outcome.analysis.metrics
        print(
            f"pass {step}: eta = {100 * metrics.efficiency:6.3f} %   "
            f"H = {metrics.height:6.1f}   B = {metrics.width:6.1f}   "
            f"feasible = {outcome.feasible}"
        )
        if outcome.feasible:
            best = design

    # An augmented Lagrangian drives the active constraints onto their bounds,
    # so its final iterate often sits a whisker outside one of them -- here it
    # is g, which the optimizer pushes hard against 0.01 mm. Keep the last pass
    # that was actually feasible rather than the last pass. (The design shipped
    # as exlink.reference.REFINED_DESIGN was produced the same way, but with a
    # margin held on the active constraints so it stays feasible at any
    # crank-angle resolution.)
    if best is None:
        print("\nno pass reached a feasible design; try a larger budget")
        return

    print()
    print(format_analysis(analyse(best, samples=1440), "augmented Lagrangian result"))


if __name__ == "__main__":
    main()
