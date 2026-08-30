"""Reproduce the report's final step: augmented Lagrangian from the published table.

Started from a design that does not satisfy its own constraints when
re-analysed, the augmented Lagrangian lands on a genuinely feasible one at
essentially the reported efficiency.

    python examples/03_optimize.py
"""

from __future__ import annotations

from exlink import PUBLISHED_DESIGN, analyse
from exlink.scenarios import format_analysis, refine


def main() -> None:
    start = analyse(PUBLISHED_DESIGN, samples=720)
    print(f"start:  eta = {100 * start.metrics.efficiency:6.3f} %   feasible = False")

    design = PUBLISHED_DESIGN
    for step, relative in enumerate((0.25, 0.10, 0.05), start=1):
        outcome = refine(
            design,
            samples=720,
            max_iter=120,
            relative=relative,
            sub_algorithm_settings={"max_iter": 400},
        )
        design = outcome.design
        metrics = outcome.analysis.metrics
        print(
            f"pass {step}: eta = {100 * metrics.efficiency:6.3f} %   "
            f"H = {metrics.height:6.1f}   B = {metrics.width:6.1f}   "
            f"feasible = {outcome.feasible}"
        )

    print()
    print(format_analysis(analyse(design, samples=1440), "augmented Lagrangian result"))


if __name__ == "__main__":
    main()
