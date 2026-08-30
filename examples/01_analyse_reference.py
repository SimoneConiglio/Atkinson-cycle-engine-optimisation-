"""Analyse the two reference designs and print their tables.

python examples/01_analyse_reference.py
"""

from __future__ import annotations

from exlink import PUBLISHED_DESIGN, PUBLISHED_METRICS, analyse
from exlink.reference import REFINED_DESIGN
from exlink.scenarios import format_analysis


def main() -> None:
    print(format_analysis(analyse(PUBLISHED_DESIGN, samples=1440), "2015 published table"))
    print()
    print(
        "Re-analysed as printed, the published table does not satisfy its own\n"
        "constraints -- most tellingly g, the quantity the optimizer drove to\n"
        "zero. See docs/theory.md section 8.\n"
    )

    refined = analyse(REFINED_DESIGN, samples=1440)
    print(format_analysis(refined, "refined by this framework"))
    print()

    print("against the properties the report reports:")
    print(f"  {'quantity':<22}{'here':>12}{'report':>12}")
    for name, reported in PUBLISHED_METRICS.items():
        if name == "torque_pressure_ratio":
            continue  # reported as a percentage with no stated normalisation
        print(f"  {name:<22}{getattr(refined.metrics, name):>12.4f}{reported:>12.4f}")


if __name__ == "__main__":
    main()
