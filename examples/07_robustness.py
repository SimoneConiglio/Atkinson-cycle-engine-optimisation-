"""Ask whether the optimum survives the tolerances of the parts that make it.

The central finding of this study is about conditioning -- the mechanism sits
near a singularity -- so a deterministic optimum without a tolerance study
would be negligent. This prints one, propagated both from the exact Jacobians
and by sampling, and then asks which ISO grade would be needed to hold each
constraint.

Takes about a minute.

    python examples/07_robustness.py
"""

from __future__ import annotations

from exlink.cli import configure_logging
from exlink.formulations import (
    compare_formulations,
    coupling_curve,
    format_coupling,
    format_formulations,
)
from exlink.reference import COUPLED_DESIGN, REFINED_DESIGN
from exlink.robustness import (
    CONSTRAINT_NAMES,
    IT_FACTORS,
    format_report,
    required_grade,
    tolerance_report,
)


def main() -> None:
    configure_logging(verbose=False)

    for title, design in (
        ("near-singular (W = 0.981)", REFINED_DESIGN),
        ("backed off (W = 0.937)", COUPLED_DESIGN),
    ):
        print(title)
        print("-" * len(title))
        report = tolerance_report(design, samples=1000, crank_samples=360)
        print(format_report(report))
        print()
        print("  first-order sigma over sampled sigma (1.0 would be exactly linear):")
        ratios = report.linearity_error()
        print("   ", {name: round(ratios[name], 2) for name in CONSTRAINT_NAMES})
        print()

    print("Which ISO grade would each constraint need?")
    print("=" * 43)
    print(f"  tightest grade in the table is IT6, at {min(IT_FACTORS.values()):.0f}i")
    print()
    for name in ("tdc_gap", "expansion_stroke", "compatibility", "side_load"):
        grade, factor = required_grade(REFINED_DESIGN, name, crank_samples=360)
        verdict = f"IT{grade}" if grade is not None else "OFF THE LADDER"
        print(f"  {name:<20}needs {factor:>8.2f}i  ->  {verdict}")
    print()
    print("  The top-dead-centre gap needs a tolerance no machining grade offers.")
    print("  That is a defect in the specification, not in any design meeting it:")
    print("  the gap has to be taken up by adjustment at assembly, or the bound")
    print("  relaxed. No optimizer can fix it.")

    print()
    print(format_coupling(coupling_curve(COUPLED_DESIGN)))
    print()
    print("  At rest rho is exactly zero: with no inertia there is no path from")
    print("  mass to load, so the quasi-static problem is recovered and there is")
    print("  nothing to iterate. The gain grows with omega^2.")

    print()
    print(
        format_formulations(
            compare_formulations(
                initial=COUPLED_DESIGN, speed_rpm=1000.0, max_iter=15, relative=0.20
            )
        )
    )
    print()
    print("  IDF is not slower here, it is unavailable: the coupling variables are")
    print("  load histories, not scalars, so it would carry four orders of magnitude")
    print("  more variables than the problem has degrees of freedom.")


if __name__ == "__main__":
    main()
