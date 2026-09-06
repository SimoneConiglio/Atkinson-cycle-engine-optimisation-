"""Ask whether the result is an artefact of the speed it was measured at.

Every figure in sections 6.1 to 6.4 is reported at one crankshaft speed. The
test is to score the same designs over a schedule of speeds -- and, crucially,
to score each of them as *one engine*, sized once for the whole schedule rather
than re-sized at every point.

Takes about a minute.

    python examples/08_drive_cycle.py
"""

from __future__ import annotations

from exlink.cli import configure_logging
from exlink.drivecycle import (
    STANDARD_CYCLE,
    format_cycle,
    score_cycle,
    score_slidercrank_cycle,
)
from exlink.performance import evaluate
from exlink.reference import RELIABLE_DESIGN
from exlink.slidercrank import SliderCrank, evaluate_slidercrank

BASELINE_OBLIQUITY = 0.195
"""The r/l of section 6.3's optimised conventional engine."""

BASELINE_RPM = 2151.0
"""The speed it was optimised at [rev/min]."""

PINNED = {"module": 0.8, "teeth": 48}
"""The gear pair every run in the study pins."""


def main() -> None:
    configure_logging(verbose=False)

    print("The EX-link over a schedule of speeds")
    print("=" * 55)
    linkage = score_cycle(RELIABLE_DESIGN, STANDARD_CYCLE, **PINNED)
    print(format_cycle(linkage))
    point = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **PINNED)
    print(
        f"\n  against {point.km_per_litre:.1f} km/L at the single design point: "
        f"{100.0 * (linkage.km_per_litre / point.km_per_litre - 1.0):+.1f} %"
    )
    print()

    print("The conventional engine over the same schedule")
    print("=" * 55)
    mechanism = SliderCrank.for_compression_ratio(16.0, BASELINE_OBLIQUITY)
    baseline = score_slidercrank_cycle(mechanism, STANDARD_CYCLE)
    print(format_cycle(baseline))
    row = evaluate_slidercrank(mechanism, BASELINE_RPM)
    print(
        f"\n  against {row.km_per_litre:.1f} km/L at its own single point: "
        f"{100.0 * (baseline.km_per_litre / row.km_per_litre - 1.0):+.1f} %"
    )
    print()

    print("Where the schedule's cost goes")
    print("=" * 55)
    for name, scored in (("EX-link", linkage), ("slider-crank", baseline)):
        structure = scored.engine_mass_kg - scored.flywheel_kg
        share = scored.flywheel_kg / scored.engine_mass_kg
        print(
            f"  {name:<13} structure {structure:5.2f} kg   "
            f"flywheel {scored.flywheel_kg:5.2f} kg ({share:.0%} of the engine)"
        )
    print(
        "\n  The lighter mechanism is the heavier engine: a schedule makes the"
        "\n  flywheel most of the mass, and the flywheel is what a flat torque"
        "\n  curve buys."
    )
    print()

    print("The advantage, at a point and over the schedule")
    print("=" * 55)
    print(
        f"  single point   {point.km_per_litre:7.1f} vs {row.km_per_litre:7.1f} km/L   "
        f"{100.0 * (point.km_per_litre / row.km_per_litre - 1.0):+.1f} %"
    )
    print(
        f"  schedule       {linkage.km_per_litre:7.1f} vs {baseline.km_per_litre:7.1f} km/L   "
        f"{100.0 * (linkage.km_per_litre / baseline.km_per_litre - 1.0):+.1f} %"
    )
    print(
        "\n  Neither design was tuned for the schedule, so the second figure is"
        "\n  a lower bound on what a cycle-aware EX-link would reach."
    )


if __name__ == "__main__":
    main()
