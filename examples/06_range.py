"""Carry one design all the way to kilometres per litre, and see what it costs.

The geometric problem stops at the engine and has three objectives with no
exchange rate between them. This prints the chain that closes that gap: the
mass budget the envelope and the torque ripple actually buy, the friction the
side loads actually cost, and the range that prices both.

Takes about a minute.

    python examples/06_range.py
"""

from __future__ import annotations

from exlink.cli import configure_logging
from exlink.performance import evaluate, speed_sweep
from exlink.reference import COUPLED_DESIGN, REFINED_DESIGN
from exlink.slidercrank import (
    SliderCrank,
    evaluate_slidercrank,
    firing_frequency_sensitivity,
    optimise_slidercrank,
)


def main() -> None:
    configure_logging(verbose=False)

    print("Where the mass actually is, at 1000 rpm")
    print("=" * 38)
    outcome = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
    assert outcome.coupled is not None and outcome.friction is not None
    for name, mass in outcome.budget.kilograms().items():
        share = 100.0 * outcome.budget.shares()[name]
        print(f"  {name:<16}{mass:8.3f} kg{share:7.1f} %")
    print(f"  {'TOTAL':<16}{outcome.engine_mass_kg:8.3f} kg")
    print()
    print(
        f"  The sizing loop alone would have reported {outcome.coupled.total_mass_kg:.3f} kg."
    )
    print("  The flywheel follows from the torque ripple; the crankcase from H and B.")

    print()
    print("Where the work goes")
    print("=" * 19)
    for name, share in outcome.friction.breakdown_kj().items():
        print(f"  {name:<12}{100.0 * share:7.1f} % of indicated work")
    print(f"  indicated thermal efficiency  {100.0 * outcome.indicated_efficiency:.1f} %")
    print(f"  brake thermal efficiency      {100.0 * outcome.brake_efficiency:.1f} %")

    print()
    print("The speed trade the geometric problem could not see")
    print("=" * 50)
    print(f"  {'rpm':>6}{'mass kg':>10}{'flywheel':>10}{'eta_b':>8}{'km/L':>9}   feasible")
    for item in speed_sweep(COUPLED_DESIGN, speeds=(600.0, 800.0, 1000.0, 1250.0, 1500.0)):
        flywheel = 1000.0 * item.budget.items.get("flywheel", 0.0)
        print(
            f"  {item.speed_rpm:>6.0f}{item.engine_mass_kg:>10.2f}{flywheel:>10.2f}"
            f"{item.brake_efficiency:>8.3f}{item.km_per_litre:>9.0f}   {item.feasible}"
        )
    print()
    print("  Flywheel mass falls as 1/omega^2; structural mass climbs with omega^2.")
    print("  The optimum is interior, and the near-singular geometry cannot reach it:")
    near = evaluate(REFINED_DESIGN, speed_rpm=1250.0)
    print(f"    the W = 0.981 design at 1250 rpm -- {near.reason() or 'feasible'}")

    print()
    print("Is the linkage worth seven members?")
    print("=" * 35)
    # The baseline has to be optimised too, or the comparison measures the
    # optimization rather than the topology -- and the sign depends on it.
    hand_set = evaluate_slidercrank(SliderCrank.for_compression_ratio(16.0), 2000.0)
    best = optimise_slidercrank()
    otto = best.comparison
    sensitivity = firing_frequency_sensitivity(outcome)
    print(f"  {'mechanism':<38}{'eta_i':>8}{'eta_m':>8}{'km/L':>9}")
    print(
        f"  {'slider-crank, hand-set (r/l = 0.30)':<38}"
        f"{hand_set.indicated_efficiency:>8.3f}{hand_set.mechanical_efficiency:>8.3f}"
        f"{hand_set.km_per_litre:>9.0f}"
    )
    print(
        f"  {f'slider-crank, optimised (r/l = {best.mechanism.obliquity:.2f})':<38}"
        f"{otto.indicated_efficiency:>8.3f}{otto.mechanical_efficiency:>8.3f}"
        f"{otto.km_per_litre:>9.0f}"
    )
    print(
        f"  {'EX-link, 1 cycle per revolution':<38}"
        f"{outcome.indicated_efficiency:>8.3f}"
        f"{outcome.friction.mechanical_efficiency:>8.3f}"
        f"{sensitivity['km_per_litre']:>9.0f}"
    )
    # Twice the sliding and twice the rotation for the same indicated work,
    # so twice the friction -- not half the efficiency.
    doubled = 1.0 - 2.0 * (1.0 - outcome.friction.mechanical_efficiency)
    print(
        f"  {'EX-link, if it were a 4-stroke':<38}{outcome.indicated_efficiency:>8.3f}"
        f"{doubled:>8.3f}{sensitivity['km_per_litre_four_stroke']:>9.0f}"
    )
    print()
    print(f"  {'against':<38}{'as modelled':>14}{'as a 4-stroke':>16}")
    for label, baseline in (
        ("the hand-set baseline", hand_set.km_per_litre),
        ("the optimised baseline", otto.km_per_litre),
    ):
        modelled = 100.0 * (sensitivity["km_per_litre"] / baseline - 1.0)
        four_stroke = 100.0 * (sensitivity["km_per_litre_four_stroke"] / baseline - 1.0)
        print(f"  {label:<38}{modelled:>13.1f} %{four_stroke:>15.1f} %")
    print()
    print("  So the advantage is firing frequency, not extended expansion --")
    print("  and against a baseline that was optimised too, the sign changes.")


if __name__ == "__main__":
    main()
