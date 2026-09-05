"""Carry one design all the way to kilometres per litre, and see what it costs.

The geometric problem stops at the engine and has three objectives with no
exchange rate between them. This prints the chain that closes that gap: the
mass budget the envelope and the torque ripple actually buy, the friction the
side loads actually cost, and the range that prices both.

Takes a few minutes.

    python examples/06_range.py
"""

from __future__ import annotations

from exlink.cli import configure_logging
from exlink.constants import DEFAULT_TARGETS
from exlink.performance import evaluate, speed_sweep
from exlink.reference import COUPLED_DESIGN, REFINED_DESIGN
from exlink.slidercrank import (
    CAP_SPEED_RPM,
    SliderCrank,
    evaluate_slidercrank,
    optimise_slidercrank,
    optimise_slidercrank_to_specification,
    side_load_ratio,
    slidercrank_reliability,
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
    # Both engines take 720 deg of their output shaft per cycle, so equal
    # output speed is equal fuel per minute and the comparison needs no
    # correction.  The baseline has to be optimised too, or what is measured is
    # the optimization rather than the topology.
    best = optimise_slidercrank()
    otto = best.comparison
    matched = evaluate_slidercrank(best.mechanism, outcome.output_speed_rpm)
    print(
        f"  the linkage runs at {outcome.output_speed_rpm:.0f} rpm on its output shaft "
        f"({outcome.cycles_per_minute:.0f} cycles/min)"
    )
    print(f"  {'mechanism':<40}{'rpm':>7}{'eta_i':>8}{'eta_m':>8}{'km/L':>9}")
    print(
        f"  {'EX-link':<40}{outcome.output_speed_rpm:>7.0f}"
        f"{outcome.indicated_efficiency:>8.3f}"
        f"{outcome.friction.mechanical_efficiency:>8.3f}{outcome.km_per_litre:>9.0f}"
    )
    print(
        f"  {f'slider-crank, optimised (r/l = {best.mechanism.obliquity:.2f})':<40}"
        f"{best.speed_rpm:>7.0f}{otto.indicated_efficiency:>8.3f}"
        f"{otto.mechanical_efficiency:>8.3f}{otto.km_per_litre:>9.0f}"
    )
    print(
        f"  {'the same, at the linkage output speed':<40}"
        f"{matched.speed_rpm:>7.0f}{matched.indicated_efficiency:>8.3f}"
        f"{matched.mechanical_efficiency:>8.3f}{matched.km_per_litre:>9.0f}"
    )
    print()
    for label, baseline in (
        ("at each engine own best speed", otto.km_per_litre),
        ("at equal output speed", matched.km_per_litre),
    ):
        print(f"  {label:<40}{100.0 * (outcome.km_per_litre / baseline - 1.0):>+7.1f} %")

    print()
    print("Held to the EX-link's own rod-angle and side-load limits")
    print("=" * 55)
    # The baseline's own optimum misses both caps, so meeting the same
    # specification costs it range.  The cap is the quasi-static side-load
    # ratio, which is the quantity the EX-link's own cap is applied to.
    held = optimise_slidercrank_to_specification(grid=(13, 7), refinements=12)
    gamma = side_load_ratio(held.mechanism, CAP_SPEED_RPM)
    print(
        f"  r/l = {held.mechanism.obliquity:.5f}, {held.speed_rpm:.0f} rpm, "
        f"gamma = {gamma:.4f} against a cap of {DEFAULT_TARGETS.max_side_load}"
    )
    print(f"  range {held.comparison.km_per_litre:.0f} km/L against {otto.km_per_litre:.0f}")
    print(
        f"  EX-link (COUPLED_DESIGN) advantage over it: "
        f"{100.0 * (outcome.km_per_litre / held.comparison.km_per_litre - 1.0):.1f} %"
    )
    print()
    print("  That optimum sits on its cap, which is where reliability goes:")
    for obliquity in (held.mechanism.obliquity, held.mechanism.obliquity - 0.0004):
        mechanism = SliderCrank.for_compression_ratio(16.0, obliquity=obliquity)
        row = evaluate_slidercrank(mechanism, held.speed_rpm)
        reliability = slidercrank_reliability(mechanism, band=0.15)
        print(
            f"    r/l = {obliquity:.5f}  {row.km_per_litre:7.1f} km/L  "
            f"P_f = {reliability.system:.3e}  beta = {reliability.system_beta:+.2f}"
        )


if __name__ == "__main__":
    main()
