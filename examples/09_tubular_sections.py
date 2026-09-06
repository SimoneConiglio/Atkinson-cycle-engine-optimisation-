"""Ask what the study's largest stated modelling conservatism is actually worth.

Section 7.2 listed solid round bars as the largest single conservatism in the
model and guessed 30 % of mass.  This measures it, and the answer is not the
one the guess implies: a bore is worth almost nothing at the design speed and a
great deal above it, because what it removes is not weight but *inertia in the
load path*.

Takes about two minutes.

    python examples/09_tubular_sections.py
"""

from __future__ import annotations

from exlink.cli import configure_logging
from exlink.performance import evaluate
from exlink.reference import RELIABLE_DESIGN
from exlink.sections import SOLID, TUBULAR, Section
from exlink.sizing import MEMBER_KINDS
from exlink.slidercrank import SliderCrank, evaluate_slidercrank

PINNED = {"module": 0.8, "teeth": 48}
SPEEDS = (600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0)
"""Analysis speeds; the crankshaft turns twice as fast."""


def main() -> None:
    configure_logging(verbose=False)

    print("Where the mass of this engine actually is")
    print("=" * 62)
    solid = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **PINNED)
    for name, share in sorted(solid.budget.shares().items(), key=lambda i: -i[1]):
        if share > 0.0:
            print(f"  {name:<15}{1000.0 * solid.budget.items[name]:7.2f} kg   {share:6.1%}")
    print(
        "\n  The seven members are the item a bore acts on, and they are under"
        "\n  2 % of the engine.  Any argument for tubes has to come from"
        "\n  somewhere other than their own weight."
    )
    print()

    print("What a bore buys, against speed")
    print("=" * 62)
    print(
        f"  {'crank rpm':>10}{'solid':>10}{'tubular':>10}{'gain':>8}"
        f"{'members g':>11}{'members g':>11}"
    )
    for speed in SPEEDS:
        bar = evaluate(RELIABLE_DESIGN, speed_rpm=speed, **PINNED)
        tube = evaluate(RELIABLE_DESIGN, speed_rpm=speed, section=TUBULAR, **PINNED)
        print(
            f"  {2.0 * speed:>10.0f}{bar.km_per_litre:>10.1f}{tube.km_per_litre:>10.1f}"
            f"{tube.km_per_litre / bar.km_per_litre - 1.0:>7.1%}"
            f"{1.0e6 * bar.budget.items['linkage']:>11.0f}"
            f"{1.0e6 * tube.budget.items['linkage']:>11.0f}"
        )
    print(
        "\n  Nothing at the bottom, everything at the top.  Above the design"
        "\n  speed the members are sized by the inertia of their own mass, so"
        "\n  taking mass out of them is compound interest; below it they are"
        "\n  sized by the gas load, which a bore does not touch."
    )
    print()

    print("Why the bore ratio has an interior optimum")
    print("=" * 62)
    print(f"  {'k':>6}{'floor mm':>10}{'km/L':>10}{'members g':>11}")
    for ratio in (0.0, 0.3, 0.5, 0.6, 0.7, 0.8):
        section = Section(bore_ratio=ratio) if ratio else SOLID
        floor = float(max(section.minimum_diameter(MEMBER_KINDS)))
        outcome = evaluate(RELIABLE_DESIGN, speed_rpm=1200.0, section=section, **PINNED)
        print(
            f"  {ratio:>6.2f}{floor:>10.1f}{outcome.km_per_litre:>10.1f}"
            f"{1.0e6 * outcome.budget.items['linkage']:>11.0f}"
        )
    print(
        "\n  The wall floor is what turns it round: a bored member may not be"
        f"\n  drawn below 2 t / (1 - k), which at k = 0.8 is {2 * 1.5 / 0.2:.0f} mm -- larger"
        "\n  than most of this linkage's members are when solid."
    )
    print()

    print("The same bore, given to the conventional engine too")
    print("=" * 62)
    for name, ratio, obliquity, rpm in (
        ("solid", 0.0, 0.1954, 2151.0),
        ("tubular", 0.6, 0.1755, 2232.0),
    ):
        mechanism = SliderCrank.for_compression_ratio(16.0, obliquity)
        section = Section(bore_ratio=ratio) if ratio else SOLID
        row = evaluate_slidercrank(mechanism, rpm, section=section)
        print(
            f"  {name:<9} r/l = {obliquity:.4f} at {rpm:.0f} rpm   "
            f"{row.km_per_litre:7.1f} km/L   {row.engine_mass:5.2f} kg"
        )
    print(
        "\n  Both re-optimised over their own freedoms, both bored under the"
        "\n  same rule.  The bore is worth +1.4 % to the conventional engine"
        "\n  and +0.9 % to the linkage, so the advantage of section 6.3 goes"
        "\n  from +17.6 % to +17.0 %."
    )


if __name__ == "__main__":
    main()
