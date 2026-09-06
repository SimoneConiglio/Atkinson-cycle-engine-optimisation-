"""Ask what the reliability figure looks like when it is not only dimensional.

Sections 6.2 and 6.4 price manufacturing scatter on the eleven dimensions, and
five of the thirteen constraints get no probability at all because they are not
functions of those eleven.  This widens the uncertain vector to seventeen --
adding material strength, stiffness, density, friction and the gas load -- and
asks two questions of the answer: does it change what was already priced, and
what does it find that was invisible?

Takes about three minutes.

    python examples/10_widened_uncertainty.py
"""

from __future__ import annotations

import dataclasses

from exlink.cli import configure_logging
from exlink.constants import DEFAULT_TARGETS
from exlink.performance import evaluate
from exlink.reference import RELIABLE_DESIGN
from exlink.robustness import constraint_moments, failure_probability, format_reliability
from exlink.uncertainty import (
    DEFAULT_SCATTER,
    PARAMETER_NAMES,
    format_shares,
    widened_moments,
    widened_reliability,
)

BAND = {"expansion_stroke": 0.15, "compression_ratio": 0.15}
"""The relaxed bounds §6.2 settles on."""

TARGETS = dataclasses.replace(DEFAULT_TARGETS, max_tdc_gap=0.1)
"""With the top-dead-centre gap widened, as §6.2 settles it."""

RUN = {"band": BAND, "targets": TARGETS, "speed_rpm": 1000.0}


def main() -> None:
    configure_logging(verbose=False)

    print("The six parameters and their scatter")
    print("=" * 62)
    for name, value in zip(PARAMETER_NAMES, DEFAULT_SCATTER.as_array(), strict=True):
        print(f"  {name:<22}CoV {value:.2f}")
    print()

    print("Thirteen constraints, seventeen uncertain inputs")
    print("=" * 62)
    moments, shares = widened_moments(RELIABLE_DESIGN, module=0.8, teeth=48, **RUN)
    print(format_shares(shares))
    print(
        "\n  The eight geometric constraints are untouched -- 100 % dimensions"
        "\n  to the last figure -- so nothing in section 6.2 or 6.4 moves.  The"
        "\n  five new ones are the mirror image: near enough 0 % dimensions,"
        "\n  which is exactly why the narrow model could not price them."
    )
    print()

    print("What was already priced is priced identically")
    print("=" * 62)
    narrow = constraint_moments(RELIABLE_DESIGN, targets=TARGETS, band=BAND)
    lookup = dict(zip(moments.names, moments.sigma, strict=True))
    print(f"  {'constraint':<15}{'sigma, narrow':>15}{'sigma, widened':>16}")
    for name, sigma in zip(narrow.names, narrow.sigma, strict=True):
        print(f"  {name:<15}{sigma:>15.5f}{lookup[name]:>16.5f}")
    print(
        "\n  Exactly, not approximately: the widened model takes these rows"
        "\n  from the same analytic Jacobian rather than differencing them."
    )
    print()

    print("What the widened model finds")
    print("=" * 62)
    print(format_reliability(widened_reliability(RELIABLE_DESIGN, module=0.8, teeth=48, **RUN)))
    print()

    print("The gear pair is the whole of it")
    print("=" * 62)
    for label, pair in (
        ("pinned, m=0.8 z=48", {"module": 0.8, "teeth": 48}),
        ("m=1.0 z=39", {"module": 1.0, "teeth": 39}),
    ):
        outcome = evaluate(RELIABLE_DESIGN, speed_rpm=1000.0, **pair)
        reliability = widened_reliability(RELIABLE_DESIGN, **pair, **RUN)
        worst = min(
            zip(reliability.moments.names, reliability.moments.beta, strict=True),
            key=lambda item: item[1],
        )
        print(
            f"  {label:<20}{outcome.km_per_litre:8.1f} km/L   "
            f"P_f = {reliability.system:9.3e}   binding: {worst[0]} at beta {worst[1]:.2f}"
        )
    narrow_reliability = failure_probability(RELIABLE_DESIGN, targets=TARGETS, band=BAND)
    print(f"\n  the dimensional-only figure §6.4 reports: {narrow_reliability.system:9.3e}")
    print(
        "\n  The pinned pair sits on its face-width limit, so whether it fits"
        "\n  is a coin flip -- and no amount of dimensional tolerance would"
        "\n  have said so, because the pair is chosen by the mixed-integer"
        "\n  master and its width is set by the torque, not by the linkage."
    )


if __name__ == "__main__":
    main()
