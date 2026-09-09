"""Ask the optimizer for the mechanism, instead of giving it one.

Section 5.4 measures the EXlink topology against a slider-crank and finds it
worth 17.6 %.  Two topologies establish a contrast and not a trend, and choosing
a third by hand would move the arbitrariness rather than remove it.  This runs
the second family of section 2.2 instead -- the spring-connected model of Kim et
al. (2007), with the gear-linkage extension of Yim et al. (2019) -- on the
question the engine actually poses: turn the circular motion of one input shaft
into the piston motion an extended-expansion cycle needs, with two equal top
dead centres and two different bottom ones.

Nothing about the linkage is assumed.  A design domain is filled with candidate
members, every candidate is a spring whose stiffness is a design variable, and
the stiffnesses are penalised towards rigid or absent.  What comes out is a
topology.

One start takes twenty minutes on the bars-only domain and around an hour on
the geared one.  ``STARTS`` below is small for that reason; sections 5.6 and 5.7
report six each.

    python examples/12_topology_synthesis.py
"""

from __future__ import annotations

import time

import numpy as np

from exlink.cli import configure_logging
from exlink.synthesis import target_motion
from exlink.topology import (
    GEAR_CATALOGUE,
    assess,
    engine_ground_structure,
    format_assessment,
    geared_ground_structure,
    random_start,
    synthesise_one,
)

DOMAIN = "geared"
"""Which candidate set to search.

``"bars"`` is section 5.6's -- binary members only, and the 2:1 relation handed
to it.  ``"geared"`` is section 5.7's, which adds three-cornered bodies as
single candidates and lets the search choose the gear ratio, or no gear at all.
"""

STARTS = 2
"""Random starts to run.  Each is independent; sections 5.6 and 5.7 used six."""

SEARCH_SAMPLES = 48
"""Input angles per equilibrium sweep during the search."""

REPORT_SAMPLES = 720
"""Input angles the answer is finally judged on."""


def main() -> None:
    """Synthesise linkages from random starts and report what they realise."""
    configure_logging()
    build = geared_ground_structure if DOMAIN == "geared" else engine_ground_structure
    structure, layout = build()
    target = target_motion(samples=SEARCH_SAMPLES)

    bars = len(structure.members)
    print(f"the ground structure ({DOMAIN})")
    print(f"  nodes    {' '.join(structure.names)}")
    print(
        f"  elements {bars} bars, {structure.n_members - bars} bodies, "
        f"{len(structure.gears)} gear pairs"
    )
    if structure.gears:
        print(f"  ratios   {', '.join(f'{r:g}:1' for r in GEAR_CATALOGUE)}")
    print(
        f"  search   {layout.size} design variables, {structure.n_reduced} unknowns per angle"
    )
    print(
        f"  target   STE {target.expansion_stroke:.1f} mm, "
        f"epsilon {target.compression_ratio:.1f}, "
        f"first harmonic {abs(target.skew):.2f} mm, "
        f"second {abs(target.amplitude):.2f} mm"
    )
    print()

    rng = np.random.default_rng(3)
    for start in range(STARTS):
        began = time.time()
        candidate = synthesise_one(
            layout, random_start(layout, rng), target.lam, samples=SEARCH_SAMPLES
        )
        verdict = assess(layout, candidate.x, target.lam, samples=REPORT_SAMPLES)
        print(f"start {start}  [{time.time() - began:.0f} s]")
        print(format_assessment(verdict))
        print()


if __name__ == "__main__":
    main()
