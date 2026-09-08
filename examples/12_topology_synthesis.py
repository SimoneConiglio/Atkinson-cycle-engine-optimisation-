"""Ask the optimizer for the mechanism, instead of giving it one.

Section 5.4 measures the EXlink topology against a slider-crank and finds it
worth 17.6 %.  Two topologies establish a contrast and not a trend, and choosing
a third by hand would move the arbitrariness rather than remove it.  This runs
the second family of section 2.2 instead -- the spring-connected model of Kim et
al. (2007) -- on the question the engine actually poses: turn the circular
motion of one input shaft into the piston motion an extended-expansion cycle
needs, with two equal top dead centres and two different bottom ones.

Nothing about the linkage is assumed.  A design domain is filled with candidate
members, every candidate is a spring whose stiffness is a design variable, and
the stiffnesses are penalised towards rigid or absent.  What comes out is a
topology.

One start takes twenty minutes or so.  ``STARTS`` below is small for that
reason; section 5.6 reports six.

    python examples/12_topology_synthesis.py
"""

from __future__ import annotations

import time

import numpy as np

from exlink.cli import configure_logging
from exlink.synthesis import target_motion
from exlink.topology import (
    assess,
    engine_ground_structure,
    format_assessment,
    random_start,
    synthesise_one,
)

STARTS = 2
"""Random starts to run.  Each is independent; section 5.6 used six."""

SEARCH_SAMPLES = 48
"""Input angles per equilibrium sweep during the search."""

REPORT_SAMPLES = 720
"""Input angles the answer is finally judged on."""


def main() -> None:
    """Synthesise linkages from random starts and report what they realise."""
    configure_logging()
    structure, layout = engine_ground_structure()
    target = target_motion(samples=SEARCH_SAMPLES)

    print("the ground structure")
    print(f"  nodes    {' '.join(structure.names)}")
    print(f"  members  {structure.n_members} candidates, {layout.size} design variables")
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
