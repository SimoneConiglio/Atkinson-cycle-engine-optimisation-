"""Check the first-order reliability estimate against exact sampled builds.

Every probability in sections 6.2, 6.4 and 6.7 is FORM: the constraints are
linearised at the nominal design and the failure probability read off a
multivariate-normal orthant.  This runs the control -- random builds drawn from
the same covariance, each analysed in full, no linearisation anywhere -- and
then explains the discrepancy it finds, which is not the one anybody expected.

The default sample count here is small enough to run in a couple of minutes and
too small to resolve 1e-2 tightly; section 6.8's figures are from 150 000
draws.  Raise ``DRAWS`` to reproduce them, at about half an hour.

    python examples/11_sampling_check.py
"""

from __future__ import annotations

import dataclasses

import numpy as np

from exlink.cli import configure_logging
from exlink.constants import DEFAULT_SPEC, DEFAULT_TARGETS
from exlink.jacobian import kinematic_jacobian, metric_jacobian
from exlink.model import analyse
from exlink.reference import COUPLED_DESIGN, RELIABLE_DESIGN
from exlink.robustness import covariance, failure_probability
from exlink.sampling import sampled_reliability
from exlink.secondorder import differencing_scaling, format_scaling

DRAWS = 20_000
"""Builds to draw.  Section 6.8 used 150 000."""

BAND = {"expansion_stroke": 0.15, "compression_ratio": 0.15}
TARGETS = dataclasses.replace(DEFAULT_TARGETS, max_tdc_gap=0.1)


def main() -> None:
    configure_logging(verbose=False)

    print("FORM against sampling, at two designs")
    print("=" * 66)
    for label, design, targets, band in (
        ("COUPLED_DESIGN, spec as written", COUPLED_DESIGN, DEFAULT_TARGETS, None),
        ("RELIABLE_DESIGN, §5.2's bounds", RELIABLE_DESIGN, TARGETS, BAND),
    ):
        form = failure_probability(design, targets=targets, band=band)
        sampled = sampled_reliability(
            design, draws=DRAWS, targets=targets, band=band, importance=False
        )
        low, high = sampled.interval
        print(
            f"  {label:<34}FORM {form.system:9.3e}   sampled {sampled.system:9.3e}"
            f"   [{low:.3e}, {high:.3e}]"
        )
    print(
        "\n  Accurate at the first design and, before the fix below, optimistic"
        "\n  by a factor of seven at the second.  The difference between them is"
        "\n  the whole finding."
    )
    print()

    print("Why: the stroke is a maximum, and one design sits on the tie")
    print("=" * 66)
    print(f"  {'design':<20}{'TDC 1':>12}{'TDC 2':>12}{'gap':>12}")
    for label, design in (
        ("COUPLED_DESIGN", COUPLED_DESIGN),
        ("RELIABLE_DESIGN", RELIABLE_DESIGN),
    ):
        phases = analyse(design, samples=1440).require_solved().thermodynamics.phases
        first, second = phases.expansion_strokes
        print(f"  {label:<20}{first:>12.5f}{second:>12.5f}{abs(first - second):>12.2e}")
    print(
        "\n  The expansion stroke is measured from the higher of the two, and a"
        "\n  maximum of two smooth functions is not differentiable at the tie."
        "\n  Manufacturing scatter on these dimensions is about 7 microns, so"
        "\n  the second design straddles its tie in every build."
    )
    print()

    print("Which branch FORM used, and which one binds")
    print("=" * 66)
    analysis = analyse(RELIABLE_DESIGN, samples=360)
    kinematic = kinematic_jacobian(
        RELIABLE_DESIGN, analysis.require_solved().kinematics, DEFAULT_SPEC
    )
    rows = metric_jacobian(RELIABLE_DESIGN, analysis, kinematic, DEFAULT_SPEC)
    sigma_matrix = covariance(RELIABLE_DESIGN)
    phases = analysis.require_solved().thermodynamics.phases
    for index in (1, 2):
        gradient = np.asarray(rows[f"expansion_stroke_{index}"], dtype=float)
        sigma = float(np.sqrt(gradient @ sigma_matrix @ gradient))
        value = phases.expansion_strokes[index - 1] - 74.0 - BAND["expansion_stroke"]
        print(f"  branch {index}:  sigma {sigma:.5f}   beta {-value / sigma:6.3f}")
    print(
        "\n  FORM linearised the branch that attains the maximum, which is the"
        "\n  stiffer of the two.  Carrying both -- one extra Jacobian row -- is"
        "\n  what brings the estimate to within 1 % of sampling."
    )
    print()

    print("And why second derivatives would not have helped")
    print("=" * 66)
    print(format_scaling(differencing_scaling(RELIABLE_DESIGN)))


if __name__ == "__main__":
    main()
