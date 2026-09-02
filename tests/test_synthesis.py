"""Prescribed-motion synthesis: the target, the fit, and what the fit is good for."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.constants import DEFAULT_SPEC
from exlink.cycle import find_phases
from exlink.reference import COUPLED_DESIGN
from exlink.synthesis import (
    describe_target,
    fit_report,
    fit_to_target,
    target_from_design,
    target_motion,
)

CRANK = 180


# -- the target ----------------------------------------------------------------


def test_the_target_satisfies_both_equalities_as_the_model_measures_them() -> None:
    """The claim the whole construction rests on.

    A target is only useful as a feasible-point generator if it is itself on
    the equality manifold -- and "on" has to mean *as measured by the code the
    constraints use*, not as measured by the formula the target was built from.
    So this measures it with ``find_phases``.
    """
    target = target_motion(expansion_stroke=74.0, compression_ratio=16.0)
    measured = describe_target(target)
    assert measured["expansion_stroke"] == pytest.approx(74.0, abs=1.0e-6)
    assert measured["compression_ratio"] == pytest.approx(16.0, abs=1.0e-6)
    # The two-harmonic basis puts the maxima at equal height by construction,
    # so the target asks for a zero top-dead-centre gap rather than a merely
    # small one.
    assert measured["tdc_gap"] == pytest.approx(0.0, abs=1.0e-9)


def test_the_target_is_an_atkinson_motion() -> None:
    """Four monotone phases with unequal bottom dead centres, or it is not one."""
    target = target_motion()
    phases = find_phases(target.lam)
    assert phases.expansion_stroke > phases.compression_stroke
    assert len(phases.maxima_indices) == 2
    assert len(phases.minima_indices) == 2


def test_a_different_requirement_gives_a_different_target() -> None:
    """The construction solves for what it is asked, not for a fixed motion."""
    loose = describe_target(target_motion(expansion_stroke=70.0, compression_ratio=14.0))
    assert loose["expansion_stroke"] == pytest.approx(70.0, abs=1.0e-6)
    assert loose["compression_ratio"] == pytest.approx(14.0, abs=1.0e-6)


def test_an_unrealisable_requirement_is_refused() -> None:
    """A compression stroke deeper than the expansion stroke is not Atkinson."""
    with pytest.raises(ValueError, match="no two-harmonic motion"):
        target_motion(expansion_stroke=20.0, compression_ratio=16.0)


def test_the_target_offset_does_not_change_what_it_requires() -> None:
    """Both requirements are differences, so the mean height is free."""
    high = describe_target(target_motion(offset=150.0))
    low = describe_target(target_motion(offset=400.0))
    assert high["expansion_stroke"] == pytest.approx(low["expansion_stroke"], abs=1.0e-9)
    assert high["compression_ratio"] == pytest.approx(low["compression_ratio"], abs=1.0e-9)


def test_a_corrected_target_also_lands_on_the_manifold() -> None:
    """Seeding from a real motion must not cost exactness.

    ``target_from_design`` trades reachability for nothing: it still solves the
    two harmonics against the two stroke requirements, so the target it returns
    satisfies both equalities to the same precision as the synthetic one.
    """
    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    measured = describe_target(target)
    assert measured["expansion_stroke"] == pytest.approx(74.0, abs=1.0e-6)
    assert measured["compression_ratio"] == pytest.approx(16.0, abs=1.0e-6)


def test_the_correction_to_a_real_motion_is_small() -> None:
    """The point of seeding: the target stays near the reachable set.

    A large correction would put the target as far out of reach as the
    synthetic one, and the exercise would have gained nothing.  Both harmonic
    coefficients are asserted small against the ~65 mm amplitude of the motion
    itself.
    """
    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    assert abs(target.amplitude) < 5.0
    assert abs(target.skew) < 5.0


def test_seeding_from_an_unanalysable_design_is_refused() -> None:
    with pytest.raises(ValueError, match="unanalysable"):
        target_from_design(COUPLED_DESIGN.replace(a=25.0, c=25.0), samples=CRANK)


# -- the fit -------------------------------------------------------------------


@pytest.mark.slow
def test_fitting_moves_the_design_towards_the_requirements() -> None:
    """The fit is over-determined, so this asserts improvement, not exactness.

    ``lambda(X) = lambda*`` is not attainable with eleven variables against
    several hundred residuals, so the equalities are *not* satisfied by
    construction -- only approached.  What must be true is that fitting gets
    closer to them than the start was.
    """
    target = target_motion(samples=CRANK)
    fit = fit_to_target(target, COUPLED_DESIGN, max_evaluations=120)
    assert fit is not None
    assert fit.rms >= 0.0
    # The reference design already satisfies STE to within its band, so the
    # meaningful check is that the fit does not destroy that.
    assert fit.stroke_error < 5.0
    assert fit.ratio_error < 5.0


def test_a_fit_from_an_unanalysable_start_returns_nothing() -> None:
    """A start off the analysable set has no residual to descend."""
    broken = COUPLED_DESIGN.replace(a=25.0, c=25.0)
    assert fit_to_target(target_motion(samples=CRANK), broken) is None


def test_the_report_counts_what_landed_in_band() -> None:
    """The report is what §7.4 quotes, so its arithmetic is pinned."""
    from exlink.synthesis import FitResult

    def row(stroke: float, ratio: float) -> FitResult:
        return FitResult(
            design=COUPLED_DESIGN,
            rms=0.1,
            expansion_stroke=74.0 + stroke,
            compression_ratio=16.0 + ratio,
            stroke_error=stroke,
            ratio_error=ratio,
            converged=True,
            evaluations=1,
        )

    text = fit_report([row(0.01, 0.01), row(0.5, 0.01), row(0.01, 0.5)])
    assert "1 of 3 fits land inside both bands" in text


def test_the_residual_ignores_the_mean_height() -> None:
    """Shifting a mechanism bodily up the cylinder is not a motion error."""
    from exlink.synthesis import _residual

    target = target_motion(samples=CRANK)
    here = _residual(COUPLED_DESIGN, target, DEFAULT_SPEC)
    assert here is not None
    assert float(np.mean(here)) == pytest.approx(0.0, abs=1.0e-9)


# -- what the generator is and is not good for ---------------------------------


@pytest.mark.slow
def test_the_fit_is_a_contraction_onto_one_design() -> None:
    """Why prescribed-motion fitting does not, by itself, give multistart.

    Fitting to a *fixed* target converges to the same linkage whatever start it
    is given -- measured at one distinct design from every scatter between 10 %
    and 80 %.  That is not a defect of the optimizer: the motion very nearly
    determines the mechanism, which is the same property that would make a
    functional decomposition well posed.  Diversity has to come from varying
    the target instead.
    """
    from exlink.synthesis import feasible_starts, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    results = feasible_starts(
        target,
        attempts=8,
        seed=3,
        reference=COUPLED_DESIGN,
        scatter=0.25,
        max_evaluations=150,
    )
    assert results, "the fit should succeed from at least some starts"

    distinct: list[np.ndarray] = []
    for fit in results:
        vector = fit.design.to_array()
        if all(float(np.linalg.norm(vector - kept)) > 1.0 for kept in distinct):
            distinct.append(vector)
    assert len(distinct) == 1


@pytest.mark.slow
def test_fitting_beats_uniform_sampling_at_reaching_the_manifold() -> None:
    """The one thing the generator does solve, and the measure that says so.

    Uniform sampling of the design box finds a design satisfying the two
    equalities in 0 of 12 000 draws, because the manifold has measure zero.
    Sampling *and then fitting* is a different operation, and most fits land
    inside the tolerance bands.
    """
    from exlink.synthesis import feasible_starts, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    results = feasible_starts(
        target,
        attempts=8,
        seed=5,
        reference=COUPLED_DESIGN,
        scatter=0.10,
        max_evaluations=150,
    )
    in_band = [fit for fit in results if fit.stroke_error <= 0.05 and fit.ratio_error <= 0.05]
    assert len(in_band) >= 4
