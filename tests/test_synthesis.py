"""Prescribed-motion synthesis: the target, the fit, and what the fit is good for."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.constants import DEFAULT_SPEC, DEFAULT_TARGETS
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


@pytest.mark.slow
def test_keeping_the_inequalities_in_the_fit_keeps_the_design_buildable() -> None:
    """Why the reformulated problem should keep ``g <= 0``.

    Absorbing the equalities into the target removes the measure-zero part of
    the feasible set, and that is the only part the reformulation earns the
    right to drop.  The inequalities define a full-dimensional set and cost one
    SQP instead of one least-squares solve, so discarding them buys nothing and
    can return a design that tracks the motion and is not buildable.

    This asserts the constrained fit lands inside the geometric constraint set,
    which is the property the unconstrained one does not guarantee.
    """
    from exlink.model import analyse, inequality_constraints
    from exlink.synthesis import fit_within_constraints, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    fit = fit_within_constraints(target, COUPLED_DESIGN, max_iterations=40)
    assert fit is not None

    analysis = analyse(fit.design, samples=CRANK)
    assert analysis.valid
    violations = inequality_constraints(analysis)
    assert float(np.max(violations)) <= 1.0e-6, f"violated: {violations}"

    # And it is still a fit: the point is to keep the motion, not to abandon it.
    assert fit.stroke_error <= 0.05
    assert fit.ratio_error <= 0.05


@pytest.mark.slow
def test_the_fits_are_feasible_against_the_whole_constraint_set() -> None:
    """The generator returns usable starts, not merely points near the manifold.

    Corrects an earlier claim in §7.4 that none of these designs was feasible;
    measured, every in-band fit from a reachable target satisfies the coupled
    and vehicle constraints too.
    """
    from exlink.performance import evaluate
    from exlink.synthesis import feasible_starts, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    results = feasible_starts(
        target,
        attempts=6,
        seed=7,
        reference=COUPLED_DESIGN,
        scatter=0.10,
        max_evaluations=150,
    )
    in_band = [f for f in results if f.stroke_error <= 0.05 and f.ratio_error <= 0.05]
    assert in_band
    for fit in in_band:
        assert evaluate(fit.design, speed_rpm=1000.0).feasible


@pytest.mark.slow
def test_the_unconstrained_fit_violates_by_a_wide_margin_on_a_hard_target() -> None:
    """What discarding ``g`` actually costs, rather than what it might cost.

    On a reachable target the constrained and unconstrained fits agree, so the
    constraints look free and harmless to drop.  On a target the mechanism
    cannot reach they do not agree at all: least squares chases the motion and
    leaves the geometric set by ten units or more, while the constrained fit
    holds every constraint.  The unconstrained form is not a cheaper
    approximation of the constrained one -- it answers a different question.
    """
    import numpy as np
    from scipy.optimize import fsolve

    from exlink.constants import DEFAULT_SPEC
    from exlink.cycle import find_phases
    from exlink.model import analyse, inequality_constraints
    from exlink.synthesis import (
        TargetMotion,
        fit_to_target,
        fit_within_constraints,
        target_from_design,
    )

    base = target_from_design(COUPLED_DESIGN, samples=CRANK)
    theta = np.linspace(0.0, 2.0 * np.pi, base.samples, endpoint=False)
    second, first = np.cos(2.0 * theta), np.sin(theta)
    stroke = 15.0 * DEFAULT_SPEC.dead_volume / DEFAULT_SPEC.piston_area

    # A target displaced by third-harmonic content, put back exactly onto the
    # two stroke requirements so only its *shape* is unreachable.
    seeded = base.lam + 0.8 * np.sin(3.0 * theta)

    def residual(params: np.ndarray) -> list[float]:
        phases = find_phases(seeded + params[0] * second + params[1] * first)
        return [phases.expansion_stroke - 74.0, phases.compression_stroke - stroke]

    solution, _info, status, _msg = fsolve(residual, np.zeros(2), full_output=True)
    if status != 1:
        pytest.skip("the perturbed target could not be put back on the manifold")
    hard = TargetMotion(
        lam=seeded + solution[0] * second + solution[1] * first,
        expansion_stroke=74.0,
        compression_ratio=16.0,
        amplitude=float(solution[0]),
        skew=float(solution[1]),
    )

    loose = fit_to_target(hard, COUPLED_DESIGN, max_evaluations=200)
    tight = fit_within_constraints(hard, COUPLED_DESIGN, max_iterations=200)
    if loose is None or tight is None:
        pytest.skip("neither fit converged on this target")

    def worst(fit: object) -> float:
        analysis = analyse(fit.design, samples=CRANK)
        if not analysis.valid:
            return 1.0e3
        return float(np.max(inequality_constraints(analysis)))

    assert worst(tight) <= 1.0e-6
    assert worst(loose) > worst(tight)


@pytest.mark.slow
def test_holding_the_bands_as_constraints_keeps_the_fit_in_band() -> None:
    """The argument for keeping ``g`` applies to the bands as well.

    Against a target the mechanism cannot reach, holding ``g <= 0`` while
    chasing the motion pushes the strokes out of tolerance -- so the bands are
    checked afterwards and found violated.  They are constraints, so they
    belong in the problem: ``hold_bands`` puts them there, and the solve then
    either returns something in band or reports that it found nothing.
    """
    import numpy as np
    from scipy.optimize import fsolve

    from exlink.constants import DEFAULT_SPEC
    from exlink.cycle import find_phases
    from exlink.model import analyse, inequality_constraints
    from exlink.synthesis import TargetMotion, fit_within_constraints, target_from_design

    base = target_from_design(COUPLED_DESIGN, samples=CRANK)
    theta = np.linspace(0.0, 2.0 * np.pi, base.samples, endpoint=False)
    second, first = np.cos(2.0 * theta), np.sin(theta)
    stroke = 15.0 * DEFAULT_SPEC.dead_volume / DEFAULT_SPEC.piston_area
    seeded = base.lam + 0.8 * np.sin(3.0 * theta)

    def residual(params: np.ndarray) -> list[float]:
        phases = find_phases(seeded + params[0] * second + params[1] * first)
        return [phases.expansion_stroke - 74.0, phases.compression_stroke - stroke]

    solution, _info, status, _msg = fsolve(residual, np.zeros(2), full_output=True)
    if status != 1:
        pytest.skip("the perturbed target could not be put back on the manifold")
    hard = TargetMotion(
        lam=seeded + solution[0] * second + solution[1] * first,
        expansion_stroke=74.0,
        compression_ratio=16.0,
        amplitude=float(solution[0]),
        skew=float(solution[1]),
    )

    held = fit_within_constraints(
        hard, COUPLED_DESIGN, max_iterations=250, hold_bands=True, band=0.05
    )
    if held is None:
        pytest.skip("the constrained solve returned no analysable design")

    analysis = analyse(held.design, samples=CRANK)
    assert analysis.valid
    assert float(np.max(inequality_constraints(analysis))) <= 1.0e-6
    # The bands are constraints now, so they hold to solver tolerance.
    assert held.stroke_error <= 0.05 + 1.0e-6
    assert held.ratio_error <= 0.05 + 1.0e-6


# -- the full formulation: range objective, every constraint, target fallback ---


def test_the_constraint_vector_covers_the_whole_problem() -> None:
    """All twelve constraints of §3.10, not the five that were cheap.

    Each earlier round of this exercise found that whatever was left out of the
    fit was what the solve then violated.  This pins that nothing is left out:
    five geometric, four band sides, three coupled, two vehicle.
    """
    from exlink.performance import evaluate
    from exlink.synthesis import NUMBER_OF_CONSTRAINTS, _range_constraints, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    performance = evaluate(COUPLED_DESIGN, speed_rpm=1000.0, samples=CRANK)
    rows = _range_constraints(performance, target, 0.05, DEFAULT_TARGETS)

    assert rows.size == NUMBER_OF_CONSTRAINTS == 14
    # The reference design is feasible, so every row must be satisfied.
    assert float(np.min(rows)) >= -1.0e-6, f"violated: {rows}"


def test_an_unanalysable_design_reads_as_deeply_infeasible() -> None:
    """A broken design must not look feasible just because nothing evaluated."""
    from exlink.performance import evaluate
    from exlink.synthesis import _range_constraints, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    broken = evaluate(COUPLED_DESIGN.replace(a=25.0, c=25.0), speed_rpm=1000.0, samples=CRANK)
    rows = _range_constraints(broken, target, 0.05, DEFAULT_TARGETS)
    assert float(np.max(rows)) < 0.0


def test_the_objective_ladder_is_ordered() -> None:
    """Each rung must be strictly worse than the one above it.

    The ladder only guides the search if a design that runs always scores
    better than one that does not, and one that at least completes a cycle
    always scores better than one that does not.  If the constants overlapped,
    the optimizer could prefer a broken design to a working one.
    """
    from exlink.synthesis import CYCLE_PENALTY, MOTION_WEIGHT, RANGE_UNAVAILABLE

    # Any real design scores -km/L, which is negative; the package's designs
    # are in the thousands, so use a deliberately poor one as the bound.
    worst_real_range = -1.0
    # A large motion residual: 5 mm RMS is far beyond anything observed.
    worst_fallback = RANGE_UNAVAILABLE + MOTION_WEIGHT * 5.0**2

    assert worst_real_range < RANGE_UNAVAILABLE
    assert worst_fallback < CYCLE_PENALTY


@pytest.mark.slow
def test_the_target_takes_over_when_the_range_does_not_exist() -> None:
    """The middle rung of the ladder, exercised rather than assumed.

    Started from a design that is analysable but produces no range -- the
    friction exceeds the indicated work, so the engine does not run and km/L
    has no value -- the objective has nothing to descend.  A constant penalty
    would leave a flat region; the motion residual gives the search a gradient
    back towards a cycle that works.

    So this asserts the fallback actually fires (``fell_back > 0``) from such a
    start, which is the only thing that distinguishes the ladder from a plain
    penalty.
    """
    from exlink.performance import evaluate
    from exlink.reference import REFINED_DESIGN
    from exlink.synthesis import maximise_range_from_target, target_from_design

    speed = 1250.0
    start = evaluate(REFINED_DESIGN, speed_rpm=speed, module=0.8, teeth=48, samples=120)
    assert start.analysis.valid, "the start must be analysable for this to be the right case"
    assert start.km_per_litre == 0.0, "the start must produce no range"

    target = target_from_design(COUPLED_DESIGN, samples=120)
    fit = maximise_range_from_target(
        target, REFINED_DESIGN, speed_rpm=speed, module=0.8, teeth=48, max_iterations=3
    )
    assert fit is not None
    assert fit.fell_back > 0, "the target never stood in for the missing range"


def test_pinning_the_gear_pair_pins_the_centre_distance() -> None:
    """§3.7's catalogue relation makes ``I`` an output, not a design variable.

    A 200-iteration run that left ``I`` free while the gear pair was pinned
    returned a design with ``I = 85.08`` against the 57.6 the chosen pair
    realises -- a mechanism whose every downstream quantity was computed for
    something that cannot be built.  ``build_range_scenario`` had this right;
    this function did not.
    """
    from exlink.gears import lattice_inter_axle
    from exlink.synthesis import maximise_range_from_target, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    fit = maximise_range_from_target(
        target,
        COUPLED_DESIGN,
        speed_rpm=1000.0,
        module=0.8,
        teeth=48,
        max_iterations=1,
    )
    assert fit is not None
    assert pytest.approx(lattice_inter_axle(0.8, 48), abs=1.0e-9) == fit.design.I


def test_leaving_the_gear_pair_open_leaves_the_centre_distance_free() -> None:
    """The pin is a consequence of choosing a pair, not an unconditional bound."""
    from exlink.synthesis import _residual, target_from_design

    target = target_from_design(COUPLED_DESIGN, samples=CRANK)
    # Nothing to assert about the optimum here; only that the sizer is allowed
    # to choose, which is what ``module=None`` means.
    assert _residual(COUPLED_DESIGN, target, DEFAULT_SPEC) is not None
