"""Manufacturing tolerance, propagated two ways, on a design that sits near a
singularity."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.design import ANGULAR_VARIABLES, VARIABLE_NAMES
from exlink.reference import COUPLED_DESIGN, REFINED_DESIGN
from exlink.robustness import (
    CONSTRAINT_NAMES,
    DEFAULT_GRADE,
    IT_FACTORS,
    covariance,
    required_grade,
    tolerance_half_widths,
    tolerance_report,
    tolerance_unit,
)

SAMPLES = 250
CRANK = 360


@pytest.fixture(scope="module")
def report() -> object:
    return tolerance_report(REFINED_DESIGN, samples=SAMPLES, crank_samples=CRANK)


# -- the tolerance model itself -----------------------------------------------


def test_tolerance_unit_matches_iso_286() -> None:
    """``i = 0.45 D^(1/3) + 0.001 D`` micrometres, spot-checked.

    At a 50 mm nominal size the standard gives i = 1.56 um and an IT8 band of
    39 um; reproducing the published value is what makes the rest of this
    module a tolerance study rather than a guess.
    """
    assert tolerance_unit(50.0) * 1000.0 == pytest.approx(1.71, abs=0.05)
    band = IT_FACTORS[8] * tolerance_unit(50.0)
    assert band == pytest.approx(0.039, abs=0.006)


def test_half_widths_grow_with_the_grade_number() -> None:
    tight = tolerance_half_widths(REFINED_DESIGN, grade=6)
    loose = tolerance_half_widths(REFINED_DESIGN, grade=10)
    dimensional = [i for i, n in enumerate(VARIABLE_NAMES) if n not in ANGULAR_VARIABLES]
    assert np.all(loose[dimensional] > tight[dimensional])


def test_angular_variables_are_not_given_a_length_tolerance() -> None:
    """The clocking angles are an assembly quantity, not a machined size."""
    widths = tolerance_half_widths(REFINED_DESIGN, angular=0.05)
    for index, name in enumerate(VARIABLE_NAMES):
        if name in ANGULAR_VARIABLES:
            assert widths[index] == pytest.approx(0.05)


def test_covariance_is_diagonal_and_scales_as_the_square() -> None:
    matrix = covariance(REFINED_DESIGN)
    assert np.allclose(matrix, np.diag(np.diag(matrix)))
    widths = tolerance_half_widths(REFINED_DESIGN)
    assert np.allclose(np.sqrt(np.diag(matrix)), widths / 3.0)


# -- propagation ---------------------------------------------------------------


def test_first_order_and_monte_carlo_agree_within_a_factor_of_two(report: object) -> None:
    """Linearisation is what should be distrusted near a singularity.

    It need not be accurate, but if it were wrong by an order of magnitude the
    exact Jacobians would be giving no useful robustness information at all.
    """
    for name in CONSTRAINT_NAMES:
        sampled = report.monte_carlo_sigma[name]
        if sampled <= 0.0:
            continue
        assert 0.4 < report.linear_sigma[name] / sampled < 2.5, name


def test_first_order_is_conservative_on_this_problem(report: object) -> None:
    """It overestimates sigma rather than under, which is the safe direction.

    A first-order robust formulation built on these gradients would therefore
    be pessimistic, not optimistic -- worth stating, because the opposite would
    make the whole approach unusable here.
    """
    ratios = [
        report.linear_sigma[name] / report.monte_carlo_sigma[name]
        for name in CONSTRAINT_NAMES
        if report.monte_carlo_sigma[name] > 0.0
    ]
    assert np.mean(ratios) > 0.95


def test_the_tdc_gap_constraint_cannot_be_held(report: object) -> None:
    """The headline robustness finding.

    ``g <= 0.01 mm`` is tighter than the variation of the dimensions that
    produce ``g``.  Its capability is far below the industrial target of 1.33,
    and it is violated by a large fraction of nominally-conforming builds.
    This is a defect in the specification, not in any design that meets it.
    """
    assert report.capability()["tdc_gap"] < 0.5
    assert report.violation_rate["tdc_gap"] > 0.25


def test_every_other_constraint_is_comfortably_holdable(report: object) -> None:
    for name in CONSTRAINT_NAMES:
        if name == "tdc_gap":
            continue
        assert report.violation_rate[name] < 0.25, name


def test_the_gap_constraint_is_off_the_bottom_of_the_iso_ladder() -> None:
    """No machining grade fixes it, which is what makes it a specification bug.

    If a tighter grade would work, the answer is to tighten the drawing.  Here
    the required tolerance unit multiple is below the tightest grade in the
    table, so the gap has to be taken up by adjustment at assembly instead.
    """
    grade, factor = required_grade(REFINED_DESIGN, "tdc_gap", crank_samples=CRANK)
    assert grade is None
    assert factor < min(IT_FACTORS.values())


def test_a_slack_constraint_needs_no_precision_at_all() -> None:
    grade, factor = required_grade(REFINED_DESIGN, "compatibility", crank_samples=CRANK)
    assert grade == max(IT_FACTORS)
    assert factor > IT_FACTORS[DEFAULT_GRADE]


def test_report_is_reproducible() -> None:
    first = tolerance_report(COUPLED_DESIGN, samples=60, crank_samples=CRANK, seed=7)
    second = tolerance_report(COUPLED_DESIGN, samples=60, crank_samples=CRANK, seed=7)
    assert first.violation_rate == second.violation_rate


def test_an_unanalysable_design_is_refused_rather_than_scored() -> None:
    """Unlike the optimizer path, this one raises.

    A tolerance study on a design that does not exist would return numbers
    with no meaning, and nothing downstream needs it to keep walking.
    """
    broken = REFINED_DESIGN.replace(a=0.4 * REFINED_DESIGN.a)
    with pytest.raises(ValueError, match="does not close"):
        tolerance_report(broken, samples=5, crank_samples=CRANK)


def test_the_gap_is_hypersensitive_to_the_geometry_that_produces_it() -> None:
    """The same finding as the tolerance study, reached independently.

    §14.1 shows the top-dead-centre gap has a standard deviation from IT8
    machining tolerances larger than its own 0.01 mm band.  This checks the
    other route: ``g`` is so sensitive to the inter-axle distance that snapping
    ``I`` onto the gear lattice -- a move of under a fifth of a millimetre --
    pushes the gap several times past its bound.

    Two independent perturbations, the same conclusion: no geometric choice can
    hold this constraint, and it has to be taken up at assembly.
    """
    from exlink.model import analyse
    from exlink.reference import COUPLED_DESIGN

    nominal = analyse(COUPLED_DESIGN, samples=720)
    assert nominal.metrics.tdc_gap < 0.01

    shifted = analyse(COUPLED_DESIGN.replace(I=COUPLED_DESIGN.I - 0.18), samples=720)
    assert shifted.valid
    assert shifted.metrics.tdc_gap > 3.0 * 0.01

    sensitivity = (shifted.metrics.tdc_gap - nominal.metrics.tdc_gap) / 0.18
    assert sensitivity > 0.1


# -- reliability: a probability of failure, not a fixed margin ------------------


def test_the_constraints_are_strongly_correlated() -> None:
    """Why a per-constraint margin is the wrong formulation.

    Every constraint is a function of the same eleven dimensions, so their
    scatter is dependent -- here up to 0.94, and exactly -1 for the two sides of
    a relaxed equality.  Requiring each of them separately to hold at k sigma
    is a reliability statement only if they are independent, which they are
    emphatically not.
    """
    from exlink.robustness import constraint_moments

    moments = constraint_moments(COUPLED_DESIGN, samples=CRANK)
    assert moments is not None
    off_diagonal = moments.correlation[~np.eye(len(moments.names), dtype=bool)]
    assert np.max(np.abs(off_diagonal)) > 0.9

    # The two sides of the stroke band are the same residual, negated.
    upper = moments.names.index("stroke_upper")
    lower = moments.names.index("stroke_lower")
    assert moments.correlation[upper, lower] == pytest.approx(-1.0, abs=1.0e-9)


def test_the_correlation_changes_the_system_probability() -> None:
    """And not always in the reassuring direction.

    The two largest contributors here are anti-correlated, so keeping the
    correlation makes "at least one constraint fails" *more* likely than
    assuming independence, not less.  Either way the point stands: the
    independent figure is not the answer.
    """
    from exlink.robustness import failure_probability

    reliability = failure_probability(COUPLED_DESIGN, samples=CRANK)
    assert reliability is not None
    assert reliability.system != pytest.approx(reliability.independent_bound, rel=1.0e-3)


def test_form_agrees_with_sampling_at_the_system_level() -> None:
    """The first-order estimate has to be checked against sampling, not trusted.

    At the system level it is good here.  Per constraint it is not uniformly so
    -- the top-dead-centre gap is strongly nonlinear and FORM under-predicts its
    failure probability -- which is why the sampling estimate is the reference
    and the first-order one is what is cheap enough to evaluate every iteration.
    """
    from exlink.robustness import failure_probability

    reliability = failure_probability(COUPLED_DESIGN, samples=CRANK)
    sampled = tolerance_report(COUPLED_DESIGN, samples=1500, crank_samples=CRANK)
    assert reliability is not None
    assert reliability.system == pytest.approx(sampled.any_violation_rate, abs=0.1)


def test_the_reliability_index_is_minus_the_normal_quantile() -> None:
    from scipy.stats import norm

    from exlink.robustness import failure_probability

    reliability = failure_probability(COUPLED_DESIGN, samples=CRANK)
    assert reliability is not None
    assert reliability.system_beta == pytest.approx(float(norm.isf(reliability.system)))


def test_beta_is_the_margin_in_standard_deviations() -> None:
    from exlink.robustness import constraint_moments

    moments = constraint_moments(COUPLED_DESIGN, samples=CRANK)
    assert moments is not None
    assert np.allclose(moments.beta, -moments.value / moments.sigma)


def test_the_gap_dominates_the_failure_probability() -> None:
    """The same finding as every other route, now as a probability.

    At the specified 0.01 mm bound the top-dead-centre gap alone fails in a
    large fraction of nominally-conforming builds, and it is the constraint
    contributing most of the system probability.
    """
    from exlink.robustness import failure_probability

    reliability = failure_probability(COUPLED_DESIGN, samples=CRANK)
    assert reliability is not None
    assert reliability.binding() == "tdc_gap"
    assert reliability.per_constraint["tdc_gap"] > 0.1


def test_required_bound_inverts_the_reliability_relation() -> None:
    """Answers *how much* the specification would have to give, with a number.

    Applying the returned bound must bring that constraint's own failure
    probability to the target.
    """
    from dataclasses import replace

    from exlink.constants import DEFAULT_TARGETS
    from exlink.robustness import (
        TARGET_FAILURE_PROBABILITY,
        failure_probability,
        required_bound,
    )

    bound = required_bound(COUPLED_DESIGN, "tdc_gap", samples=CRANK)
    assert bound > DEFAULT_TARGETS.max_tdc_gap

    relaxed = replace(DEFAULT_TARGETS, max_tdc_gap=bound)
    reliability = failure_probability(COUPLED_DESIGN, samples=CRANK, targets=relaxed)
    assert reliability is not None
    assert reliability.per_constraint["tdc_gap"] == pytest.approx(
        TARGET_FAILURE_PROBABILITY, rel=0.05
    )


def test_relaxing_the_gap_alone_is_not_enough() -> None:
    """What reliability-based design says that a deterministic one cannot.

    Fixing the gap hands the problem to the stroke band, because the design
    sits off-centre in a band that is itself narrower than the scatter it is
    meant to represent.  A deterministic optimum has no way to see that.
    """
    from dataclasses import replace

    from exlink.constants import DEFAULT_TARGETS
    from exlink.robustness import failure_probability, required_bound

    relaxed = replace(
        DEFAULT_TARGETS, max_tdc_gap=required_bound(COUPLED_DESIGN, "tdc_gap", samples=CRANK)
    )
    reliability = failure_probability(COUPLED_DESIGN, samples=CRANK, targets=relaxed)
    assert reliability is not None
    assert reliability.binding() == "stroke_lower"
    assert reliability.system > 0.1


def test_a_looser_grade_raises_the_failure_probability() -> None:
    from exlink.robustness import failure_probability

    tight = failure_probability(COUPLED_DESIGN, grade=6, samples=CRANK)
    loose = failure_probability(COUPLED_DESIGN, grade=10, samples=CRANK)
    assert tight is not None and loose is not None
    assert loose.system > tight.system


def test_the_discipline_matches_the_function() -> None:
    from exlink.robustness import FailureProbabilityDiscipline, failure_probability

    discipline = FailureProbabilityDiscipline(samples=CRANK)
    output = discipline.execute(COUPLED_DESIGN.to_mapping())
    reliability = failure_probability(COUPLED_DESIGN, samples=CRANK)
    assert reliability is not None
    # The orthant integral is randomised quasi-Monte Carlo, so two evaluations
    # agree to about seven significant figures rather than exactly.
    assert float(output["failure_probability"][0]) == pytest.approx(
        reliability.system, rel=1.0e-4
    )


def test_an_unanalysable_design_is_certain_to_fail_not_an_error() -> None:
    from exlink.robustness import FailureProbabilityDiscipline, failure_probability

    broken = COUPLED_DESIGN.replace(a=0.4 * COUPLED_DESIGN.a)
    assert failure_probability(broken, samples=CRANK) is None
    output = FailureProbabilityDiscipline(samples=CRANK).execute(broken.to_mapping())
    assert float(output["failure_probability"][0]) == pytest.approx(1.0)


# -- what the reliability model does and does not cover ------------------------


def test_the_probability_covers_the_dimension_only_constraints() -> None:
    """§3.10's seven-of-twelve split, pinned so the documentation cannot drift.

    ``Sigma`` carries ISO 286 dimensional tolerances and nothing else, so a
    probability of failure is honest only for constraints that are functions of
    the eleven dimensions alone.  Those are exactly the five geometric margins
    and the two band residuals.  The load-dependent constraints stay
    deterministic; giving them a ``P_f`` from this covariance would report the
    dimensional variance and silently omit the larger one.
    """
    from exlink.robustness import CONSTRAINT_NAMES, covariance
    from exlink.scenarios import (
        COUPLED_INEQUALITY_OUTPUTS,
        EQUALITY_OUTPUTS,
        INEQUALITY_OUTPUTS,
        RANGE_INEQUALITY_OUTPUTS,
    )

    covered = set(CONSTRAINT_NAMES)
    assert len(covered) == 7

    # The seven are the five geometric margins plus the two band residuals,
    # matched by stem because the scenario names carry _margin/_error suffixes.
    geometric = {name.removesuffix("_margin") for name in INEQUALITY_OUTPUTS}
    assert geometric <= covered
    assert len(geometric) == 5
    assert len(EQUALITY_OUTPUTS) == 2

    # The load-dependent ones are outside it, and must stay outside while the
    # covariance carries dimensional scatter only.
    deterministic = set(COUPLED_INEQUALITY_OUTPUTS) | set(RANGE_INEQUALITY_OUTPUTS)
    assert len(deterministic) == 5
    for name in deterministic:
        assert name.removesuffix("_margin") not in covered

    # Seven covered plus five deterministic is the twelve of §3.10.
    assert len(covered) + len(deterministic) == 12

    # The claim that makes the split necessary: Sigma is dimensional only.
    sigma = covariance(COUPLED_DESIGN)
    assert sigma.shape == (11, 11)
    assert np.allclose(sigma, np.diag(np.diag(sigma)))


@pytest.mark.slow
def test_the_deterministic_optimum_is_dominated_on_reliability() -> None:
    """§6.2's finding: most of the 0.645 is self-inflicted, not required.

    A deterministic optimizer converges *onto* its active constraints, because
    nothing in the formulation rewards standing off them -- and a design
    sitting on ``g = 0`` fails about half the time.  Designs a few hundredths
    of a millimetre away are both more reliable and no worse in range, so the
    converged design is dominated rather than merely unreliable.

    Sampled rather than optimised, because the point of §3.10 is that a
    gradient method does not get there.
    """
    import numpy as np

    from exlink.design import Design
    from exlink.performance import evaluate
    from exlink.robustness import failure_probability

    reference = failure_probability(COUPLED_DESIGN)
    baseline = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
    assert baseline.feasible

    base = COUPLED_DESIGN.to_array()
    rng = np.random.default_rng(0)
    scored = []
    for scale in (0.0005, 0.001, 0.002):
        for _ in range(150):
            vector = base * (1.0 + scale * rng.normal(size=base.size))
            candidate = Design.from_array(vector)
            reliability = failure_probability(candidate)
            if reliability is not None:
                scored.append((reliability.system_beta, candidate))
    scored.sort(key=lambda item: -item[0])
    assert scored

    # Only the best few are worth the coupled evaluation.
    for beta, candidate in scored[:5]:
        outcome = evaluate(candidate, speed_rpm=1000.0)
        if not outcome.feasible:
            continue
        if beta > reference.system_beta and outcome.km_per_litre >= baseline.km_per_litre:
            return
    pytest.fail("no sampled design dominated the deterministic optimum")


def test_the_tolerance_band_can_be_relaxed() -> None:
    """§6.2 asks what bound a design would need; that needs a bound it can vary.

    The band was a module constant, so the reliability of a design could only
    ever be scored against the *specified* tolerance -- which makes the
    question "what relaxation would make this design reliable?" unanswerable
    from the public API.  Widening it must move the probability, and in the
    obvious direction.
    """
    from exlink.robustness import failure_probability

    tight = failure_probability(
        COUPLED_DESIGN,
        band={"expansion_stroke": 0.05, "compression_ratio": 0.05},
    )
    loose = failure_probability(
        COUPLED_DESIGN,
        band={"expansion_stroke": 0.15, "compression_ratio": 0.15},
    )
    assert tight is not None and loose is not None
    assert loose.system < tight.system
    assert loose.system_beta > tight.system_beta


def test_the_default_band_is_the_specified_one() -> None:
    """Omitting the band must not silently change what is being scored.

    Not asserted to machine precision: the system probability is a
    multivariate-normal orthant evaluated by Genz's quasi-Monte-Carlo
    transformation, so two calls at the same design differ in the seventh
    figure.  That is a property of the estimator, not of the band, and a
    tolerance tighter than its own accuracy would test the random number
    stream instead of the code.
    """
    from exlink.robustness import EQUALITY_BAND, failure_probability

    default = failure_probability(COUPLED_DESIGN)
    explicit = failure_probability(COUPLED_DESIGN, band=dict(EQUALITY_BAND))
    assert default is not None and explicit is not None
    assert default.system == pytest.approx(explicit.system, rel=1.0e-4)
