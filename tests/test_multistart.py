"""Restarting a gradient solve on a feasible set of measure zero."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.design import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design
from exlink.model import analyse
from exlink.reference import GRADIENT_DESIGN, REFINED_DESIGN
from exlink.scenarios import is_feasible, multistart


def test_uniform_sampling_never_finds_a_feasible_point() -> None:
    """Why GEMSEO's MultiStart cannot be used here, as a property rather than a hunch.

    The feasible set contains two *equality* constraints, so it is a
    codimension-two manifold with measure zero, and uniform sampling hits it
    with probability zero.  ``MultiStart`` draws its restarts from an LHS over
    the design box, so on this problem it returns its starting point unchanged.

    The test uses a modest sample here; the figure quoted in the README is zero
    feasible points in 12 000 samples, including 4 000 drawn from within 10 % of
    a design known to be feasible.
    """
    rng = np.random.default_rng(0)
    box = Bounds.around(REFINED_DESIGN, relative=0.1)
    points = rng.uniform(box.lower, box.upper, size=(300, len(VARIABLE_NAMES)))
    assert not any(is_feasible(analyse(Design.from_array(row), samples=180)) for row in points)


def test_most_of_the_global_box_is_not_even_analysable() -> None:
    rng = np.random.default_rng(1)
    points = rng.uniform(
        GLOBAL_BOUNDS.lower, GLOBAL_BOUNDS.upper, size=(300, len(VARIABLE_NAMES))
    )
    valid = sum(analyse(Design.from_array(row), samples=180).valid for row in points)
    assert valid / len(points) < 0.25


def test_restarts_land_on_the_equality_manifold() -> None:
    """The projection is what makes a restart usable.

    A perturbation of a feasible design is off both equalities; projecting it
    back puts the restart on the manifold, so the optimizer only has to restore
    the inequalities.  Without that step a restart begins outside a set it
    generally cannot re-enter.
    """
    from exlink.model import equality_constraints
    from exlink.scenarios import project_onto_equalities

    rng = np.random.default_rng(2)
    base = REFINED_DESIGN.to_array()
    moved = Design.from_array(base + rng.normal(0.0, 0.03, base.size) * np.abs(base))
    before = np.max(np.abs(equality_constraints(analyse(moved, samples=360))))
    after = np.max(
        np.abs(
            equality_constraints(
                analyse(project_onto_equalities(moved, samples=360), samples=360)
            )
        )
    )
    assert before > 1.0e-3
    assert after < 1.0e-6


def test_multistart_returns_the_incumbent_when_nothing_beats_it() -> None:
    """A search that finds nothing must not degrade the answer."""
    outcome = multistart(
        "neg_efficiency",
        initial=REFINED_DESIGN,
        bounds=Bounds.around(REFINED_DESIGN, relative=0.05),
        n_start=2,
        spread=0.02,
        max_iter=20,
        samples=360,
    )
    assert outcome.starts == 2
    assert is_feasible(analyse(outcome.design, samples=360))


def test_spread_is_zero_for_a_single_local_optimum() -> None:
    from exlink.scenarios import MultiStartOutcome

    single = MultiStartOutcome(REFINED_DESIGN, -0.28, 3, 1, [-0.28])
    assert single.spread == 0.0
    several = MultiStartOutcome(REFINED_DESIGN, -0.30, 3, 3, [-0.30, -0.27])
    assert several.spread == pytest.approx(0.1)


@pytest.mark.slow
def test_the_single_start_efficiency_optimum_is_only_local() -> None:
    """The answer to "did you use multistart?": no, and it mattered.

    ``GRADIENT_DESIGN`` was produced by one SLSQP run from one starting point
    and reports 30.9 % efficiency.  Restarting from projected perturbations
    finds better than 35 %, so it was a local optimum.

    What the better point is *not* is a better engine: it stands 443 mm tall
    against 320, and it sits exactly on the top-dead-centre gap bound.  The
    single-objective efficiency problem is unbounded in mechanism size -- which
    is the reason this package does not use it -- so a stronger search simply
    exploits that harder.
    """
    outcome = multistart(
        "neg_efficiency",
        initial=GRADIENT_DESIGN,
        bounds=Bounds.around(GRADIENT_DESIGN, relative=0.5),
        n_start=10,
        spread=0.03,
        max_iter=200,
        samples=720,
    )
    incumbent = analyse(GRADIENT_DESIGN, samples=720).metrics.efficiency
    assert outcome.feasible_starts >= 1
    assert -outcome.value > incumbent + 0.03
    assert analyse(outcome.design, samples=720).metrics.height > 400.0
