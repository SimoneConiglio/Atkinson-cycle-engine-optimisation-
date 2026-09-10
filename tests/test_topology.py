"""The spring-connected synthesis, checked against results computed without it.

The load-bearing test is :func:`test_a_single_member_reproduces_a_slider_crank`:
the whole model is an energy minimisation, and if a network whose only member is
a connecting rod does not reproduce the connecting rod's own closed form, nothing
built on it means anything.
"""

from __future__ import annotations

import numpy as np
import pytest

from exlink.synthesis import target_motion
from exlink.topology import (
    COUNT_WEIGHT,
    DISCRETENESS_SCHEDULE,
    GEAR_CATALOGUE,
    HARMONIC_WEIGHT,
    PENALTY_SCHEDULE,
    STIFFNESS_FLOOR,
    STRAIN_WEIGHT,
    TRAVEL_WEIGHT,
    _energy_terms,
    assess,
    best_datum,
    default_bounds,
    engine_ground_structure,
    exlink_in_the_domain,
    gear_radii,
    geared_ground_structure,
    harmonic,
    harmonic_error,
    initial_configuration,
    member_lengths,
    motion_error,
    objective,
    random_start,
    stiffnesses,
    sweep,
    travel_shortfall,
)


def _slider_crank_design(radius: float = 30.0, rod: float = 120.0, axis: float = 0.0):
    """A design vector whose only present member is a connecting rod."""
    structure, layout = engine_ground_structure()
    x = np.zeros(layout.size)
    x[0:2] = [200.0, 200.0]  # the second shaft, parked out of the way
    x[2:4] = [-200.0, 200.0]  # the spare ground pivot, likewise
    x[4:6] = [radius, 0.0]  # input crank: radius, phase
    x[6:8] = [10.0, 0.0]
    x[8:10] = [-150.0, -150.0]  # free nodes, parked
    x[10:12] = [150.0, -150.0]
    x[12] = axis
    x[13] = np.sqrt(rod**2 - (radius - axis) ** 2)
    names = [(structure.names[i], structure.names[j]) for i, j in structure.members]
    rho = np.zeros(structure.n_members)
    rho[names.index(("P1", "S"))] = 1.0
    x[layout.presence_offset :] = rho
    return layout, x


def _geared_slider_crank(structure, layout, ratio: float):
    """A design on the geared domain: one rod from the fast shaft to the piston."""
    x = np.zeros(layout.size)
    x[0:2] = [180.0, 0.0]  # the geared shaft's centre, 180 mm from the input
    x[2:4] = [-150.0, 120.0]
    x[4:6] = [20.0, 0.0]
    x[6:8] = [35.0, 0.0]
    x[8:10] = [-120.0, -120.0]
    x[10:12] = [120.0, -120.0]
    x[12], x[13] = 180.0, 300.0
    x[14] = 0.0
    names = [tuple(structure.names[n] for n in e.nodes) for e in structure.elements]
    rho = np.zeros(structure.n_presences)
    rho[names.index(("P2", "S"))] = 1.0
    pair = next(g for g, gear in enumerate(structure.gears) if gear.ratio == ratio)
    rho[structure.n_members + pair] = 1.0
    x[layout.presence_offset :] = rho
    return x


def test_a_single_member_reproduces_a_slider_crank() -> None:
    """The energy minimisation must agree with the linkage's own closed form.

    With one rod from the crank pin to a piston on the cylinder axis, the piston
    height is ``y_pin + sqrt(L^2 - (x_pin - axis)^2)`` and nothing about springs
    enters it.  Agreement to well under a micron is the evidence that a stiff
    spring is standing in for a rigid link rather than approximating one.
    """
    radius, rod, axis = 30.0, 120.0, 0.0
    layout, x = _slider_crank_design(radius, rod, axis)
    motion = sweep(layout, x, samples=360, penalty=1.0)

    theta = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False)
    pin_x, pin_y = radius * np.cos(theta), radius * np.sin(theta)
    exact = pin_y + np.sqrt(rod**2 - (pin_x - axis) ** 2)

    assert motion.converged
    assert np.max(np.abs(motion.lam - exact)) < 1.0e-3
    assert motion.strain < 1.0e-5
    assert np.ptp(motion.lam) == pytest.approx(2.0 * radius, abs=1.0e-3)


def test_the_energy_derivatives_agree_with_differences() -> None:
    """Gradient and Hessian of the spring energy, against central differences."""
    pos = np.array([[0.0, 0.0], [3.0, 4.0], [10.0, 0.0]])
    ends = (np.array([0, 1]), np.array([1, 2]))
    k, rest = np.array([2.0, 3.0]), np.array([4.0, 6.0])

    def energy(v: np.ndarray) -> float:
        return _energy_terms(v.reshape(3, 2), ends, k, rest)[0]

    def gradient(v: np.ndarray) -> np.ndarray:
        return _energy_terms(v.reshape(3, 2), ends, k, rest)[1]

    _, grad, hess = _energy_terms(pos, ends, k, rest)
    flat, step, basis = pos.reshape(-1), 1.0e-6, np.eye(6)
    fd_grad = np.array(
        [(energy(flat + step * e) - energy(flat - step * e)) / (2.0 * step) for e in basis]
    )
    fd_hess = np.array(
        [(gradient(flat + step * e) - gradient(flat - step * e)) / (2.0 * step) for e in basis]
    )
    assert np.max(np.abs(grad - fd_grad)) < 1.0e-6
    assert np.max(np.abs(hess - fd_hess)) < 1.0e-6
    assert np.max(np.abs(hess - hess.T)) == 0.0


def test_the_objective_is_assembled_from_its_terms() -> None:
    """Motion error, harmonic error, travel shortfall, strain and count.

    Worth pinning because the travel term exists to repair a real failure.  With
    ``COUNT_WEIGHT`` at its original 1.0 mm and no travel term, the cheapest
    design was very nearly no design: a start converged to a fully discrete
    answer whose motion tracked nothing, scoring the target's own standard
    deviation and paying nothing for the members it had dropped.
    """
    _structure, layout = engine_ground_structure()
    x = random_start(layout, np.random.default_rng(0))
    x[layout.presence_offset :] = 0.0
    target = 150.0 + 32.0 * np.cos(2.0 * np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False))

    motion = sweep(layout, x, samples=48, penalty=1.0)
    expected = (
        motion_error(target, motion.lam)
        + HARMONIC_WEIGHT * harmonic_error(target, motion.lam)
        + TRAVEL_WEIGHT * travel_shortfall(target, motion.lam)
        + STRAIN_WEIGHT * motion.strain
    )
    assert objective(layout, x, target, samples=48, penalty=1.0) == pytest.approx(expected)


def test_dropping_every_member_is_not_worth_a_millimetre() -> None:
    """The member count is a tie-break and must never be a reason to build nothing."""
    assert COUNT_WEIGHT * 1.0 < 0.1


def test_absent_members_still_drag_the_piston() -> None:
    """The stiffness floor is not zero, so "no members" is not "no motion".

    A detail with a consequence: an emptied structure does not sit still, it
    follows the floor springs, so a test that asserted a motionless piston would
    be asserting something false.  What makes emptiness unattractive is the
    travel shortfall of the motion it does produce, not the absence of one.
    """
    _structure, layout = engine_ground_structure()
    x = random_start(layout, np.random.default_rng(0))
    x[layout.presence_offset :] = 0.0
    motion = sweep(layout, x, samples=48, penalty=1.0)
    assert motion.converged
    assert motion.strain == pytest.approx(0.0), "an absent member is not charged strain"
    assert np.ptp(motion.lam) > 0.0


def test_travel_shortfall_charges_only_the_shortfall() -> None:
    """Under-travel is charged, exact travel and over-travel are not."""
    target = np.array([0.0, 74.0, 10.0, 60.0])
    assert travel_shortfall(target, np.full(4, 5.0)) == pytest.approx(74.0)
    assert travel_shortfall(target, target) == pytest.approx(0.0)
    assert travel_shortfall(target, 2.0 * target) == pytest.approx(0.0)


def test_motion_error_ignores_an_offset_and_a_datum_but_not_a_shape() -> None:
    """Two nuisances are removed; the thing being measured is not.

    Where the cylinder is bolted and where the cycle is drawn from are both
    choices.  A *different motion* is not, and the distinction is real for this
    target: rotating the datum by ``s`` turns the first harmonic by ``s`` and
    the second by ``2s``, so no single rotation flips both, and an inverted
    motion is charged.  (For a lone sinusoid it would not be -- negation is a
    half-turn of the datum -- which is why this uses the real two-harmonic
    target rather than a sine.)
    """
    target = target_motion(samples=360).lam
    assert motion_error(target, target + 1000.0) == pytest.approx(0.0, abs=1.0e-9)
    assert motion_error(target, np.roll(target, 91)) == pytest.approx(0.0, abs=1.0e-9)
    assert motion_error(target, -target) > 1.0, "an inverted motion is a different one"
    assert motion_error(target, 0.5 * target) > 1.0, "and so is a shallower one"


def test_stiffness_runs_from_the_floor_to_one() -> None:
    """SIMP, with the floor that keeps an abandoned node's Hessian invertible."""
    rho = np.array([0.0, 0.5, 1.0])
    for penalty in PENALTY_SCHEDULE:
        k = stiffnesses(rho, penalty)
        assert k[0] == pytest.approx(STIFFNESS_FLOOR)
        assert k[2] == pytest.approx(1.0)
        assert k[0] < k[1] < k[2]
    assert stiffnesses(rho, 4.0)[1] < stiffnesses(rho, 1.0)[1]


def test_the_ground_structure_excludes_fully_prescribed_members() -> None:
    """A member between two prescribed nodes constrains the design, not the motion."""
    structure, layout = engine_ground_structure()
    movable = {"F1", "F2", "S"}
    for i, j in structure.members:
        pair = {structure.names[i], structure.names[j]}
        assert pair & movable, f"{pair} joins two prescribed nodes"
    assert structure.free_mask().sum() == 5
    assert layout.size == 14 + structure.n_members


def test_member_lengths_come_from_the_initial_pose() -> None:
    """No design is prestressed, because every rest length is measured on it."""
    structure, layout = engine_ground_structure()
    x = random_start(layout, np.random.default_rng(3))
    pos = initial_configuration(layout, x)
    rest = member_lengths(layout, x)
    for m, (i, j) in enumerate(structure.members):
        assert rest[m] == pytest.approx(float(np.hypot(*(pos[i] - pos[j]))))


def test_a_random_start_lies_in_the_box_with_every_member_half_present() -> None:
    """The start is admissible and takes no decision the search should take."""
    _structure, layout = engine_ground_structure()
    lower, upper = default_bounds(layout)
    x = random_start(layout, np.random.default_rng(7))
    assert np.all(x >= lower) and np.all(x <= upper)
    assert np.all(layout.presences(x) == 0.5)


def test_the_schedules_are_paired_and_ramp() -> None:
    """The continuation raises both pressures together, and starts with neither."""
    assert len(DISCRETENESS_SCHEDULE) == len(PENALTY_SCHEDULE)
    assert DISCRETENESS_SCHEDULE[0] == 0.0
    assert list(PENALTY_SCHEDULE) == sorted(PENALTY_SCHEDULE)
    assert list(DISCRETENESS_SCHEDULE) == sorted(DISCRETENESS_SCHEDULE)
    assert STRAIN_WEIGHT > TRAVEL_WEIGHT > COUNT_WEIGHT


def test_a_piston_driven_from_one_shaft_alone_has_no_odd_harmonic() -> None:
    """The structural reason a synthesis can produce a four-stroke Otto engine.

    This is not an observation, it is an identity.  Reach the piston from the
    geared shaft alone and its height is a function of that shaft's angle, which
    is :math:`-2\\theta_1 + \\varphi`; every term is therefore periodic in
    :math:`\\theta_1` with period :math:`\\pi`, so *every odd harmonic vanishes*.
    The two up-and-downs are there and the two halves of the revolution are
    identical, which is exactly a plain Otto motion: four monotone phases, and
    ``STE`` equal to ``STC``.

    It is what the synthesis of §5.6 converged to, and it says where extended
    expansion has to come from -- the piston must be reached from *both* shafts,
    which is what EXlink's trigonal link does.

    The identity is exact; what is measured here is the equilibrium solver's
    precision, so the first harmonic is compared against the second rather than
    against an absolute figure.  It comes out seven orders of magnitude down.
    """
    structure, layout = engine_ground_structure()
    x = np.zeros(layout.size)
    x[0:2] = [-140.0, -215.0]  # the geared shaft's centre
    x[2:4] = [150.0, 40.0]
    x[4:6] = [15.0, 2.1]
    x[6:8] = [36.0, 1.45]  # its crank: radius, phase
    x[8:10] = [-150.0, -150.0]
    x[10:12] = [150.0, -150.0]
    x[12] = -107.0
    names = [(structure.names[i], structure.names[j]) for i, j in structure.members]
    rho = np.zeros(structure.n_members)
    rho[names.index(("P2", "S"))] = 1.0
    x[layout.presence_offset :] = rho
    pos = initial_configuration(layout, x)
    x[13] = pos[structure.index("P2")][1] + 300.0

    target = target_motion(samples=720)
    verdict = assess(layout, x, target.lam, samples=720)

    assert verdict.four_phases, "it is a perfectly good four-stroke motion"
    assert verdict.second_harmonic > 20.0, "the two up-and-downs are there"
    assert verdict.first_harmonic / verdict.second_harmonic < 1.0e-5, (
        "and the asymmetry is absent to the solver's precision"
    )
    assert verdict.asymmetry == pytest.approx(0.0, abs=1.0e-3)
    assert not verdict.is_extended_expansion


def test_a_harmonic_the_target_does_not_have_is_not_charged() -> None:
    """A relative error needs something to be relative to.

    A target with no content at an order still carries a coefficient of order
    1e-15 there rather than exactly zero.  Dividing by that turned a bounded
    objective into 7.8e16 the first time this ran, against a pure
    second-harmonic target; the floor is why it no longer does.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, 96, endpoint=False)
    target = 150.0 + 32.0 * np.cos(2.0 * theta)
    assert abs(harmonic(target, 1)[0]) < 1.0e-12, "the target has no first harmonic"

    assert harmonic_error(target, target) == pytest.approx(0.0, abs=1.0e-9)
    charge = harmonic_error(target, 150.0 + 24.0 * np.cos(2.0 * theta))
    assert 0.0 < charge < float(np.std(target)), "and the second is charged in proportion"


def test_a_fitted_gear_imposes_its_ratio() -> None:
    """The geared shaft must turn at the catalogue ratio, and the other way.

    A pair is carried as slip at the pitch point rather than as an angular
    error, so this also checks the radii: they have to sum to the centre
    distance and stand in the ratio, or the mesh the search chose is not one a
    pair could be cut for.
    """
    structure, layout = geared_ground_structure()
    for ratio in GEAR_CATALOGUE:
        x = _geared_slider_crank(structure, layout, ratio)
        radii = gear_radii(layout, x)
        pair = next(g for g, gear in enumerate(structure.gears) if gear.ratio == ratio)
        assert radii[pair].sum() == pytest.approx(180.0), "the pair fits its centres"
        assert radii[pair][0] / radii[pair][1] == pytest.approx(ratio)

        motion = sweep(layout, x, samples=720, penalty=PENALTY_SCHEDULE[-1])
        assert motion.converged
        assert motion.strain < 1.0e-4, "a fitted pair rolls without slipping"
        assert motion.gear_ratios == (ratio,)
        slope = np.sign(np.diff(np.concatenate([motion.lam, motion.lam[:1]])))
        reversals = int((slope[:-1] != slope[1:]).sum() + (slope[-1] != slope[0]))
        assert reversals // 2 == int(ratio), "one up-and-down per turn of the shaft"


def test_only_a_whole_ratio_closes_the_cycle() -> None:
    """The catalogue is integers, and that is physics rather than convenience.

    A geared shaft turning ``r`` times per input revolution is back where it
    started only when ``r`` is whole.  At three-to-two the mechanism's period is
    *two* input revolutions, so successive cycles of the engine differ and there
    is no four-stroke to speak of.
    """
    assert all(float(r).is_integer() for r in GEAR_CATALOGUE)


def test_a_body_is_one_presence_driving_three_rigid_sides() -> None:
    """The change §5.6's stall asked for: a triangle switched on by one variable.

    Assembled from three separate bars a triangle is unreachable by descent --
    each bar alone does nothing, so the objective is flat in its geometry until
    all three happen to switch on together.  One presence removes that barrier,
    and the three sides have to stay rigid for it to be a body at all.
    """
    structure, _layout = geared_ground_structure()
    bodies = [e for e in structure.elements if e.kind == "body"]
    assert bodies, "the geared domain offers three-cornered links"
    assert all(len(e.springs) == 3 for e in bodies)
    assert all(len(e.springs) == 1 for e in structure.elements if e.kind == "bar")

    owner = structure.spring_owner()
    assert owner.size == len(structure.springs())
    for m, element in enumerate(structure.elements):
        assert int((owner == m).sum()) == len(element.springs)


def test_a_body_stays_rigid_through_the_revolution() -> None:
    """A switched-on triangle holds all three of its sides, not just two."""
    structure, layout = geared_ground_structure()
    x = _geared_slider_crank(structure, layout, 2.0)
    names = [tuple(structure.names[n] for n in e.nodes) for e in structure.elements]
    rho = layout.presences(x).copy()
    rho[names.index(("P2", "F1", "S"))] = 1.0
    x[layout.presence_offset :] = rho

    rest = member_lengths(layout, x)
    body = next(
        m
        for m, e in enumerate(structure.elements)
        if e.nodes == tuple(structure.index(n) for n in ("P2", "F1", "S"))
    )
    sides = np.flatnonzero(structure.spring_owner() == body)
    assert sides.size == 3
    assert np.all(rest[sides] > 0.0)

    motion = sweep(layout, x, samples=360, penalty=PENALTY_SCHEDULE[-1])
    assert motion.converged
    assert motion.strain < 1.0e-4, "every side of the body holds its length"


def test_the_geared_domain_contains_the_older_one() -> None:
    """§5.7 adds freedom; it does not take any away."""
    old, old_layout = engine_ground_structure()
    new, new_layout = geared_ground_structure()
    assert set(old.members) <= set(new.members)
    assert new.n_members > old.n_members
    assert new.gears and not old.gears
    assert old.free_shafts == () and new.free_shafts == (1,)
    assert new_layout.size > old_layout.size


def test_the_domain_can_express_the_studied_mechanism() -> None:
    """The check that says whether a negative result is the domain or the search.

    §5.6 and §5.7 searched and found no extended-expansion mechanism.  That
    means nothing until the domain is shown to contain one, and it does: EXlink
    is three of its candidates and a gear pair -- the swing rod, the trigonal
    link as a *body*, the piston rod, and a two-to-one mesh -- and the spring
    model reproduces its analytic kinematics to a couple of microns.

    Every corner is a design variable, which is what makes it expressible: the
    two free nodes are the trigonal link's own corners, and the body's three
    side lengths follow from where they are put.
    """
    from exlink.kinematics import solve
    from exlink.reference import RELIABLE_DESIGN

    _structure, layout, x = exlink_in_the_domain(samples=720)
    motion = sweep(layout, x, samples=720, penalty=PENALTY_SCHEDULE[-1])
    exact = solve(RELIABLE_DESIGN, samples=720).lam

    difference = (motion.lam - motion.lam.mean()) - (exact - exact.mean())
    assert np.max(np.abs(difference)) < 0.05, "the spring model is the same mechanism"
    assert motion.converged
    assert motion.strain < 1.0e-3, "it runs as a linkage, not as a deforming structure"
    assert motion.output_slack == pytest.approx(0.0), "and the input determines it"
    assert motion.gear_ratios == (2.0,)

    verdict = assess(layout, x, target_motion(samples=720).lam, samples=720)
    assert verdict.is_extended_expansion
    assert verdict.expansion_stroke == pytest.approx(74.0, abs=0.2)
    assert verdict.compression_ratio == pytest.approx(16.0, abs=0.2)
    assert len(verdict.present) == 3
    assert any(len(nodes) == 3 for nodes in verdict.present), "one of them is a body"


def test_the_objective_does_not_charge_for_the_crank_datum() -> None:
    """Where theta_1 = 0 sits is a choice, and it used to cost 10 mm.

    Turn a design about its input axis and every crank phase moves with it: the
    same engine, its cycle starting somewhere else.  Charging for that inverted
    the ranking between the studied mechanism and a degenerate one -- EXlink
    sits 170 degrees from the target's datum and scored 13.97 mm as posed
    against 4.05 mm aligned, which put it behind an Otto engine that had tuned
    its phase to 7.19 mm.
    """
    target = target_motion(samples=360).lam
    for shift in (0, 37, 180, 259):
        rolled = np.roll(target - target.mean(), shift) + 100.0
        assert motion_error(target, rolled) == pytest.approx(0.0, abs=1.0e-9)
        assert harmonic_error(target, rolled) == pytest.approx(0.0, abs=1.0e-6)
        assert best_datum(target, rolled) == (target.size - shift) % target.size


def test_the_studied_mechanism_now_scores_best() -> None:
    """With the datum free, the objective ranks EXlink ahead of an Otto engine.

    The point of the repair.  A synthesis that scores the answer it is looking
    for *behind* a degenerate one cannot find it however long it searches, and
    under the objective as first posed that is exactly what would have happened.
    """
    target = target_motion(samples=720)
    _structure, layout, x = exlink_in_the_domain(samples=720)
    ours = objective(layout, x, target.lam, samples=720, penalty=PENALTY_SCHEDULE[-1])

    otto = np.roll(_two_up_and_downs(target.lam), 0)
    theirs = (
        motion_error(target.lam, otto)
        + HARMONIC_WEIGHT * harmonic_error(target.lam, otto)
        + TRAVEL_WEIGHT * travel_shortfall(target.lam, otto)
    )
    assert ours < theirs, "the mechanism with the asymmetry wins"
    assert ours < 10.0


def _two_up_and_downs(target: np.ndarray) -> np.ndarray:
    """A pure second harmonic of the target's own amplitude: an Otto motion.

    Two equal up-and-downs, no asymmetry at all -- what a piston reached from
    the geared shaft alone produces, and what both searches converged to.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, target.size, endpoint=False)
    amplitude = float(np.hypot(*harmonic(target, 2)))
    return amplitude * np.cos(2.0 * theta)
