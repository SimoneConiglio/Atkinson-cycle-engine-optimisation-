"""Synthesising the mechanism itself, rather than the dimensions of a given one.

What this is for
----------------
§5.4 measures the EXlink topology against a slider-crank and finds it worth
17.6 %.  Two topologies establish a contrast and not a trend, and §6.3 records
the missing third -- with the warning that choosing one by hand moves the
arbitrariness rather than removing it.  This module removes it instead, by
asking the optimizer for the topology.

The method is the second family of §2.2: the spring-connected model of
:cite:t:`kim2007spring`, in the reduced *link* form of
:cite:t:`tran2024slm`.  A design domain is populated with candidate members,
every candidate is a spring whose stiffness is a design variable, and the
stiffnesses are penalised towards their extremes.  A spring driven stiff is a
rigid binary link; a spring driven soft is no link at all.  One continuous
parameterisation therefore covers every linkage the domain can hold, and the
topology falls out of a gradient solve rather than being assumed.

The synthesis problem
---------------------
The input is the **half-speed shaft**, turning once per cycle, and the output is
the piston height.  Over one input revolution the target motion has two equal
maxima and two unequal minima -- :func:`exlink.synthesis.target_motion` builds
the one that realises ``STE = 74`` mm and ``epsilon = 16`` exactly -- so the
synthesis is asked for a linkage that turns one circular input into two
unequal up-and-downs of a point on the cylinder axis.

That period is the whole difficulty, and it is worth being explicit about why.
A linkage driven by a shaft is periodic in that shaft's angle, so a mechanism
driven from the *crankshaft* alone cannot produce a motion that differs between
its two revolutions: extended expansion needs something in the mechanism that
turns once per cycle.  Every four-stroke already has one, and the candidate set
below offers the search both ways of using it -- a second grounded shaft geared
to the input, which is what EXlink does, or a longer chain of bars off the input
shaft alone, which is what a six-bar function generator would do.  Which it
takes is an outcome rather than an assumption.

The model
---------
Nodes carry planar coordinates.  Some are grounded, two are *driven* -- pins on
cranks whose angles are prescribed multiples of the input angle -- one is the
piston, whose abscissa is the cylinder axis and whose ordinate is the output,
and the rest are free.  Every candidate member :math:`m` joining nodes
:math:`i` and :math:`j` contributes

.. math:: E_m = \\tfrac12 k_m \\bigl(\\lVert x_i - x_j \\rVert - L_m\\bigr)^2,
          \\qquad k_m = k_{\\min} + (1 - k_{\\min})\\,\\rho_m^{\\,p}

with :math:`\\rho_m \\in [0, 1]` its presence and :math:`L_m` its length.  The
configuration at an input angle is the minimiser of :math:`\\sum_m E_m` over the
free coordinates, found by a damped Newton iteration warm-started from the
previous angle -- the incremental analysis the published method uses, and the
reason the branch a mechanism runs on is followed rather than re-chosen at every
angle.

The penalty exponent :math:`p` is raised on a schedule (:data:`PENALTY_SCHEDULE`).
At :math:`p = 1` the problem is nearly convex in the stiffnesses and the search
moves freely; raising it makes intermediate stiffnesses uneconomic, so the
members separate into present and absent.

Why the objective can afford to be only the motion
--------------------------------------------------
Nothing here prices a part, and that is a real restriction rather than an
oversight: §2.6 records that the topology-synthesis literature targets a
kinematic quantity and therefore sizes nothing.  What this module produces is a
*candidate*, and a candidate is exactly what the rest of the package is for --
the range chain of §3.3 scores it afterwards.  Synthesising against range
directly would need the sizing fixed point inside every equilibrium solve at
every angle, which is a study of its own (§6.3).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]

STIFFNESS_FLOOR = 1.0e-6
"""Stiffness of a fully absent member, relative to a fully present one.

Not zero, for the reason SIMP keeps a floor: a node reachable only through
absent members would otherwise leave the Hessian singular, and the damped
Newton would have nothing to say about where it goes.
"""

PENALTY_SCHEDULE: tuple[float, ...] = (1.0, 2.0, 3.0, 4.0)
"""Exponents :math:`p` the continuation walks through."""

GROUND = -1
"""Node kind: coordinates fixed by the design vector."""
DRIVEN = -2
"""Node kind: coordinates prescribed by the input angle."""
SLIDER = -3
"""Node kind: abscissa fixed, ordinate free.  The piston."""
FREE = -4
"""Node kind: both coordinates free."""


@dataclass(frozen=True)
class GroundStructure:
    """The candidate nodes and members a linkage may be assembled from.

    Instances are pure topology -- which nodes exist, what kind each is, and
    which pairs a member may join.  Everything dimensional lives in the design
    vector, so one ground structure serves every design drawn on it.
    """

    names: tuple[str, ...]
    """Node names, in coordinate order."""

    kinds: tuple[int, ...]
    """One of :data:`GROUND`, :data:`DRIVEN`, :data:`SLIDER`, :data:`FREE`."""

    members: tuple[tuple[int, int], ...]
    """Candidate members, as pairs of node indices."""

    speed_ratios: dict[int, float] = field(default_factory=dict)
    """For each driven node, its crank speed as a multiple of the input's.

    ``1`` is a pin on the input shaft itself.  ``-2`` is a pin on a shaft geared
    to it two to one and turning the other way, which is the relation
    :mod:`exlink.kinematics` gives the EXlink.
    """

    @property
    def n_nodes(self) -> int:
        """Nodes in the structure."""
        return len(self.names)

    @property
    def n_members(self) -> int:
        """Candidate members."""
        return len(self.members)

    def free_mask(self) -> NDArray[np.bool_]:
        """Which of the ``2 n_nodes`` coordinates the equilibrium solves for."""
        mask = np.zeros(2 * self.n_nodes, dtype=bool)
        for i, kind in enumerate(self.kinds):
            if kind == FREE:
                mask[2 * i : 2 * i + 2] = True
            elif kind == SLIDER:
                mask[2 * i + 1] = True
        return mask

    def ends(self) -> tuple[IntArray, IntArray]:
        """The members' first and second node indices, as arrays."""
        pairs = np.asarray(self.members, dtype=np.int64).reshape(-1, 2)
        return pairs[:, 0], pairs[:, 1]

    def index(self, name: str) -> int:
        """Index of a node by name."""
        return self.names.index(name)


@dataclass(frozen=True)
class Layout:
    """Where each quantity sits in the design vector.

    The vector holds the *initial* configuration and the member presences, and
    nothing else.  Member lengths are deliberately absent: each is taken as the
    distance between its endpoints in the initial pose, so every design is
    assembled without prestress and an entire class of structures that could
    never be built is never visited.
    """

    structure: GroundStructure
    ground_slots: tuple[tuple[int, int], ...]
    """``(node, offset)`` for every grounded node whose position is free."""
    driven_slots: tuple[tuple[int, int], ...]
    """``(node, offset)`` for every driven pin: radius then phase."""
    free_slots: tuple[tuple[int, int], ...]
    """``(node, offset)`` for every node whose initial position is chosen."""
    slider_slot: tuple[int, int]
    """``(node, offset)`` of the piston: cylinder abscissa then initial height."""
    presence_offset: int
    """First index of the member presences."""
    centre_of: dict[int, int]
    """Grounded centre each driven pin turns about."""

    @property
    def size(self) -> int:
        """Length of the design vector."""
        return self.presence_offset + self.structure.n_members

    def presences(self, x: FloatArray) -> FloatArray:
        """The presence variables of a design vector."""
        return np.asarray(x[self.presence_offset :], dtype=float)


def _fixed_positions(layout: Layout, x: FloatArray) -> FloatArray:
    """Node coordinates that do not depend on the input angle."""
    pos = np.zeros((layout.structure.n_nodes, 2), dtype=float)
    for node, offset in layout.ground_slots:
        pos[node] = x[offset : offset + 2]
    return pos


def initial_configuration(layout: Layout, x: FloatArray) -> FloatArray:
    """The full node layout at zero input angle, free nodes included."""
    pos = _fixed_positions(layout, x)
    for node, offset in layout.driven_slots:
        radius, phase = float(x[offset]), float(x[offset + 1])
        ratio = layout.structure.speed_ratios[node]
        angle = ratio * 0.0 + phase
        pos[node] = pos[layout.centre_of[node]] + radius * np.array(
            [np.cos(angle), np.sin(angle)]
        )
    for node, offset in layout.free_slots:
        pos[node] = x[offset : offset + 2]
    node, offset = layout.slider_slot
    pos[node] = np.array([x[offset], x[offset + 1]])
    return pos


def member_lengths(layout: Layout, x: FloatArray) -> FloatArray:
    """Rest length of every candidate member, from the initial pose."""
    pos = initial_configuration(layout, x)
    out = np.empty(layout.structure.n_members, dtype=float)
    for m, (i, j) in enumerate(layout.structure.members):
        out[m] = float(np.hypot(*(pos[i] - pos[j])))
    return out


def stiffnesses(presence: FloatArray, penalty: float) -> FloatArray:
    """SIMP stiffness of each member: a floor plus its penalised presence."""
    rho = np.clip(np.asarray(presence, dtype=float), 0.0, 1.0)
    return STIFFNESS_FLOOR + (1.0 - STIFFNESS_FLOOR) * rho**penalty


def _energy_terms(
    pos: FloatArray,
    ends: tuple[IntArray, IntArray],
    k: FloatArray,
    rest: FloatArray,
) -> tuple[float, FloatArray, FloatArray]:
    """Total spring energy and its first two derivatives in the coordinates.

    The blocks are the textbook ones for an axial spring, assembled over the
    full coordinate set with the caller restricting to the free ones.  It is
    written over all members at once because it is the innermost thing in the
    method -- every Newton step of every angle of every finite difference -- and
    a Python loop over members costs more than the arithmetic does.

    Args:
        pos: ``(n_nodes, 2)`` coordinates.
        ends: The members' first and second node indices.
        k: Member stiffnesses.
        rest: Member rest lengths.

    Returns:
        Energy, a ``2 n_nodes`` gradient, and the ``2n x 2n`` Hessian.
    """
    i, j = ends
    n = pos.shape[0]
    d = pos[i] - pos[j]
    r = np.sqrt(np.einsum("mi,mi->m", d, d))
    r = np.maximum(r, 1.0e-12)
    stretch = r - rest
    total = float(0.5 * np.dot(k, stretch**2))

    unit = d / r[:, None]
    f = (k * stretch)[:, None] * unit
    grad = np.zeros((n, 2), dtype=float)
    np.add.at(grad, i, f)
    np.add.at(grad, j, -f)

    outer = unit[:, :, None] * unit[:, None, :]
    eye = np.eye(2)[None, :, :]
    block = k[:, None, None] * (outer + (stretch / r)[:, None, None] * (eye - outer))
    hess = np.zeros((2 * n, 2 * n), dtype=float)
    rows = np.stack([2 * i, 2 * i + 1], axis=1)
    cols = np.stack([2 * j, 2 * j + 1], axis=1)
    for first, second, sign in (
        (rows, rows, 1.0),
        (cols, cols, 1.0),
        (rows, cols, -1.0),
        (cols, rows, -1.0),
    ):
        np.add.at(hess, (first[:, :, None], second[:, None, :]), sign * block)
    return total, grad.reshape(-1), hess


def driven_positions(layout: Layout, x: FloatArray, theta: float) -> FloatArray:
    """Where each driven pin sits at an input angle."""
    pos = _fixed_positions(layout, x)
    out = np.zeros((layout.structure.n_nodes, 2), dtype=float)
    for node, offset in layout.driven_slots:
        radius, phase = float(x[offset]), float(x[offset + 1])
        angle = layout.structure.speed_ratios[node] * theta + phase
        out[node] = pos[layout.centre_of[node]] + radius * np.array(
            [np.cos(angle), np.sin(angle)]
        )
    return out


NEWTON_STEPS = 30
"""Newton iterations allowed per input angle."""
NEWTON_TOLERANCE = 1.0e-9
"""Convergence test on the free-coordinate gradient, relative to its scale.

The gradient of a spring energy carries the stiffness and a length, and both
vary by orders of magnitude across a penalisation schedule, so an absolute
threshold would be met instantly at one end of it and never at the other.
"""
RIDGE = 1.0e-9
"""Fraction of the mean diagonal added to the Hessian before solving.

The Hessian of a spring network is singular exactly when some node is held only
by members the search has switched off -- which, in a topology optimization, is
most of the time.  The ridge leaves such a node where the warm start put it
instead of sending it to infinity, and is small enough not to disturb a node
that is properly held.
"""
BACKTRACK_STEPS = 20
"""Halvings allowed before a Newton step is abandoned."""


def _equilibrate(
    layout: Layout,
    x: FloatArray,
    theta: float,
    guess: FloatArray,
    k: FloatArray,
    rest: FloatArray,
) -> tuple[FloatArray, bool]:
    """Minimise the spring energy over the free coordinates at one angle.

    Newton with a ridge and a backtracking line search.  The ridge handles the
    singular Hessian of a partly switched-off network; the line search handles
    the fact that a spring energy is not convex in the coordinates once a member
    is in compression, which is where a bare Newton step overshoots.

    Args:
        layout: The design vector's layout.
        x: Design vector.
        theta: Input angle [rad].
        guess: Free coordinates to start from, normally the previous angle's.
        k: Member stiffnesses.
        rest: Member rest lengths.

    Returns:
        The free coordinates and whether the gradient test was met.
    """
    structure = layout.structure
    mask = structure.free_mask()
    ends = structure.ends()
    base = _fixed_positions(layout, x) + driven_positions(layout, x, theta)
    slider_node, slider_offset = layout.slider_slot
    base[slider_node, 0] = float(x[slider_offset])

    scale = max(1.0, float(np.max(k)) * float(np.max(rest, initial=1.0)))
    tol = NEWTON_TOLERANCE * scale
    n_free = int(mask.sum())
    eye = np.eye(n_free)

    def energy_at(free: FloatArray) -> tuple[float, FloatArray, FloatArray]:
        pos = base.copy()
        pos.reshape(-1)[mask] = free
        return _energy_terms(pos, ends, k, rest)

    u = np.array(guess, dtype=float)
    for _ in range(NEWTON_STEPS):
        energy, grad, hess = energy_at(u)
        g = grad[mask]
        if float(np.max(np.abs(g))) < tol:
            return u, True
        h = hess[np.ix_(mask, mask)]
        ridge = RIDGE * max(float(np.mean(np.abs(np.diag(h)))), 1.0e-30)
        try:
            step = np.linalg.solve(h + ridge * eye, -g)
        except np.linalg.LinAlgError:
            return u, False
        alpha = 1.0
        for _back in range(BACKTRACK_STEPS):
            trial = u + alpha * step
            if energy_at(trial)[0] <= energy:
                u = trial
                break
            alpha *= 0.5
        else:
            return u, False
    return u, float(np.max(np.abs(energy_at(u)[1][mask]))) < tol


@dataclass(frozen=True)
class Motion:
    """What a candidate structure does over one input revolution."""

    lam: FloatArray
    """Piston height at each input angle [mm]."""
    converged: bool
    """Whether every angle reached equilibrium."""
    strain: float
    """Largest member strain over the revolution, as a fraction of rest length.

    The honest measure of whether the springs are behaving as rigid links.  A
    linkage the penalisation has resolved runs at a strain of order the
    stiffness floor; anything larger means the answer is a deforming structure
    rather than a mechanism, and is not a linkage at all.
    """


def sweep(
    layout: Layout,
    x: FloatArray,
    samples: int = 180,
    penalty: float = 1.0,
) -> Motion:
    """Run a candidate through one revolution of the input shaft.

    The angles are walked in order from zero and each equilibrium is warm-started
    from the last, which is what makes the result a *mechanism's* motion: the
    configuration follows one branch of the loop-closure solution rather than
    being re-chosen at every angle.

    Args:
        layout: The design vector's layout.
        x: Design vector.
        samples: Input angles over ``[0, 2 pi)``.
        penalty: SIMP exponent applied to the presences.

    Returns:
        The piston motion and two diagnostics.
    """
    structure = layout.structure
    mask = structure.free_mask()
    first, second = structure.ends()
    presence = layout.presences(x)
    k = stiffnesses(presence, penalty)
    rest = member_lengths(layout, x)
    safe_rest = np.where(rest > 0.0, rest, 1.0)

    start = initial_configuration(layout, x)
    u = start.reshape(-1)[mask].copy()

    angles = np.linspace(0.0, 2.0 * np.pi, int(samples), endpoint=False)
    lam = np.empty(angles.size, dtype=float)
    slider_node, slider_offset = layout.slider_slot
    ok = True
    worst = 0.0
    for step, theta in enumerate(angles):
        u, converged = _equilibrate(layout, x, float(theta), u, k, rest)
        ok = ok and converged
        pos = _fixed_positions(layout, x) + driven_positions(layout, x, float(theta))
        pos[slider_node, 0] = float(x[slider_offset])
        pos.reshape(-1)[mask] = u
        lam[step] = pos[slider_node, 1]
        gap = pos[first] - pos[second]
        length = np.sqrt(np.einsum("mi,mi->m", gap, gap))
        strain = np.abs(length - rest) / safe_rest
        worst = max(worst, float(np.max(strain * presence, initial=0.0)))
    return Motion(lam=lam, converged=ok, strain=worst)


def engine_ground_structure() -> tuple[GroundStructure, Layout]:
    """The candidate set this study searches.

    Eight nodes and twelve candidate members.  The input shaft is grounded at
    the origin -- a gauge choice, since translating the whole mechanism changes
    nothing -- and carries the pin ``P1``.  A second grounded shaft ``B`` carries
    ``P2`` at twice the input speed and the opposite sense, which is the relation
    a 2:1 gear pair imposes and the one EXlink uses; the search is free to leave
    it out by switching off every member that touches it, in which case the
    answer is a linkage driven from the input shaft alone.  ``C`` is a plain
    grounded pivot, ``F1`` and ``F2`` are free nodes, and ``S`` is the piston,
    held on the cylinder axis.

    Members between two prescribed nodes are excluded: their length is a
    function of the input angle alone, so such a member is a constraint on the
    design rather than a degree of freedom for the mechanism.

    Returns:
        The structure and the layout of the design vector drawn on it.
    """
    names = ("A", "P1", "B", "P2", "C", "F1", "F2", "S")
    kinds = (GROUND, DRIVEN, GROUND, DRIVEN, GROUND, FREE, FREE, SLIDER)
    idx = {name: i for i, name in enumerate(names)}
    movable = ("P1", "P2", "C", "F1", "F2", "S")
    members = tuple(
        (idx[a], idx[b])
        for i, a in enumerate(movable)
        for b in movable[i + 1 :]
        if kinds[idx[a]] in (FREE, SLIDER) or kinds[idx[b]] in (FREE, SLIDER)
    )
    structure = GroundStructure(
        names=names,
        kinds=kinds,
        members=members,
        speed_ratios={idx["P1"]: 1.0, idx["P2"]: -2.0},
    )
    layout = Layout(
        structure=structure,
        ground_slots=((idx["B"], 0), (idx["C"], 2)),
        driven_slots=((idx["P1"], 4), (idx["P2"], 6)),
        free_slots=((idx["F1"], 8), (idx["F2"], 10)),
        slider_slot=(idx["S"], 12),
        presence_offset=14,
        centre_of={idx["P1"]: idx["A"], idx["P2"]: idx["B"]},
    )
    return structure, layout


def default_bounds(layout: Layout) -> tuple[FloatArray, FloatArray]:
    """Box the search runs in, in millimetres and radians.

    Deliberately loose.  The point of a topology search is that the answer is
    not known in advance, and a box tight enough to be informative would be
    assuming it.
    """
    n = layout.size
    lower = np.empty(n, dtype=float)
    upper = np.empty(n, dtype=float)
    span = 260.0
    for _node, offset in layout.ground_slots:
        lower[offset : offset + 2] = -span
        upper[offset : offset + 2] = span
    for _node, offset in layout.driven_slots:
        lower[offset], upper[offset] = 5.0, 120.0
        lower[offset + 1], upper[offset + 1] = -np.pi, np.pi
    for _node, offset in layout.free_slots:
        lower[offset : offset + 2] = -span
        upper[offset : offset + 2] = span
    _node, offset = layout.slider_slot
    lower[offset], upper[offset] = -span, span
    lower[offset + 1], upper[offset + 1] = -span, span
    lower[layout.presence_offset :] = 0.0
    upper[layout.presence_offset :] = 1.0
    return lower, upper


COUNT_WEIGHT = 0.05
"""Millimetres of motion error a whole extra member is worth.

A tie-break and nothing more.  It began at 1.0 mm, which was a mistake worth
recording: a design with no members at all leaves the piston free, the piston
then stays where the warm start put it, and a *constant* motion scores the
target's own standard deviation -- 23.7 mm here.  At 1.0 mm the empty structure
was therefore cheaper than any partly-connected one, and one start duly
converged to it, fully discrete and completely useless.  The travel floor below
is the real repair; this weight is now small enough not to matter either way.
"""

TRAVEL_WEIGHT = 2.0
"""Millimetres charged per millimetre of stroke the piston fails to deliver.

Without it "build nothing" is a local minimum of the motion error, for the
reason above.  With it a static piston is charged twice the required stroke,
which no moving mechanism can be worse than, so the search is obliged to
connect the piston to something before it can start improving the shape.
"""

STRAIN_WEIGHT = 800.0
"""Millimetres of motion error a unit of member strain is worth.

This is the term that decides whether the answer is a mechanism, and it began
an order of magnitude too small.  At 50 a start reached a motion error of
2.0 mm while carrying 3 % strain, and paid 1.5 mm for it -- so the search was
being *paid* to leave members half present, because a half-stiff member lets an
overconstrained network move by giving rather than by articulating.  The gray
design fitted the target well and rounded to nothing.  At 800 the same 3 % costs
25 mm and softness stops being a way to buy motion.
"""

DISCRETENESS_SCHEDULE: tuple[float, ...] = (0.0, 2.0, 6.0, 18.0)
r"""Millimetres charged for a fully undecided member, one per penalty rung.

The SIMP exponent alone did not separate the members here.  Ramping an explicit
:math:`4\rho(1-\rho)` charge alongside it is the standard repair, and the ramp
matters: charged from the start it decides the topology before the geometry is
any good, and never charged at all it leaves an answer that cannot be rounded.
"""


@dataclass(frozen=True)
class Candidate:
    """A linkage the synthesis produced, and how well it does."""

    x: FloatArray
    """Design vector."""
    rms: float
    """Motion error against the target, mean removed [mm]."""
    strain: float
    """Largest strain in a present member (:attr:`Motion.strain`)."""
    discreteness: float
    """Mean of :math:`4\\rho(1-\\rho)`: zero when every member has resolved."""
    present: tuple[tuple[str, str], ...]
    """Members that survived, by node name."""
    motion: FloatArray
    """The piston motion it realises [mm]."""
    converged: bool
    """Whether every equilibrium solve met its gradient test."""

    @property
    def uses_second_shaft(self) -> bool:
        """Whether any surviving member touches the geared shaft's pin."""
        return any("P2" in pair for pair in self.present)


def _centred(values: FloatArray) -> FloatArray:
    """Values with their mean removed, so an offset is not an error."""
    arr = np.asarray(values, dtype=float)
    return arr - float(arr.mean())


def motion_error(target: FloatArray, lam: FloatArray) -> float:
    """Root-mean-square distance between two motions, offset removed.

    The target's mean height is arbitrary -- it is set by where the cylinder is
    bolted, not by the mechanism -- so comparing raw heights would charge a
    design for a quantity the objective does not care about.
    """
    a, b = _centred(target), _centred(lam)
    if a.size != b.size:
        grid = np.linspace(0.0, 1.0, a.size, endpoint=False)
        b = np.interp(grid, np.linspace(0.0, 1.0, b.size, endpoint=False), b, period=1.0)
    return float(np.sqrt(np.mean((a - b) ** 2)))


def travel_shortfall(target: FloatArray, lam: FloatArray) -> float:
    """Stroke the piston fails to deliver, in millimetres, never negative.

    Overshoot is not charged here: a design that moves too far is wrong in a way
    the motion error already sees, and charging it twice would bias the search
    towards under-travel -- which is the failure this term exists to prevent.
    """
    required = float(np.ptp(np.asarray(target, dtype=float)))
    return max(0.0, required - float(np.ptp(np.asarray(lam, dtype=float))))


def objective(
    layout: Layout,
    x: FloatArray,
    target: FloatArray,
    samples: int,
    penalty: float,
    discreteness: float = 0.0,
) -> float:
    """What the search minimises: the motion error, and four deterrents.

    Args:
        layout: The design vector's layout.
        x: Design vector.
        target: Target piston motion.
        samples: Input angles per sweep.
        penalty: SIMP exponent on the presences.
        discreteness: Millimetres charged for a fully undecided member.

    Returns:
        The objective, in millimetres of equivalent motion error.
    """
    motion = sweep(layout, x, samples=samples, penalty=penalty)
    if not np.all(np.isfinite(motion.lam)):
        return 1.0e6
    rho = layout.presences(x)
    return (
        motion_error(target, motion.lam)
        + TRAVEL_WEIGHT * travel_shortfall(target, motion.lam)
        + STRAIN_WEIGHT * motion.strain
        + discreteness * float(np.mean(4.0 * rho * (1.0 - rho)))
        + COUNT_WEIGHT * float(np.mean(rho))
    )


def random_start(layout: Layout, rng: np.random.Generator) -> FloatArray:
    """A design drawn from the box, with every member half present.

    Starting the presences at one half rather than at random is deliberate: at
    :math:`p = 1` the stiffnesses enter almost linearly, so a symmetric start
    lets the first rung of the continuation decide which members matter instead
    of inheriting a decision from the draw.
    """
    lower, upper = default_bounds(layout)
    x = lower + rng.random(layout.size) * (upper - lower)
    x[layout.presence_offset :] = 0.5
    return x


def _minimise(
    layout: Layout,
    x0: FloatArray,
    target: FloatArray,
    samples: int,
    penalty: float,
    iterations: int,
    frozen: NDArray[np.bool_] | None = None,
    discreteness: float = 0.0,
) -> FloatArray:
    """One rung: bound-constrained descent on :func:`objective`.

    L-BFGS-B rather than SLSQP because the only constraints are the bounds, and
    because a quasi-Newton method tolerates the finite-difference gradient this
    objective is stuck with -- an equilibrium sweep is differentiable in
    principle but the branch it follows is chosen by a warm start, so an
    analytic derivative would have to differentiate the warm start too.
    """
    from scipy.optimize import minimize

    lower, upper = default_bounds(layout)
    if frozen is not None:
        lower = np.where(frozen, x0, lower)
        upper = np.where(frozen, x0, upper)
    result = minimize(
        lambda v: objective(layout, v, target, samples, penalty, discreteness),
        np.clip(x0, lower, upper),
        method="L-BFGS-B",
        bounds=list(zip(lower, upper, strict=True)),
        options={"maxiter": int(iterations), "eps": 1.0e-4},
    )
    return np.asarray(result.x, dtype=float)


PRESENCE_THRESHOLD = 0.5
"""Presence above which a member is kept when the answer is made discrete."""


def _describe(layout: Layout, x: FloatArray, target: FloatArray, samples: int) -> Candidate:
    """Score a design at full penalty and name the members it kept."""
    motion = sweep(layout, x, samples=samples, penalty=PENALTY_SCHEDULE[-1])
    rho = layout.presences(x)
    names = layout.structure.names
    present = tuple(
        (names[i], names[j])
        for m, (i, j) in enumerate(layout.structure.members)
        if rho[m] > PRESENCE_THRESHOLD
    )
    return Candidate(
        x=np.array(x, dtype=float),
        rms=motion_error(target, motion.lam),
        strain=motion.strain,
        discreteness=float(np.mean(4.0 * rho * (1.0 - rho))),
        present=present,
        motion=motion.lam,
        converged=motion.converged,
    )


def synthesise_one(
    layout: Layout,
    x0: FloatArray,
    target: FloatArray,
    samples: int = 90,
    iterations: int = 120,
    schedule: tuple[float, ...] = PENALTY_SCHEDULE,
) -> Candidate:
    """Walk one start up the penalisation schedule and make the answer discrete.

    The last step is the one that turns a topology optimization into a
    mechanism.  Rounding the presences to zero and one gives a linkage; the
    geometry is then re-optimised with that topology frozen, so what is reported
    is a real linkage's best fit rather than a blurred one's.

    Args:
        layout: The design vector's layout.
        x0: Starting design.
        target: Target motion, on any number of samples.
        samples: Input angles per sweep during the search.
        iterations: L-BFGS-B iterations per rung.
        schedule: Penalisation exponents.

    Returns:
        The discrete linkage and its score.
    """
    x = np.array(x0, dtype=float)
    for penalty, sharpness in zip(schedule, DISCRETENESS_SCHEDULE, strict=False):
        x = _minimise(layout, x, target, samples, penalty, iterations, discreteness=sharpness)
    rho = layout.presences(x)
    x[layout.presence_offset :] = (rho > PRESENCE_THRESHOLD).astype(float)
    frozen = np.zeros(layout.size, dtype=bool)
    frozen[layout.presence_offset :] = True
    x = _minimise(layout, x, target, samples, schedule[-1], iterations, frozen=frozen)
    return _describe(layout, x, target, samples=len(target))


def harmonic(lam: FloatArray, order: int) -> tuple[float, float]:
    """The cosine and sine coefficients of one harmonic of a motion.

    The diagnostic that matters for extended expansion.  A motion's *even*
    harmonics give the two up-and-downs per input revolution; its *odd* ones are
    what make the two halves differ, so the whole Atkinson asymmetry lives in
    the first harmonic and nowhere else.  A synthesised mechanism can therefore
    be read off directly: a first harmonic of zero is an Otto engine whatever
    else it does.
    """
    arr = _centred(lam)
    theta = np.linspace(0.0, 2.0 * np.pi, arr.size, endpoint=False)
    return (
        float(2.0 * np.mean(arr * np.cos(order * theta))),
        float(2.0 * np.mean(arr * np.sin(order * theta))),
    )


@dataclass(frozen=True)
class Assessment:
    """What a synthesised linkage actually delivers, against the requirements."""

    rms: float
    """Motion error against the target [mm]."""
    strain: float
    """Largest strain in a present member; a linkage runs at the stiffness floor."""
    four_phases: bool
    """Whether the motion is a four-stroke one at all."""
    expansion_stroke: float
    """``STE`` realised [mm]; ``nan`` when the motion is not a four-stroke."""
    compression_stroke: float
    """``STC`` realised [mm]; ``nan`` likewise."""
    compression_ratio: float
    """``epsilon`` implied by ``STC`` and the clearance volume."""
    first_harmonic: float
    """Amplitude of the first harmonic [mm] -- the asymmetry, and nothing else."""
    second_harmonic: float
    """Amplitude of the second harmonic [mm] -- the two up-and-downs."""
    present: tuple[tuple[str, str], ...]
    """Members that survived."""

    @property
    def asymmetry(self) -> float:
        """``STE - STC``: what makes the cycle extended-expansion rather than Otto."""
        return self.expansion_stroke - self.compression_stroke

    @property
    def is_extended_expansion(self) -> bool:
        """Whether the two strokes actually differ."""
        return bool(self.four_phases and self.asymmetry > 1.0)


def assess(
    layout: Layout,
    x: FloatArray,
    target: FloatArray,
    samples: int = 720,
) -> Assessment:
    """Score a synthesised design against what the engine needs.

    Deliberately separate from the objective.  The objective is a distance to a
    prescribed motion, which is what a gradient method can descend; this is the
    requirement, which is what a design has to meet.  Keeping them apart is what
    lets a run report that it minimised the one without satisfying the other.
    """
    from .constants import DEFAULT_SPEC
    from .cycle import PhaseError, find_phases

    motion = sweep(layout, x, samples=samples, penalty=PENALTY_SCHEDULE[-1])
    rho = layout.presences(x)
    names = layout.structure.names
    present = tuple(
        (names[i], names[j])
        for m, (i, j) in enumerate(layout.structure.members)
        if rho[m] > PRESENCE_THRESHOLD
    )
    try:
        phases = find_phases(motion.lam)
        ste, stc = phases.expansion_stroke, phases.compression_stroke
        four = True
    except PhaseError:
        ste = stc = float("nan")
        four = False
    ratio = 1.0 + stc * DEFAULT_SPEC.piston_area / DEFAULT_SPEC.dead_volume
    return Assessment(
        rms=motion_error(target, motion.lam),
        strain=motion.strain,
        four_phases=four,
        expansion_stroke=ste,
        compression_stroke=stc,
        compression_ratio=ratio,
        first_harmonic=float(np.hypot(*harmonic(motion.lam, 1))),
        second_harmonic=float(np.hypot(*harmonic(motion.lam, 2))),
        present=present,
    )


def format_assessment(assessment: Assessment) -> str:
    """A readable summary of a synthesised linkage."""
    lines = [
        f"members     {' '.join('-'.join(pair) for pair in assessment.present) or 'none'}",
        f"motion rms  {assessment.rms:8.3f} mm",
        f"peak strain {assessment.strain:8.2e}",
        f"harmonic 1  {assessment.first_harmonic:8.3f} mm   (the asymmetry)",
        f"harmonic 2  {assessment.second_harmonic:8.3f} mm   (the two up-and-downs)",
    ]
    if assessment.four_phases:
        lines += [
            f"STE         {assessment.expansion_stroke:8.3f} mm",
            f"STC         {assessment.compression_stroke:8.3f} mm",
            f"asymmetry   {assessment.asymmetry:8.3f} mm",
            f"epsilon     {assessment.compression_ratio:8.3f}",
        ]
    else:
        lines.append("            not a four-stroke motion")
    lines.append(
        "verdict     extended expansion"
        if assessment.is_extended_expansion
        else "verdict     symmetric strokes -- an Otto engine"
    )
    return "\n".join(lines)
