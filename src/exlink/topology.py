"""Synthesising the mechanism itself, rather than the dimensions of a given one.

What this is for
----------------
§5.4 measures the EXlink topology against a slider-crank and finds it worth
17.6 %.  Two topologies establish a contrast and not a trend, and §6.3 recorded
the missing third -- with the warning that choosing one by hand moves the
arbitrariness rather than removing it.  This module removes it instead, by
asking the optimizer for the topology.

The method is the second family of §2.2: the spring-connected model of
:cite:t:`kim2007spring`, in the reduced *link* form of :cite:t:`tran2024slm`,
with the gear-linkage extension of :cite:t:`yim2019gearlinkage`.  A design
domain is populated with candidate elements, every candidate is a spring whose
stiffness is a design variable, and the stiffnesses are penalised towards their
extremes.  A spring driven stiff is a rigid connection; a spring driven soft is
no connection at all.  One continuous parameterisation therefore covers every
mechanism the domain can hold, and the topology falls out of a gradient solve
rather than being assumed.

The synthesis problem
---------------------
The input is the **half-speed shaft**, turning once per cycle, and the output is
the piston height.  Over one input revolution the target motion has two equal
maxima and two unequal minima -- :func:`exlink.synthesis.target_motion` builds
the one that realises ``STE = 74`` mm and ``epsilon = 16`` exactly -- so the
synthesis is asked for a mechanism that turns one circular input into two
unequal up-and-downs of a point on the cylinder axis.

That period is the whole difficulty.  A mechanism driven by a shaft is periodic
in that shaft's angle, so one driven from the *crankshaft* alone cannot produce
a motion that differs between its two revolutions: extended expansion needs
something that turns once per cycle.  §5.6 sharpens this into an identity --
reach the piston from one shaft only and every odd harmonic vanishes exactly,
so the asymmetry, which lives entirely in the first harmonic, is *zero* rather
than small.  The piston must be reached from both shafts.

What the domain can hold
------------------------
Three kinds of candidate, each with one presence variable :math:`\\rho`:

**Bars.**  A spring between two nodes.

**Bodies.**  A rigid triangle on three nodes, carried as three springs sharing a
single presence.  §5.6 measured why this matters: assembled from three separate
bars a triangle is unreachable by descent, because each bar alone does nothing
and the objective is flat in its geometry until all three switch on together.
One variable turns the whole three-cornered link on, which is the shape the
identity above says extended expansion needs.

**Gears.**  A pair between two shafts, at a ratio from a small catalogue.  Every
shaft carries its own angle; the input's is prescribed and the rest are *free
coordinates of the equilibrium*, so a shaft with no gear on it is simply loose
and one with two is over-constrained -- both visible in the strain the answer
runs at.  A pair between shafts at centres :math:`c_i, c_j` meshes externally,
so its pitch radii are fixed by the centre distance and the ratio,

.. math:: r_i = \\frac{d\\,r}{1 + r}, \\qquad r_j = \\frac{d}{1 + r},
          \\qquad d = \\lVert c_i - c_j \\rVert,

and rolling without slip at the pitch point costs

.. math:: E = \\tfrac12 k(\\rho)\\,
          \\bigl(r_i \\vartheta_i + r_j \\vartheta_j - \\varphi\\bigr)^2 .

Writing it as slip rather than as an angular error is what keeps the mesh
geometrically consistent for free -- no separate centre-distance constraint is
needed, because the radii are read off the centres the search chose -- and keeps
the term in the same units as every other spring.  Choosing the ratio from a
catalogue of one presence per (pair, ratio) needs no new machinery either: two
ratios switched on at once is over-constraint, which the strain already detects.

The model
---------
Nodes carry planar coordinates.  Some are grounded, some are *pins* on shafts,
one is the piston -- abscissa on the cylinder axis, ordinate the output -- and
the rest are free.  Every spring :math:`m` joining nodes :math:`i` and :math:`j`
contributes

.. math:: E_m = \\tfrac12 k_m \\bigl(\\lVert x_i - x_j \\rVert - L_m\\bigr)^2,
          \\qquad k_m = k_{\\min} + (1 - k_{\\min})\\,\\rho_m^{\\,p}

with :math:`L_m` measured on the initial pose, so no design is prestressed.  The
configuration at an input angle minimises the total over the free coordinates --
free node positions *and* free shaft angles -- by a ridged Newton with a
backtracking line search, warm-started from the previous angle.  That is the
incremental analysis the published method uses, and the reason the branch a
mechanism runs on is followed rather than re-chosen at every angle.

The penalty exponent :math:`p` is raised on a schedule
(:data:`PENALTY_SCHEDULE`).  At :math:`p = 1` the problem is nearly convex in the
stiffnesses and the search moves freely; raising it makes intermediate
stiffnesses uneconomic, so the candidates separate into present and absent.

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
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_TARGETS

if TYPE_CHECKING:
    from .design import Design

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]

STIFFNESS_FLOOR = 1.0e-6
"""Stiffness of a fully absent candidate, relative to a fully present one.

Not zero, for the reason SIMP keeps a floor: a node reachable only through
absent springs would otherwise leave the Hessian singular, and the damped
Newton would have nothing to say about where it goes.
"""

PENALTY_SCHEDULE: tuple[float, ...] = (1.0, 2.0, 3.0, 4.0)
"""Exponents :math:`p` the continuation walks through."""

GROUND = -1
"""Node kind: coordinates fixed by the design vector."""
PIN = -2
"""Node kind: a point on a shaft's crank, carried round by the shaft's angle."""
SLIDER = -3
"""Node kind: abscissa fixed, ordinate free.  The piston."""
FREE = -4
"""Node kind: both coordinates free."""

DRIVEN = PIN
"""Old name for :data:`PIN`, kept so §5.6's vocabulary still reads."""

PRESENCE_THRESHOLD = 0.5
"""Presence above which a candidate is kept when the answer is made discrete."""


@dataclass(frozen=True)
class Shaft:
    """A grounded axis with an angle, which may be prescribed or solved for."""

    name: str
    centre_slot: int | None
    """Design-vector offset of its centre, or ``None`` when pinned at the origin."""
    ratio: float | None
    """Its speed as a multiple of the input's, or ``None`` when the angle is free.

    ``1`` is the input shaft itself.  ``-2`` is a shaft geared to it two to one
    and turning the other way, which is the relation :mod:`exlink.kinematics`
    gives the EXlink -- prescribed here rather than discovered, which is what
    §5.6's domain does and §5.7's does not.
    """
    angle_slot: int | None = None
    """Design-vector offset of its starting angle, for a free shaft."""

    @property
    def is_free(self) -> bool:
        """Whether the equilibrium solves for this shaft's angle."""
        return self.ratio is None


@dataclass(frozen=True)
class Pin:
    """A node carried round by a shaft at a radius and a phase."""

    node: int
    shaft: int
    slot: int
    """Design-vector offset of the radius, with the phase next to it."""


@dataclass(frozen=True)
class Element:
    """A candidate the search may switch on: one presence variable, one shape."""

    kind: str
    """``"bar"`` for two nodes, ``"body"`` for a rigid triangle on three."""
    nodes: tuple[int, ...]

    @property
    def springs(self) -> tuple[tuple[int, int], ...]:
        """The springs it is carried as: one for a bar, three for a body."""
        if self.kind == "bar":
            return ((self.nodes[0], self.nodes[1]),)
        i, j, k = self.nodes
        return ((i, j), (j, k), (k, i))


@dataclass(frozen=True)
class Gear:
    """A candidate gear pair between two shafts, at one catalogue ratio."""

    first: int
    second: int
    ratio: float
    """How many times faster the second shaft turns, and in the other sense."""
    phase_slot: int
    """Design-vector offset of its assembly phase."""


@dataclass(frozen=True)
class GroundStructure:
    """The candidate nodes, shafts and elements a mechanism may be built from.

    Instances are pure topology -- what exists and what may connect to what.
    Everything dimensional lives in the design vector, so one ground structure
    serves every design drawn on it.
    """

    names: tuple[str, ...]
    """Node names, in coordinate order."""
    kinds: tuple[int, ...]
    """One of :data:`GROUND`, :data:`PIN`, :data:`SLIDER`, :data:`FREE`."""
    shafts: tuple[Shaft, ...]
    pins: tuple[Pin, ...]
    elements: tuple[Element, ...]
    gears: tuple[Gear, ...] = ()

    @property
    def n_nodes(self) -> int:
        """Nodes in the structure."""
        return len(self.names)

    @property
    def n_members(self) -> int:
        """Candidate elements -- bars and bodies, not gears."""
        return len(self.elements)

    @property
    def n_presences(self) -> int:
        """Presence variables in all: one per element and one per gear."""
        return len(self.elements) + len(self.gears)

    @property
    def members(self) -> tuple[tuple[int, int], ...]:
        """The bar elements, as node pairs.  §5.6's vocabulary."""
        return tuple(e.nodes[:2] for e in self.elements if e.kind == "bar")  # type: ignore[misc]

    @property
    def free_shafts(self) -> tuple[int, ...]:
        """Indices of the shafts whose angle the equilibrium solves for."""
        return tuple(i for i, s in enumerate(self.shafts) if s.is_free)

    def springs(self) -> tuple[tuple[int, int], ...]:
        """Every spring in the domain, expanded from the elements."""
        return tuple(s for element in self.elements for s in element.springs)

    def spring_owner(self) -> IntArray:
        """Which presence variable each spring takes its stiffness from."""
        owner = [m for m, element in enumerate(self.elements) for _ in element.springs]
        return np.asarray(owner, dtype=np.int64)

    def free_mask(self) -> NDArray[np.bool_]:
        """Which of the ``2 n_nodes`` coordinates the equilibrium solves for."""
        mask = np.zeros(2 * self.n_nodes, dtype=bool)
        for i, kind in enumerate(self.kinds):
            if kind == FREE:
                mask[2 * i : 2 * i + 2] = True
            elif kind == SLIDER:
                mask[2 * i + 1] = True
        return mask

    @property
    def n_reduced(self) -> int:
        """Unknowns per equilibrium: free coordinates, then free shaft angles."""
        return int(self.free_mask().sum()) + len(self.free_shafts)

    def ends(self) -> tuple[IntArray, IntArray]:
        """The springs' first and second node indices, as arrays."""
        pairs = np.asarray(self.springs(), dtype=np.int64).reshape(-1, 2)
        return pairs[:, 0], pairs[:, 1]

    def index(self, name: str) -> int:
        """Index of a node by name."""
        return self.names.index(name)


@dataclass(frozen=True)
class Layout:
    """Where each quantity sits in the design vector.

    The vector holds the *initial* configuration, the gear phases and the
    presences, and nothing else.  Element lengths are deliberately absent: each
    is the distance between its endpoints in the initial pose, so every design is
    assembled without prestress and an entire class of structures that could
    never be built is never visited.
    """

    structure: GroundStructure
    ground_slots: tuple[tuple[int, int], ...]
    """``(node, offset)`` for every grounded node whose position is free."""
    free_slots: tuple[tuple[int, int], ...]
    """``(node, offset)`` for every node whose initial position is chosen."""
    slider_slot: tuple[int, int]
    """``(node, offset)`` of the piston: cylinder abscissa then initial height."""
    presence_offset: int
    """First index of the presences."""
    driven_slots: tuple[tuple[int, int], ...] = ()
    """``(node, offset)`` of each pin's radius, with its phase next to it."""
    centre_of: dict[int, int] = field(default_factory=dict)
    """Grounded node each pin turns about.  §5.6's vocabulary."""

    @property
    def size(self) -> int:
        """Length of the design vector."""
        return self.presence_offset + self.structure.n_presences

    def presences(self, x: FloatArray) -> FloatArray:
        """The presence variables of a design vector."""
        return np.asarray(x[self.presence_offset :], dtype=float)

    def element_presences(self, x: FloatArray) -> FloatArray:
        """The presences of the bars and bodies alone."""
        return self.presences(x)[: self.structure.n_members]

    def gear_presences(self, x: FloatArray) -> FloatArray:
        """The presences of the candidate gear pairs."""
        return self.presences(x)[self.structure.n_members :]


def shaft_centres(layout: Layout, x: FloatArray) -> FloatArray:
    """Where each shaft's axis sits."""
    centres = np.zeros((len(layout.structure.shafts), 2), dtype=float)
    for s, shaft in enumerate(layout.structure.shafts):
        if shaft.centre_slot is not None:
            centres[s] = x[shaft.centre_slot : shaft.centre_slot + 2]
    return centres


def _shaft_angles(layout: Layout, x: FloatArray, theta: float, free: FloatArray) -> FloatArray:
    """Every shaft's angle: prescribed from the input, or taken from the state."""
    angles = np.zeros(len(layout.structure.shafts), dtype=float)
    k = 0
    for s, shaft in enumerate(layout.structure.shafts):
        if shaft.ratio is None:
            angles[s] = free[k]
            k += 1
        else:
            angles[s] = shaft.ratio * theta
    return angles


def _base_positions(layout: Layout, x: FloatArray) -> FloatArray:
    """Grounded coordinates, and the piston's abscissa."""
    pos = np.zeros((layout.structure.n_nodes, 2), dtype=float)
    for node, offset in layout.ground_slots:
        pos[node] = x[offset : offset + 2]
    node, offset = layout.slider_slot
    pos[node, 0] = float(x[offset])
    return pos


def _place(
    layout: Layout, x: FloatArray, theta: float, reduced: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """Node coordinates at one input angle, and each pin's tangent.

    The tangent is how the pin moves when its shaft turns, and it is what makes a
    free shaft angle a coordinate the equilibrium can solve for rather than a
    parameter it is handed.
    """
    structure = layout.structure
    mask = structure.free_mask()
    n_coords = int(mask.sum())
    centres = shaft_centres(layout, x)
    angles = _shaft_angles(layout, x, theta, reduced[n_coords:])

    pos = _base_positions(layout, x)
    tangents = np.zeros((structure.n_nodes, 2), dtype=float)
    for pin in structure.pins:
        radius, phase = float(x[pin.slot]), float(x[pin.slot + 1])
        carried = angles[pin.shaft] + phase
        direction = np.array([np.cos(carried), np.sin(carried)])
        pos[pin.node] = centres[pin.shaft] + radius * direction
        tangents[pin.node] = radius * np.array([-np.sin(carried), np.cos(carried)])
    pos.reshape(-1)[mask] = reduced[:n_coords]
    return pos, tangents


def initial_configuration(layout: Layout, x: FloatArray) -> FloatArray:
    """The full node layout at zero input angle, free nodes included."""
    structure = layout.structure
    mask = structure.free_mask()
    reduced = np.zeros(structure.n_reduced, dtype=float)
    seed = np.zeros(2 * structure.n_nodes, dtype=float)
    for node, offset in layout.free_slots:
        seed[2 * node : 2 * node + 2] = x[offset : offset + 2]
    node, offset = layout.slider_slot
    seed[2 * node + 1] = float(x[offset + 1])
    reduced[: int(mask.sum())] = seed[mask]
    for k, s in enumerate(structure.free_shafts):
        slot = structure.shafts[s].angle_slot
        if slot is not None:
            reduced[int(mask.sum()) + k] = float(x[slot])
    return _place(layout, x, 0.0, reduced)[0]


def member_lengths(layout: Layout, x: FloatArray) -> FloatArray:
    """Rest length of every spring in the domain, from the initial pose."""
    pos = initial_configuration(layout, x)
    springs = layout.structure.springs()
    out = np.empty(len(springs), dtype=float)
    for m, (i, j) in enumerate(springs):
        out[m] = float(np.hypot(*(pos[i] - pos[j])))
    return out


def stiffnesses(presence: FloatArray, penalty: float) -> FloatArray:
    """SIMP stiffness of each candidate: a floor plus its penalised presence."""
    rho = np.clip(np.asarray(presence, dtype=float), 0.0, 1.0)
    return STIFFNESS_FLOOR + (1.0 - STIFFNESS_FLOOR) * rho**penalty


def gear_radii(layout: Layout, x: FloatArray) -> FloatArray:
    """Pitch radii of every candidate pair, from its centres and its ratio.

    An external mesh puts the pitch point between the two axes, so the radii sum
    to the centre distance and their ratio is the speed ratio.  Reading them off
    the centres the search chose is what makes the mesh consistent for free: no
    separate centre-distance constraint is needed, because a pair drawn this way
    always fits the shafts it is drawn between.
    """
    centres = shaft_centres(layout, x)
    out = np.zeros((len(layout.structure.gears), 2), dtype=float)
    for g, gear in enumerate(layout.structure.gears):
        distance = float(np.hypot(*(centres[gear.first] - centres[gear.second])))
        out[g] = (distance * gear.ratio / (1.0 + gear.ratio), distance / (1.0 + gear.ratio))
    return out


def _energy_terms(
    pos: FloatArray,
    ends: tuple[IntArray, IntArray],
    k: FloatArray,
    rest: FloatArray,
) -> tuple[float, FloatArray, FloatArray]:
    """Total spring energy and its first two derivatives in the coordinates.

    The blocks are the textbook ones for an axial spring, assembled over the
    full coordinate set with the caller restricting to the free ones.  It is
    written over all springs at once because it is the innermost thing in the
    method -- every Newton step of every angle of every finite difference -- and
    a Python loop over springs costs more than the arithmetic does.

    Args:
        pos: ``(n_nodes, 2)`` coordinates.
        ends: The springs' first and second node indices.
        k: Spring stiffnesses.
        rest: Spring rest lengths.

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


def _jacobian(layout: Layout, tangents: FloatArray) -> FloatArray:
    """How the node coordinates move with each reduced unknown."""
    structure = layout.structure
    mask = structure.free_mask()
    n_coords = int(mask.sum())
    jac = np.zeros((2 * structure.n_nodes, structure.n_reduced), dtype=float)
    jac[np.flatnonzero(mask), np.arange(n_coords)] = 1.0
    for k, s in enumerate(structure.free_shafts):
        for pin in structure.pins:
            if pin.shaft == s:
                jac[2 * pin.node : 2 * pin.node + 2, n_coords + k] = tangents[pin.node]
    return jac


def _gear_terms(
    layout: Layout,
    x: FloatArray,
    theta: float,
    reduced: FloatArray,
    stiffness: FloatArray,
    radii: FloatArray,
) -> tuple[float, FloatArray, FloatArray, FloatArray]:
    """Slip energy of every candidate pair, in the reduced coordinates.

    Returns the energy, its reduced gradient and Hessian, and the slip of each
    pair as a fraction of its centre distance -- the gear's own strain, and the
    reading that says whether a shaft has ended up with two pairs on it.
    """
    structure = layout.structure
    n = structure.n_reduced
    n_coords = n - len(structure.free_shafts)
    angles = _shaft_angles(layout, x, theta, reduced[n_coords:])
    column = {s: n_coords + k for k, s in enumerate(structure.free_shafts)}

    total = 0.0
    grad = np.zeros(n, dtype=float)
    hess = np.zeros((n, n), dtype=float)
    slip = np.zeros(len(structure.gears), dtype=float)
    for g, gear in enumerate(structure.gears):
        first, second = radii[g]
        residual = (
            first * angles[gear.first]
            + second * angles[gear.second]
            - float(x[gear.phase_slot])
        )
        total += 0.5 * stiffness[g] * residual**2
        span = first + second
        slip[g] = abs(residual) / span if span > 0.0 else 0.0
        coefficients = {}
        if gear.first in column:
            coefficients[column[gear.first]] = first
        if gear.second in column:
            coefficients[column[gear.second]] = second
        for a, ca in coefficients.items():
            grad[a] += stiffness[g] * residual * ca
            for b, cb in coefficients.items():
                hess[a, b] += stiffness[g] * ca * cb
    return total, grad, hess, slip


NEWTON_STEPS = 30
"""Newton iterations allowed per input angle."""
NEWTON_TOLERANCE = 1.0e-9
"""Convergence test on the reduced gradient, relative to its scale.

The gradient of a spring energy carries the stiffness and a length, and both
vary by orders of magnitude across a penalisation schedule, so an absolute
threshold would be met instantly at one end of it and never at the other.
"""
RIDGE = 1.0e-9
"""Fraction of the mean diagonal added to the Hessian before solving.

The Hessian is singular exactly when some node or shaft is held only by
candidates the search has switched off -- which, in a topology optimization, is
most of the time.  The ridge leaves such a coordinate where the warm start put
it instead of sending it to infinity, and is small enough not to disturb one
that is properly held.
"""
BACKTRACK_STEPS = 20
"""Halvings allowed before a Newton step is abandoned."""

SLACK_TOLERANCE = 1.0e-3
"""Component in a null direction above which the piston counts as undetermined."""

SINGULAR_TOLERANCE = 1.0e-8
"""Singular values below this fraction of the largest count as a null direction."""


def _output_slack(hess: FloatArray, slider_row: int) -> float:
    """How much of the Hessian's null space the piston's coordinate lives in.

    A spring network holds a coordinate only through the candidates the search
    switched on.  Where it holds none, the reduced Hessian is singular in that
    direction and the equilibrium is not unique -- so the motion reported is the
    one the ridge and the warm start happened to pick.  Projecting the slider's
    coordinate onto those directions says whether the *output* is among them,
    which is the only part that matters: an unused node floating off on its own
    is harmless, and a piston floating off on its own is the whole answer.
    """
    if hess.size == 0:
        return 1.0
    _u, singular, right = np.linalg.svd(hess)
    largest = float(singular[0]) if singular.size else 0.0
    if largest <= 0.0:
        return 1.0
    null = right[singular <= SINGULAR_TOLERANCE * largest]
    if null.size == 0:
        return 0.0
    return float(np.max(np.abs(null[:, slider_row])))


def _assemble(
    layout: Layout,
    x: FloatArray,
    theta: float,
    reduced: FloatArray,
    k: FloatArray,
    rest: FloatArray,
    gear_stiffness: FloatArray,
    radii: FloatArray,
    ends: tuple[IntArray, IntArray] | None = None,
) -> tuple[float, FloatArray, FloatArray, FloatArray, FloatArray]:
    """Total energy and its reduced derivatives, springs and gears together.

    The reduction is a change of variables rather than a restriction: a free
    shaft angle moves several pins at once, so its column of the Jacobian
    carries a tangent at each of them, and its diagonal picks up the curvature
    of the circle each pin runs on.
    """
    structure = layout.structure
    pos, tangents = _place(layout, x, theta, reduced)
    energy, grad_pos, hess_pos = _energy_terms(
        pos, structure.ends() if ends is None else ends, k, rest
    )

    jac = _jacobian(layout, tangents)
    grad = jac.T @ grad_pos
    hess = jac.T @ hess_pos @ jac

    n_coords = structure.n_reduced - len(structure.free_shafts)
    centres = shaft_centres(layout, x)
    for column, s in enumerate(structure.free_shafts):
        curvature = 0.0
        for pin in structure.pins:
            if pin.shaft != s:
                continue
            radial = pos[pin.node] - centres[pin.shaft]
            curvature -= float(np.dot(grad_pos[2 * pin.node : 2 * pin.node + 2], radial))
        hess[n_coords + column, n_coords + column] += curvature

    gear_energy, gear_grad, gear_hess, slip = _gear_terms(
        layout, x, theta, reduced, gear_stiffness, radii
    )
    return energy + gear_energy, grad + gear_grad, hess + gear_hess, pos, slip


def _equilibrate(
    layout: Layout,
    x: FloatArray,
    theta: float,
    guess: FloatArray,
    k: FloatArray,
    rest: FloatArray,
    gear_stiffness: FloatArray,
    radii: FloatArray,
    ends: tuple[IntArray, IntArray] | None = None,
) -> tuple[FloatArray, bool]:
    """Minimise the energy over the free coordinates and shaft angles.

    Newton with a ridge and a backtracking line search.  The ridge handles the
    singular Hessian of a partly switched-off mechanism; the line search handles
    the fact that the energy is not convex in the coordinates once a spring is in
    compression or a shaft has turned past its pin's tangent, which is where a
    bare Newton step overshoots.

    Args:
        layout: The design vector's layout.
        x: Design vector.
        theta: Input angle [rad].
        guess: Reduced coordinates to start from, normally the previous angle's.
        k: Spring stiffnesses.
        rest: Spring rest lengths.
        gear_stiffness: Stiffness of each candidate pair.
        radii: Pitch radii of each candidate pair.

    Returns:
        The reduced coordinates and whether the gradient test was met.
    """
    scale = max(1.0, float(np.max(k, initial=0.0)) * float(np.max(rest, initial=1.0)))
    tol = NEWTON_TOLERANCE * scale
    eye = np.eye(layout.structure.n_reduced)

    u = np.array(guess, dtype=float)
    for _ in range(NEWTON_STEPS):
        energy, grad, hess, _pos, _slip = _assemble(
            layout, x, theta, u, k, rest, gear_stiffness, radii, ends
        )
        if float(np.max(np.abs(grad))) < tol:
            return u, True
        ridge = RIDGE * max(float(np.mean(np.abs(np.diag(hess)))), 1.0e-30)
        try:
            step = np.linalg.solve(hess + ridge * eye, -grad)
        except np.linalg.LinAlgError:
            return u, False
        alpha = 1.0
        for _back in range(BACKTRACK_STEPS):
            trial = u + alpha * step
            if (
                _assemble(layout, x, theta, trial, k, rest, gear_stiffness, radii, ends)[0]
                <= energy
            ):
                u = trial
                break
            alpha *= 0.5
        else:
            return u, False
    final = _assemble(layout, x, theta, u, k, rest, gear_stiffness, radii, ends)[1]
    return u, float(np.max(np.abs(final))) < tol


@dataclass(frozen=True)
class Motion:
    """What a candidate mechanism does over one input revolution."""

    lam: FloatArray
    """Piston height at each input angle [mm]."""
    converged: bool
    """Whether every angle reached equilibrium."""
    strain: float
    """Largest strain over the revolution, as a fraction of a rest length.

    The honest measure of whether the springs are behaving as rigid links, and
    it covers the gears too: a pair's strain is its slip at the pitch point as a
    fraction of the centre distance.  A mechanism the penalisation has resolved
    runs at a strain of order the stiffness floor; anything larger means the
    answer is a deforming structure rather than a mechanism.
    """
    gear_ratios: tuple[float, ...] = ()
    """Ratio of every pair the design kept, in catalogue order."""
    output_slack: float = 0.0
    """How undetermined the piston is, from 0 (fixed by the input) to 1 (free).

    Strain catches an *over*-constrained answer.  This catches the opposite, and
    the opposite is the one that lies.  If the equilibrium leaves a direction the
    energy does not resist -- a shaft with no gear on it, a chain with a link
    missing -- then where the mechanism goes is chosen by the ridge and the warm
    start rather than by the mechanism, and the motion that comes back is an
    artefact that can look like anything, including extended expansion.  The
    reading is the largest component the slider's coordinate has in the null
    space of the reduced Hessian, over the revolution.  A mechanism whose motion
    the input actually determines runs at zero; anything above
    :data:`SLACK_TOLERANCE` is not a mechanism.
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
        The piston motion and its diagnostics.
    """
    structure = layout.structure
    element_rho = layout.element_presences(x)
    gear_rho = layout.gear_presences(x)
    owner = structure.spring_owner()
    k = stiffnesses(element_rho, penalty)[owner]
    spring_rho = element_rho[owner]
    gear_stiffness = stiffnesses(gear_rho, penalty)
    radii = gear_radii(layout, x)
    rest = member_lengths(layout, x)
    safe_rest = np.where(rest > 0.0, rest, 1.0)
    first, second = structure.ends()

    # Built once and threaded through, rather than rebuilt inside every Newton
    # step of every angle -- which is where an enumeration over topologies would
    # otherwise spend most of its time.
    ends = (first, second)
    solve_k, solve_rest = k, rest

    start = initial_configuration(layout, x)
    mask = structure.free_mask()
    u = np.zeros(structure.n_reduced, dtype=float)
    u[: int(mask.sum())] = start.reshape(-1)[mask]
    for column, s in enumerate(structure.free_shafts):
        slot = structure.shafts[s].angle_slot
        if slot is not None:
            u[int(mask.sum()) + column] = float(x[slot])

    angles = np.linspace(0.0, 2.0 * np.pi, int(samples), endpoint=False)
    lam = np.empty(angles.size, dtype=float)
    slider_node, _slider_offset = layout.slider_slot
    slider_row = int(np.flatnonzero(mask).tolist().index(2 * slider_node + 1))
    ok = True
    worst = 0.0
    slack = 0.0
    for step, theta in enumerate(angles):
        u, converged = _equilibrate(
            layout, x, float(theta), u, solve_k, solve_rest, gear_stiffness, radii, ends
        )
        ok = ok and converged
        hess = _assemble(
            layout, x, float(theta), u, solve_k, solve_rest, gear_stiffness, radii, ends
        )[2]
        slack = max(slack, _output_slack(hess, slider_row))
        pos, _tangents = _place(layout, x, float(theta), u)
        lam[step] = pos[slider_node, 1]
        gap = pos[first] - pos[second]
        length = np.sqrt(np.einsum("mi,mi->m", gap, gap))
        strain = np.abs(length - rest) / safe_rest
        worst = max(worst, float(np.max(strain * spring_rho, initial=0.0)))
        slip = _gear_terms(layout, x, float(theta), u, gear_stiffness, radii)[3]
        worst = max(worst, float(np.max(slip * gear_rho, initial=0.0)))
    kept = tuple(
        gear.ratio for g, gear in enumerate(structure.gears) if gear_rho[g] > PRESENCE_THRESHOLD
    )
    return Motion(lam=lam, converged=ok, strain=worst, gear_ratios=kept, output_slack=slack)


def engine_ground_structure() -> tuple[GroundStructure, Layout]:
    """§5.6's candidate set: bars only, and the gear ratio given.

    Eight nodes and twelve candidate bars.  The input shaft is grounded at the
    origin -- a gauge choice, since translating the whole mechanism changes
    nothing -- and carries the pin ``P1``.  A second grounded shaft carries
    ``P2`` at twice the input speed and the opposite sense, which is the relation
    a 2:1 gear pair imposes and the one EXlink uses; the search is free to leave
    it out by switching off every bar that touches it, but not to choose it.
    ``C`` is a plain grounded pivot, ``F1`` and ``F2`` are free nodes, and ``S``
    is the piston, held on the cylinder axis.

    Bars between two prescribed nodes are excluded: their length is a function
    of the input angle alone, so such a bar is a constraint on the design rather
    than a degree of freedom for the mechanism.

    Returns:
        The structure and the layout of the design vector drawn on it.
    """
    names = ("A", "P1", "B", "P2", "C", "F1", "F2", "S")
    kinds = (GROUND, PIN, GROUND, PIN, GROUND, FREE, FREE, SLIDER)
    idx = {name: i for i, name in enumerate(names)}
    movable = ("P1", "P2", "C", "F1", "F2", "S")
    elements = tuple(
        Element("bar", (idx[a], idx[b]))
        for i, a in enumerate(movable)
        for b in movable[i + 1 :]
        if kinds[idx[a]] in (FREE, SLIDER) or kinds[idx[b]] in (FREE, SLIDER)
    )
    structure = GroundStructure(
        names=names,
        kinds=kinds,
        shafts=(
            Shaft("input", centre_slot=None, ratio=1.0),
            Shaft("geared", centre_slot=0, ratio=-2.0),
        ),
        pins=(Pin(idx["P1"], 0, 4), Pin(idx["P2"], 1, 6)),
        elements=elements,
    )
    layout = Layout(
        structure=structure,
        ground_slots=((idx["B"], 0), (idx["C"], 2)),
        free_slots=((idx["F1"], 8), (idx["F2"], 10)),
        slider_slot=(idx["S"], 12),
        presence_offset=14,
        driven_slots=((idx["P1"], 4), (idx["P2"], 6)),
        centre_of={idx["P1"]: idx["A"], idx["P2"]: idx["B"]},
    )
    return structure, layout


GEAR_CATALOGUE: tuple[float, ...] = (1.0, 2.0, 3.0, 4.0)
"""Speed ratios a candidate pair may be cut for.

Whole numbers, and that is a requirement rather than a convenience.  The cycle
has to *close*: after one revolution of the input shaft the mechanism must be
back where it started, or the engine's motion differs between successive cycles
and there is no four-stroke to speak of.  A geared shaft turning :math:`r` times
per input revolution returns to its start only when :math:`r` is an integer, so
a three-to-two pair would give a mechanism whose period is two input
revolutions.  That is a real machine and not this one.

Small on purpose otherwise.  A ratio is one integer tooth pair and one module
away from being buildable, and §4.4 already owns that step; what this catalogue
has to settle is only whether the *synthesis* wants a two-to-one relation, which
is the question §5.6 could not ask because it was handed one.
"""


def geared_ground_structure() -> tuple[GroundStructure, Layout]:
    """§5.7's candidate set: bodies as well as bars, and the gear chosen.

    Three changes from :func:`engine_ground_structure`, each aimed at something
    §5.6 measured.

    The second shaft's angle is **free**, and four candidate pairs -- one per
    catalogue ratio -- may drive it from the input.  So the ratio is an outcome:
    a design that keeps the two-to-one pair has rediscovered EXlink's gearing,
    one that keeps another has found something else, and one that keeps none has
    a loose shaft its linkage must hold.

    Ten candidate **bodies** join the twelve bars: rigid triangles carried as
    three springs sharing a single presence.  §5.6's stall was a triangle that
    could only be assembled from three separately-switched bars, each useless
    alone; one variable now turns a whole three-cornered link on.

    Returns:
        The structure and the layout of the design vector drawn on it.
    """
    names = ("A", "P1", "B", "P2", "C", "F1", "F2", "S")
    kinds = (GROUND, PIN, GROUND, PIN, GROUND, FREE, FREE, SLIDER)
    idx = {name: i for i, name in enumerate(names)}
    movable = ("P1", "P2", "C", "F1", "F2", "S")
    loose = {idx["F1"], idx["F2"], idx["S"]}

    bars = tuple(
        Element("bar", (idx[a], idx[b]))
        for i, a in enumerate(movable)
        for b in movable[i + 1 :]
        if kinds[idx[a]] in (FREE, SLIDER) or kinds[idx[b]] in (FREE, SLIDER)
    )
    triples = tuple(
        Element("body", (idx[a], idx[b], idx[c]))
        for i, a in enumerate(movable)
        for j, b in enumerate(movable[i + 1 :], start=i + 1)
        for c in movable[j + 1 :]
        if len(loose & {idx[a], idx[b], idx[c]}) >= 2
    )
    gears = tuple(
        Gear(first=0, second=1, ratio=ratio, phase_slot=15 + g)
        for g, ratio in enumerate(GEAR_CATALOGUE)
    )
    structure = GroundStructure(
        names=names,
        kinds=kinds,
        shafts=(
            Shaft("input", centre_slot=None, ratio=1.0),
            Shaft("geared", centre_slot=0, ratio=None, angle_slot=14),
        ),
        pins=(Pin(idx["P1"], 0, 4), Pin(idx["P2"], 1, 6)),
        elements=bars + triples,
        gears=gears,
    )
    layout = Layout(
        structure=structure,
        ground_slots=((idx["B"], 0), (idx["C"], 2)),
        free_slots=((idx["F1"], 8), (idx["F2"], 10)),
        slider_slot=(idx["S"], 12),
        presence_offset=15 + len(GEAR_CATALOGUE),
        driven_slots=((idx["P1"], 4), (idx["P2"], 6)),
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
    for pin in layout.structure.pins:
        lower[pin.slot], upper[pin.slot] = 5.0, 120.0
        lower[pin.slot + 1], upper[pin.slot + 1] = -np.pi, np.pi
    for _node, offset in layout.free_slots:
        lower[offset : offset + 2] = -span
        upper[offset : offset + 2] = span
    _node, offset = layout.slider_slot
    lower[offset : offset + 2] = -span
    upper[offset : offset + 2] = span
    for shaft in layout.structure.shafts:
        if shaft.angle_slot is not None:
            lower[shaft.angle_slot], upper[shaft.angle_slot] = -np.pi, np.pi
    for gear in layout.structure.gears:
        lower[gear.phase_slot], upper[gear.phase_slot] = -span * np.pi, span * np.pi
    lower[layout.presence_offset :] = 0.0
    upper[layout.presence_offset :] = 1.0
    return lower, upper


COUNT_WEIGHT = 0.05
"""Millimetres of motion error a whole extra candidate is worth.

A tie-break and nothing more.  It began at 1.0 mm, which was a mistake worth
recording: a design with no elements at all leaves the piston free, the piston
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
"""Millimetres of motion error a unit of strain is worth.

This is the term that decides whether the answer is a mechanism, and it began
an order of magnitude too small.  At 50 a start reached a motion error of
2.0 mm while carrying 3 % strain, and paid 1.5 mm for it -- so the search was
being *paid* to leave candidates half present, because a half-stiff member lets
an overconstrained network move by giving rather than by articulating.  The gray
design fitted the target well and rounded to nothing.  At 800 the same 3 % costs
25 mm and softness stops being a way to buy motion.
"""

HARMONIC_WEIGHT = 1.0
"""Weight on the *relative* error of each harmonic, against the plain distance.

The repair for the stall §5.6 measured.  The target carries 32.33 mm at its
second harmonic and 9.02 mm at its first, and the whole Atkinson asymmetry is
the first; but a plain root-mean-square is dominated by the second, so a
mechanism that reproduces the second and none of the first already scores
7.19 mm against the 23.74 mm of doing nothing.  That is a deep local minimum
holding an Otto engine.  Charging each harmonic by its error *relative to the
target's own amplitude there* makes missing the small one cost as much,
proportionally, as missing the large one.
"""

BRIDGE_WEIGHT = 30.0
"""Millimetres charged for a *driven* pin the linkage does not carry to the piston.

The condition §5.6 *derives* rather than assumes: a piston reached from one
shaft alone is periodic in that shaft's angle, so every odd harmonic vanishes
exactly and the answer is an Otto engine whatever else it does.  Extended
expansion therefore needs both pins joined to the piston, and both searches
spent their budget on designs where one of them was not -- start 2 of §5.7 built
two three-cornered bodies and hung all of them off the geared shaft.

Imposing a proven necessary condition is not assuming the answer: nothing here
says *how* the two chains meet, how many elements it takes, or what shape the
body that joins them has.

The charge is weighted by how far each pin's shaft is *driven*
(:func:`driven_fraction`), which is what keeps it honest in both directions.  A
pin on a shaft no gear turns is not a second source, so reaching it earns
nothing -- and not reaching it costs nothing either, because a topology with no
gear at all is a legitimate single-input linkage rather than a broken two-input
one.  The enumeration of :func:`enumerate_mechanisms` contains both families.
"""

SLACK_WEIGHT = 200.0
"""Millimetres charged for a piston the equilibrium does not determine.

The other half of §5.7's finding.  Strain charges an over-constrained answer;
without this, nothing charges an under-constrained one, and an under-constrained
one can report any motion at all -- start 3 of §5.7 reported 70.9 mm of
asymmetry from a mechanism whose geared shaft nothing was driving.
"""

HARMONIC_ORDERS: tuple[int, ...] = (1, 2)
"""Harmonics the relative term is taken over: the asymmetry, and the strokes."""

HARMONIC_FLOOR = 1.0e-6
"""Amplitude, relative to the target's own scale, below which an order is skipped.

A relative error needs something to be relative *to*.  A target with no content
at an order still has a coefficient of order 1e-15 there rather than exactly
zero, and dividing by that turns a bounded objective into 1e16 -- which is what
happened the first time this was run against a pure second-harmonic target.
"""

DISCRETENESS_SCHEDULE: tuple[float, ...] = (0.0, 2.0, 6.0, 18.0)
r"""Millimetres charged for a fully undecided candidate, one per penalty rung.

The SIMP exponent alone did not separate them here.  Ramping an explicit
:math:`4\rho(1-\rho)` charge alongside it is the standard repair, and the ramp
matters: charged from the start it decides the topology before the geometry is
any good, and never charged at all it leaves an answer that cannot be rounded.
"""


@dataclass(frozen=True)
class Candidate:
    """A mechanism the synthesis produced, and how well it does."""

    x: FloatArray
    """Design vector."""
    rms: float
    """Motion error against the target, mean removed [mm]."""
    strain: float
    """Largest strain in a present candidate (:attr:`Motion.strain`)."""
    discreteness: float
    """Mean of :math:`4\\rho(1-\\rho)`: zero when every candidate has resolved."""
    present: tuple[tuple[str, ...], ...]
    """Elements that survived, by node name."""
    gear_ratios: tuple[float, ...]
    """Ratios of the gear pairs that survived."""
    motion: FloatArray
    """The piston motion it realises [mm]."""
    converged: bool
    """Whether every equilibrium solve met its gradient test."""

    @property
    def uses_second_shaft(self) -> bool:
        """Whether any surviving element touches the geared shaft's pin."""
        return any("P2" in nodes for nodes in self.present)


def _centred(values: FloatArray) -> FloatArray:
    """Values with their mean removed, so an offset is not an error."""
    arr = np.asarray(values, dtype=float)
    return arr - float(arr.mean())


def _resample(values: FloatArray, size: int) -> FloatArray:
    """Put a motion on a given number of uniformly spaced angles."""
    arr = np.asarray(values, dtype=float)
    if arr.size == size:
        return arr
    grid = np.linspace(0.0, 1.0, size, endpoint=False)
    return np.interp(grid, np.linspace(0.0, 1.0, arr.size, endpoint=False), arr, period=1.0)


def best_datum(target: FloatArray, lam: FloatArray) -> int:
    """The rotation of the input datum that lines a motion up with the target.

    Where :math:`\\theta_1 = 0` sits is a choice, not a property of a mechanism:
    turn the whole design about the input axis and every crank phase moves with
    it, giving the same engine with its cycle starting somewhere else.  A
    comparison that fixes the datum therefore charges a design for something it
    was never asked to control, and it charges a lot -- the real EXlink sits
    170 degrees from the target's datum and was scored at 13.97 mm as posed
    against 4.05 mm once aligned, which put it *behind* a degenerate answer that
    had tuned its phase.

    Minimising the squared difference over the shift is the same as maximising
    the circular cross-correlation, so one transform pair settles it.
    """
    a, b = _centred(target), _centred(lam)
    spectrum = np.fft.rfft(a) * np.conj(np.fft.rfft(b))
    return int(np.argmax(np.fft.irfft(spectrum, a.size)))


def align(target: FloatArray, lam: FloatArray) -> FloatArray:
    """A motion resampled onto the target's angles and rolled onto its datum."""
    resampled = _resample(_centred(lam), np.asarray(target).size)
    return np.roll(resampled, best_datum(target, resampled))


def motion_error(target: FloatArray, lam: FloatArray) -> float:
    """Root-mean-square distance between two motions, offset and datum removed.

    The target's mean height is arbitrary -- it is set by where the cylinder is
    bolted, not by the mechanism -- and so is the crank angle its cycle is drawn
    from, for the reason :func:`best_datum` gives.  Both are removed before
    anything is charged.
    """
    a = _centred(target)
    return float(np.sqrt(np.mean((a - align(a, lam)) ** 2)))


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


def harmonic_error(target: FloatArray, lam: FloatArray) -> float:
    """Millimetres charged for missing each harmonic, in proportion to its size.

    Each order's coefficient error is divided by the target's own amplitude
    there and the mean is scaled back into millimetres by the target's standard
    deviation, so a design that misses the first harmonic entirely is charged
    the same as one that misses the second entirely -- which a plain distance
    does not do, and which is why §5.6 stalled on an Otto engine.

    Taken on the *aligned* motion, for the reason :func:`best_datum` gives: a
    harmonic coefficient carries a phase, and rotating the input datum by
    :math:`s` turns the :math:`n`-th of them by :math:`ns`, so comparing
    coefficients at a fixed datum charges a design an order-one error per
    harmonic for a quantity that is not its business.
    """
    reference = float(np.std(_centred(target)))
    resampled = align(target, lam)
    errors = []
    for order in HARMONIC_ORDERS:
        wanted = np.array(harmonic(target, order))
        got = np.array(harmonic(resampled, order))
        amplitude = float(np.hypot(*wanted))
        if amplitude <= HARMONIC_FLOOR * reference:
            continue
        errors.append(float(np.hypot(*(got - wanted))) / amplitude)
    return reference * float(np.mean(errors)) if errors else 0.0


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
    """What the search minimises: the motion, and four deterrents.

    Args:
        layout: The design vector's layout.
        x: Design vector.
        target: Target piston motion.
        samples: Input angles per sweep.
        penalty: SIMP exponent on the presences.
        discreteness: Millimetres charged for a fully undecided candidate.

    Returns:
        The objective, in millimetres of equivalent motion error.
    """
    motion = sweep(layout, x, samples=samples, penalty=penalty)
    if not np.all(np.isfinite(motion.lam)):
        return 1.0e6
    rho = layout.presences(x)
    driven = driven_fraction(layout, x)
    unbridged = float(np.sum(driven * (1.0 - reachability(layout, x))))
    return (
        motion_error(target, motion.lam)
        + HARMONIC_WEIGHT * harmonic_error(target, motion.lam)
        + TRAVEL_WEIGHT * travel_shortfall(target, motion.lam)
        + STRAIN_WEIGHT * motion.strain
        + SLACK_WEIGHT * motion.output_slack
        + BRIDGE_WEIGHT * unbridged
        + discreteness * float(np.mean(4.0 * rho * (1.0 - rho)))
        + COUNT_WEIGHT * float(np.mean(rho))
    )


START_SPAN = 110.0
"""Half-width, in millimetres, of the box a start is drawn from.

Not the same as the box the search runs in, and the difference matters.  The
bounds are deliberately loose because the answer is not known in advance; a
*start* drawn that loosely is a different thing -- nodes scattered over half a
metre, elements three hundred millimetres long, and a piston that has to travel
74.  The first runs drew from the full box and spent their budget dragging that
back to scale.  This is a stroke and a half, which is the size of the mechanism
the specification implies without saying anything about its shape.
"""


def random_start(layout: Layout, rng: np.random.Generator) -> FloatArray:
    """A design drawn on the specification's own scale, every candidate half present.

    Starting the presences at one half rather than at random is deliberate: at
    :math:`p = 1` the stiffnesses enter almost linearly, so a symmetric start
    lets the first rung of the continuation decide which candidates matter
    instead of inheriting a decision from the draw.
    """
    lower, upper = default_bounds(layout)
    span = np.minimum(upper - lower, 2.0 * START_SPAN)
    middle = 0.5 * (np.clip(lower, -START_SPAN, None) + np.clip(upper, None, START_SPAN))
    x = middle + (rng.random(layout.size) - 0.5) * span
    for pin in layout.structure.pins:
        x[pin.slot] = lower[pin.slot] + rng.random() * min(
            60.0 - lower[pin.slot], upper[pin.slot] - lower[pin.slot]
        )
        x[pin.slot + 1] = -np.pi + rng.random() * 2.0 * np.pi
    for shaft in layout.structure.shafts:
        if shaft.angle_slot is not None:
            x[shaft.angle_slot] = -np.pi + rng.random() * 2.0 * np.pi
    for gear in layout.structure.gears:
        x[gear.phase_slot] = (rng.random() - 0.5) * 2.0 * START_SPAN * np.pi
    x[layout.presence_offset :] = 0.5
    return np.clip(x, lower, upper)


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


def _describe(layout: Layout, x: FloatArray, target: FloatArray, samples: int) -> Candidate:
    """Score a design at full penalty and name the candidates it kept."""
    motion = sweep(layout, x, samples=samples, penalty=PENALTY_SCHEDULE[-1])
    rho = layout.presences(x)
    names = layout.structure.names
    present = tuple(
        tuple(names[n] for n in element.nodes)
        for m, element in enumerate(layout.structure.elements)
        if rho[m] > PRESENCE_THRESHOLD
    )
    return Candidate(
        x=np.array(x, dtype=float),
        rms=motion_error(target, motion.lam),
        strain=motion.strain,
        discreteness=float(np.mean(4.0 * rho * (1.0 - rho))),
        present=present,
        gear_ratios=motion.gear_ratios,
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
    is a real mechanism's best fit rather than a blurred one's.

    Args:
        layout: The design vector's layout.
        x0: Starting design.
        target: Target motion, on any number of samples.
        samples: Input angles per sweep during the search.
        iterations: L-BFGS-B iterations per rung.
        schedule: Penalisation exponents.

    Returns:
        The discrete mechanism and its score.
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


@dataclass(frozen=True)
class Assessment:
    """What a synthesised mechanism actually delivers, against the requirements."""

    rms: float
    """Motion error against the target [mm]."""
    strain: float
    """Largest strain in a present candidate; a mechanism runs at the floor."""
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
    present: tuple[tuple[str, ...], ...]
    """Elements that survived."""
    gear_ratios: tuple[float, ...] = ()
    """Ratios of the gear pairs that survived."""
    output_slack: float = 0.0
    """How undetermined the piston is (:attr:`Motion.output_slack`)."""
    travel: float = 0.0
    """Peak-to-peak piston movement [mm]."""

    @property
    def delivers_the_stroke(self) -> bool:
        """Whether the piston moves far enough for the strokes to mean anything.

        A guard the results asked for.  A design whose piston wobbles by two
        millimetres can still show four monotone phases and an "asymmetry" above
        a millimetre, and one did -- start 3 of §5.7's re-run, whose harmonics
        are 1.9 and 2.5 mm against a target's 9.0 and 32.3.  Calling that
        extended expansion is an artefact of testing the *shape* of a motion
        without testing its *size*.
        """
        return self.travel >= 0.5 * DEFAULT_TARGETS.expansion_stroke

    @property
    def is_mechanism(self) -> bool:
        """Whether the input actually determines where the piston goes."""
        return self.output_slack < SLACK_TOLERANCE

    @property
    def asymmetry(self) -> float:
        """``STE - STC``: what makes the cycle extended-expansion rather than Otto."""
        return self.expansion_stroke - self.compression_stroke

    @property
    def is_extended_expansion(self) -> bool:
        """Whether this is a mechanism, and its two strokes actually differ.

        The mobility and travel tests come first on purpose.  An
        under-constrained answer can report any motion at all, extended
        expansion included, and so can a piston that barely moves.  Both
        happened.
        """
        return bool(
            self.is_mechanism
            and self.delivers_the_stroke
            and self.four_phases
            and self.asymmetry > 1.0
        )


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
        tuple(names[n] for n in element.nodes)
        for m, element in enumerate(layout.structure.elements)
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
        gear_ratios=motion.gear_ratios,
        output_slack=motion.output_slack,
        travel=float(np.ptp(motion.lam)),
    )


def format_assessment(assessment: Assessment) -> str:
    """A readable summary of a synthesised mechanism."""
    kept = " ".join("-".join(nodes) for nodes in assessment.present) or "none"
    gears = (
        " ".join(f"{r:g}:1" for r in assessment.gear_ratios)
        if assessment.gear_ratios
        else "none"
    )
    lines = [
        f"elements    {kept}",
        f"gears       {gears}",
        f"motion rms  {assessment.rms:8.3f} mm",
        f"peak strain {assessment.strain:8.2e}   (over-constraint)",
        f"slack       {assessment.output_slack:8.2e}   (under-constraint)",
        f"travel      {assessment.travel:8.3f} mm",
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
    if not assessment.is_mechanism:
        verdict = "the input does not determine the piston -- not a mechanism"
    elif not assessment.delivers_the_stroke:
        verdict = "the piston barely moves -- the strokes mean nothing"
    elif assessment.is_extended_expansion:
        verdict = "extended expansion"
    elif assessment.four_phases:
        verdict = "symmetric strokes -- an Otto engine"
    else:
        verdict = "not a four-stroke motion"
    return "\n".join([*lines, f"verdict     {verdict}"])


def exlink_in_the_domain(
    design: Design | None = None,
    samples: int = 720,
) -> tuple[GroundStructure, Layout, FloatArray]:
    """Write the studied mechanism into the synthesis domain, as a design vector.

    The check that says whether a negative result is about the *domain* or about
    the *search*, and here it says the search.  EXlink is three elements and a
    gear pair -- the swing rod ``P1-F1``, the trigonal link as a body
    ``P2-F1-F2``, the piston rod ``F2-S``, and a two-to-one pair -- so it is
    expressible in :func:`geared_ground_structure` exactly, and the spring model
    reproduces its analytic motion to about 2e-3 mm.

    Every corner is a design variable, which is what makes this possible: the
    two free nodes are the trigonal link's corners :math:`A` and :math:`E`, the
    pins carry a radius and a phase, the piston carries its axis, and each
    element's lengths follow from where those corners are put.  Nothing about
    the shape of the three-cornered body is fixed in advance.

    Args:
        design: The linkage to write in.  Defaults to the study's result.
        samples: Angles used to read its pose at zero input angle.

    Returns:
        The geared ground structure, its layout, and the design vector.
    """
    from .kinematics import solve
    from .reference import RELIABLE_DESIGN

    kinematics = solve(RELIABLE_DESIGN if design is None else design, samples=samples)
    origin = kinematics.R1[0]
    second = kinematics.R2[0] - origin
    pin_1 = kinematics.Q[0] - origin
    pin_2 = kinematics.D[0] - kinematics.R2[0]

    structure, layout = geared_ground_structure()
    x = np.zeros(layout.size, dtype=float)
    x[0:2] = second
    x[2:4] = [-250.0, -250.0]
    x[4], x[5] = float(np.hypot(*pin_1)), float(np.arctan2(pin_1[1], pin_1[0]))
    x[6], x[7] = float(np.hypot(*pin_2)), float(np.arctan2(pin_2[1], pin_2[0]))
    x[8:10] = kinematics.A[0] - origin
    x[10:12] = kinematics.E[0] - origin
    x[12:14] = kinematics.P[0] - origin

    names = [tuple(structure.names[n] for n in e.nodes) for e in structure.elements]
    rho = np.zeros(structure.n_presences, dtype=float)
    for element in (("P1", "F1"), ("P2", "F1", "F2"), ("F2", "S")):
        rho[names.index(element)] = 1.0
    pair = next(g for g, gear in enumerate(structure.gears) if gear.ratio == 2.0)
    rho[structure.n_members + pair] = 1.0
    x[layout.presence_offset :] = rho
    return structure, layout, x


def reachability(layout: Layout, x: FloatArray) -> FloatArray:
    """How well each driven pin is joined to the piston through the linkage.

    A bottleneck path: the value of a route is the weakest presence along it,
    and the value of a pin is the best route it has.  Fully present elements all
    the way give 1, no route at all gives 0, and a half-built chain gives its
    weakest link -- which is what makes it something a gradient can climb rather
    than a yes-or-no test.

    This is the path alone.  Whether a pin is worth reaching is a separate
    question -- :func:`driven_fraction` answers it -- and keeping the two apart
    matters in both directions.  Reaching a pin whose shaft nothing drives is
    vacuous, and a run duly reached 1.00 on both while keeping no gear at all;
    but *charging* for not reaching it would be wrong too, because a topology
    with no gear is a legitimate single-input linkage rather than a broken
    two-input one.

    Returns:
        One value per pin, in the order the structure lists its pins.
    """
    structure = layout.structure
    rho = layout.element_presences(x)
    weight = np.zeros((structure.n_nodes, structure.n_nodes), dtype=float)
    for m, element in enumerate(structure.elements):
        for i, j in element.springs:
            weight[i, j] = weight[j, i] = max(weight[i, j], float(rho[m]))

    slider = layout.slider_slot[0]
    out = np.zeros(len(structure.pins), dtype=float)
    for p, pin in enumerate(structure.pins):
        value = np.zeros(structure.n_nodes, dtype=float)
        value[pin.node] = 1.0
        for _sweep in range(structure.n_nodes):
            improved = np.max(np.minimum(value[:, None], weight), axis=0)
            better = np.maximum(value, improved)
            if np.allclose(better, value):
                break
            value = better
        out[p] = value[slider]
    return out


def driven_fraction(layout: Layout, x: FloatArray) -> FloatArray:
    """How far each pin's shaft is actually driven, from 0 to 1.

    The input shaft is driven by definition.  Any other is driven only as far as
    a gear presence says, because a shaft carrying no gear is a passive grounded
    pivot -- a perfectly good linkage element, and not a second *source* of
    motion.

    Returns:
        One value per pin, in the order the structure lists its pins.
    """
    structure = layout.structure
    gear_rho = layout.gear_presences(x)
    driven = np.ones(len(structure.shafts), dtype=float)
    for shaft_index, shaft in enumerate(structure.shafts):
        if not shaft.is_free:
            continue
        on_it = [
            float(gear_rho[g])
            for g, gear in enumerate(structure.gears)
            if shaft_index in (gear.first, gear.second)
        ]
        driven[shaft_index] = max(on_it) if on_it else 0.0
    return np.array([driven[pin.shaft] for pin in structure.pins], dtype=float)


@dataclass(frozen=True)
class Mechanism:
    """One discrete topology: which elements are present, and which gear."""

    elements: tuple[int, ...]
    """Indices into the structure's elements."""
    gear: int | None
    """Index of the gear pair fitted, or ``None``."""

    def presences(self, structure: GroundStructure) -> FloatArray:
        """The presence vector this topology stands for."""
        rho = np.zeros(structure.n_presences, dtype=float)
        for m in self.elements:
            rho[m] = 1.0
        if self.gear is not None:
            rho[structure.n_members + self.gear] = 1.0
        return rho

    def names(self, structure: GroundStructure) -> tuple[tuple[str, ...], ...]:
        """The elements it keeps, by node name."""
        return tuple(
            tuple(structure.names[n] for n in structure.elements[m].nodes)
            for m in self.elements
        )


def _joins_the_piston(structure: GroundStructure, elements: tuple[int, ...]) -> bool:
    """Whether both driven pins reach the piston through the chosen elements."""
    adjacency: dict[int, set[int]] = {i: set() for i in range(structure.n_nodes)}
    for m in elements:
        for i, j in structure.elements[m].springs:
            adjacency[i].add(j)
            adjacency[j].add(i)
    slider = next(i for i, kind in enumerate(structure.kinds) if kind == SLIDER)
    for pin in structure.pins:
        seen, stack = {pin.node}, [pin.node]
        while stack:
            for neighbour in adjacency[stack.pop()]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        if slider not in seen:
            return False
    return True


def is_determinate(
    layout: Layout,
    mechanism: Mechanism,
    rng: np.random.Generator,
    tries: int = 3,
) -> bool:
    """Whether the input angle fixes every unknown, at a general pose.

    The rank of the reduced Hessian is the test, taken at random geometries so
    that a topology is not rejected for being singular at one unlucky pose.  A
    topology failing it is under-constrained however its dimensions are chosen,
    which is the property §5.7 found reports whatever motion the solver's ridge
    happens to pick.

    The absent elements are given a stiffness of *exactly* zero rather than the
    usual floor.  The floor exists so that the continuation always has a Hessian
    to invert, and it leaves that Hessian nominally full rank whatever the
    topology, so a rank test taken through it answers a question about the
    tolerance rather than about the mechanism.
    """
    structure = layout.structure
    for _try in range(tries):
        x = random_start(layout, rng)
        x[layout.presence_offset :] = mechanism.presences(structure)
        k = layout.element_presences(x)[structure.spring_owner()]
        gear_k = layout.gear_presences(x)
        radii = gear_radii(layout, x)
        rest = member_lengths(layout, x)
        reduced = np.zeros(structure.n_reduced, dtype=float)
        mask = structure.free_mask()
        reduced[: int(mask.sum())] = initial_configuration(layout, x).reshape(-1)[mask]
        hessian = _assemble(layout, x, 0.0, reduced, k, rest, gear_k, radii)[2]
        if np.linalg.matrix_rank(hessian, tol=1.0e-6 * float(np.max(np.abs(hessian)))) == (
            structure.n_reduced
        ):
            return True
    return False


def enumerate_mechanisms(
    layout: Layout,
    seed: int = 0,
    max_elements: int = 4,
) -> tuple[Mechanism, ...]:
    """Every discrete topology in the domain that is a mechanism at all.

    The globalization this problem actually wants.  §5.7 leaves a search
    problem: the objective ranks the answer an order of magnitude ahead of
    anything six local runs found, and the runs commit to a short path in their
    first rung and never leave it.  The discrete part of the problem is small
    enough not to need a heuristic -- it can be *enumerated*, which is the same
    treatment §4.4 gives the gear lattice and the exhaustive baseline of
    Appendix C.3.

    Three conditions, applied in increasing cost.  The elements must supply
    exactly as many constraints as there are unknowns, counting a bar as one, a
    body as three and a fitted gear as one.  Both driven pins must reach the
    piston, which is §5.6's identity as a condition on the graph.  And the
    reduced Hessian must have full rank at a general pose, which rejects the
    topologies that are under-constrained however they are dimensioned.

    Args:
        layout: The design vector's layout.
        seed: Draws the random poses the rank test uses.
        max_elements: Largest number of elements to consider.

    Returns:
        The admissible topologies, smallest first.
    """
    from itertools import combinations

    structure = layout.structure
    cost = [len(element.springs) for element in structure.elements]
    unknowns = structure.n_reduced
    rng = np.random.default_rng(seed)

    found: list[Mechanism] = []
    for size in range(1, max_elements + 1):
        for elements in combinations(range(structure.n_members), size):
            supplied = sum(cost[m] for m in elements)
            if supplied not in (unknowns - 1, unknowns):
                continue
            if not _joins_the_piston(structure, elements):
                continue
            gears: list[int | None] = (
                list(range(len(structure.gears))) if supplied == unknowns - 1 else [None]
            )
            for gear in gears:
                candidate = Mechanism(elements=elements, gear=gear)
                if is_determinate(layout, candidate, rng):
                    found.append(candidate)
    return tuple(found)


def screen_mechanism(
    layout: Layout,
    mechanism: Mechanism,
    target: FloatArray,
    rng: np.random.Generator,
    draws: int = 4,
    iterations: int = 20,
    samples: int = 24,
) -> tuple[float, FloatArray]:
    """Give one topology a cheap chance, and score it.

    With the presences fixed the problem is no longer a topology optimization at
    all -- it is the well-conditioned continuous fit :mod:`exlink.synthesis`
    already solves reliably on a known linkage.  That is the point of enumerating:
    the hard, deceptive part of the search is the discrete part, and enumeration
    removes it, leaving each topology a problem that local descent handles.

    Args:
        layout: The design vector's layout.
        mechanism: The topology to try.
        target: Target piston motion.
        rng: Draws the starting geometries.
        draws: Random geometries to try before descending from the best.
        iterations: L-BFGS-B iterations allowed.
        samples: Input angles per sweep -- coarse, because this is a screen.

    Returns:
        The objective reached and the design vector reaching it.
    """
    presences = mechanism.presences(layout.structure)
    frozen = np.zeros(layout.size, dtype=bool)
    frozen[layout.presence_offset :] = True

    best_value, best_x = np.inf, None
    for _draw in range(draws):
        x = random_start(layout, rng)
        x[layout.presence_offset :] = presences
        value = objective(layout, x, target, samples, PENALTY_SCHEDULE[-1])
        if value < best_value:
            best_value, best_x = value, x
    assert best_x is not None

    polished = _minimise(
        layout, best_x, target, samples, PENALTY_SCHEDULE[-1], iterations, frozen=frozen
    )
    return objective(layout, polished, target, samples, PENALTY_SCHEDULE[-1]), polished
