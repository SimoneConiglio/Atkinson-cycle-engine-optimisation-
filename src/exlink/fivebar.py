r"""The geared five-bar: the alternative §5.9's synthesis found, on identical terms.

§5.9 reports that EXlink is not the only mechanism in the domain that meets the
specification.  The best thing the precision-point fit reached is **not** a
six-bar with a three-cornered link at all.  Once its dead corner is removed it
is a *floating pin held by two rods* -- one to a pin on the half-speed shaft,
one to a pin on the crankshaft -- driving the piston through a connecting rod:

.. code-block:: text

      R1 (half-speed shaft)          R2 (crankshaft, -2x)
        o--- q_1 ---> P1               o--- q_2 ---> P2
                       \                /
                    L_1 \              / L_2
                         \            /
                          `--- F ----'          the floating pin
                               |
                               | e              the connecting rod
                               |
                               S                the piston, on x = x_1

Five moving members against the EX-link's seven, and four journals against
seven, which is the whole reason it is worth pricing: a mechanism that meets the
same kinematic specification with fewer parts should win on friction and on
mass, and the question is whether it does once the loads are carried through.

A warning the kinematics alone already gives.  The fit that found this mechanism
was asked for a motion and nothing else, and it bought the motion with a 322 mm
connecting rod whose axis runs some 58 degrees from the cylinder -- a side load
of order 1.6 times the gas force, where a conventional engine keeps it below
0.3.  That is not an argument against the topology; it is an argument against
judging a mechanism by its motion, which is the argument this whole package
exists to make.  What the topology is worth is what it is worth *after*
:mod:`exlink.sizing`, :mod:`exlink.friction` and :mod:`exlink.mass_budget` have
had their say.

Conventions
-----------
Angles are measured from **+x** and counter-clockwise, which is the convention
:mod:`exlink.topology` solves in, so that the closed form here can be checked
directly against the spring model that found the mechanism.  This differs from
:mod:`exlink.kinematics`, which measures the EX-link's shaft angles from +y; the
two are never mixed.

:math:`\theta` is the angle of the **half-speed** shaft, one revolution per
four-stroke cycle.  The crankshaft turns at :math:`-2\theta`, the ratio the
2:1 pair imposes and the sense the mesh gives it, so the crank pin sits at
:math:`\theta_f - 2\theta`.  The half-speed crank's own phase is the datum and
is taken as zero: rotating it merely relabels :math:`\theta`, and every quantity
this module reports is read off the motion rather than off the datum.
"""

from __future__ import annotations

import math
from dataclasses import astuple, dataclass, fields
from typing import ClassVar

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]

#: Order of the design variables in the flat vector representation.
VARIABLE_NAMES: tuple[str, ...] = (
    "q_1",
    "q_2",
    "I",
    "theta_r",
    "theta_f",
    "L_1",
    "L_2",
    "e",
    "x_1",
)

#: Human-readable description of each design variable.
VARIABLE_DESCRIPTIONS: dict[str, str] = {
    "q_1": "crank radius on the half-speed shaft [mm]",
    "q_2": "crank radius on the crankshaft [mm]",
    "I": "distance between the two shaft axes [mm]",
    "theta_r": "orientation of the half-speed-shaft-to-crankshaft axis [deg]",
    "theta_f": "crank dephasing: crankshaft pin angle at theta = 0 [deg]",
    "L_1": "rod from the half-speed crank pin to the floating pin [mm]",
    "L_2": "rod from the crank pin to the floating pin [mm]",
    "e": "connecting rod from the floating pin to the piston [mm]",
    "x_1": "abscissa of the cylinder axis [mm]",
}

#: The design variables stored in degrees.
ANGULAR_VARIABLES: frozenset[str] = frozenset({"theta_r", "theta_f"})

#: Speed of the crankshaft relative to the half-speed shaft.
CRANK_RATIO: float = -2.0
"""Set by the 2:1 pair, sign included.

The mesh makes the two shafts turn in opposite senses, and the ratio has to be a
whole number or the cycle never closes -- a shaft turning :math:`r` times per
input revolution returns to its start only when :math:`r` is an integer.  §5.9
finds this value rather than assuming it: of the seven best mechanisms in the
family, every one took 2:1 from a catalogue of 1, 2, 3 and 4 to one.
"""


@dataclass(frozen=True)
class FiveBar:
    """One set of geared five-bar dimensions.

    Angles are stored in **degrees**, as in :class:`exlink.design.Design`,
    because that is how they are bounded and reported; consumers convert through
    :attr:`theta_r_rad` and :attr:`theta_f_rad`.

    The two branch flags are not dimensions but *assembly choices*: a pair of
    circles meets in two points and a rod reaches a cylinder axis at two
    heights.  Both are discrete properties of the machine as built, fixed once
    and never crossed, because crossing one is a different mechanism rather than
    a different design.
    """

    q_1: float
    q_2: float
    I: float
    theta_r: float
    theta_f: float
    L_1: float
    L_2: float
    e: float
    x_1: float
    branch: int = 1
    """Which intersection of the two rod circles carries the floating pin."""
    piston_branch: int = -1
    """Which side of the floating pin the piston sits on."""

    names: ClassVar[tuple[str, ...]] = VARIABLE_NAMES

    @classmethod
    def from_array(cls, values: ArrayLike, branch: int = 1, piston_branch: int = -1) -> FiveBar:
        """Build a design from a flat vector ordered as :data:`VARIABLE_NAMES`."""
        array = np.asarray(values, dtype=float).ravel()
        if array.size != len(VARIABLE_NAMES):
            msg = f"expected {len(VARIABLE_NAMES)} design variables, got {array.size}"
            raise ValueError(msg)
        return cls(*array.tolist(), branch=branch, piston_branch=piston_branch)

    @classmethod
    def from_mapping(cls, mapping: dict[str, ArrayLike]) -> FiveBar:
        """Build a design from a ``{name: value}`` mapping (GEMSEO input data)."""
        return cls(*[float(np.ravel(mapping[name])[0]) for name in VARIABLE_NAMES])

    def to_array(self) -> FloatArray:
        """Return the dimensions as a flat vector ordered as :data:`VARIABLE_NAMES`."""
        return np.array(astuple(self)[: len(VARIABLE_NAMES)], dtype=float)

    def to_mapping(self) -> dict[str, FloatArray]:
        """Return the design as GEMSEO-style ``{name: array([value])}`` data."""
        return {name: np.array([getattr(self, name)], dtype=float) for name in VARIABLE_NAMES}

    def replace(self, **changes: float) -> FiveBar:
        """Return a copy with some variables overridden."""
        current = {field.name: getattr(self, field.name) for field in fields(self)}
        unknown = set(changes) - set(current)
        if unknown:
            msg = f"unknown design variables: {sorted(unknown)}"
            raise ValueError(msg)
        current.update(changes)
        return FiveBar(**current)

    @property
    def theta_r_rad(self) -> float:
        """Shaft-axis orientation in radians."""
        return math.radians(self.theta_r)

    @property
    def theta_f_rad(self) -> float:
        """Crank dephasing in radians."""
        return math.radians(self.theta_f)

    @property
    def crank_centre(self) -> FloatArray:
        """Where the crankshaft axis sits, the half-speed axis being the origin."""
        return self.I * np.array(
            [math.cos(self.theta_r_rad), math.sin(self.theta_r_rad)], dtype=float
        )


class AssemblyError(ValueError):
    """The dimensions do not close at some crank angle.

    Raised rather than returned as ``nan`` because a five-bar that comes apart
    partway round is not a worse design, it is not a design: no continuation of
    the motion past the break has any meaning, and a cycle read from it would be
    read from a fiction.
    """


@dataclass(frozen=True)
class FiveBarKinematics:
    """Where every joint goes, over one revolution of the half-speed shaft."""

    theta: FloatArray
    """Half-speed shaft angle [rad], uniformly spaced over ``[0, 2 pi)``."""
    P1: FloatArray
    """Pin on the half-speed crank, ``(n, 2)`` [mm]."""
    P2: FloatArray
    """Pin on the crankshaft, ``(n, 2)`` [mm]."""
    F: FloatArray
    """The floating pin, ``(n, 2)`` [mm]."""
    S: FloatArray
    """The piston, ``(n, 2)`` [mm]."""
    lam: FloatArray
    """Piston height ``lambda(theta)`` [mm] -- the ordinate of ``S``."""

    @property
    def rod_angle(self) -> FloatArray:
        """Connecting-rod inclination from the cylinder axis [rad].

        The quantity the side load is set by: a piston pushed along a rod at
        angle :math:`\\psi` presses on the liner with :math:`\\tan\\psi` times
        the gas force.
        """
        offset = self.S[:, 0] - self.F[:, 0]
        along = self.S[:, 1] - self.F[:, 1]
        return np.arctan2(offset, along)

    @property
    def side_load_ratio(self) -> float:
        """Largest ``tan`` of the rod angle over the cycle, dimensionless.

        A conventional engine keeps this below about 0.3.  It is reported here
        because the synthesis that produced the first five-bar was asked for a
        motion and bought it with a rod running 58 degrees off the cylinder.
        """
        return float(np.max(np.abs(np.tan(self.rod_angle))))


def solve(
    design: FiveBar, samples: int = 720, theta: ArrayLike | None = None
) -> FiveBarKinematics:
    """Close the linkage at every crank angle, in closed form.

    Two circle intersections in sequence, which is all a five-bar with a
    prescribed shaft ratio needs.  Both crank pins are known outright from the
    shaft angles, so the floating pin is where two rods of known length reach
    from two known points, and the piston is where the connecting rod reaches
    the cylinder axis.

    Args:
        design: The dimensions.
        samples: Angles over one revolution of the half-speed shaft.
        theta: Evaluate at these angles instead, in radians.  Used to check the
            closed form against a solver on the solver's own datum, without
            rounding the offset between them to the nearest sample.

    Returns:
        Every joint position and the piston motion.

    Raises:
        AssemblyError: If the rods cannot close, or the connecting rod cannot
            reach the cylinder axis, at any angle.
    """
    theta = (
        np.linspace(0.0, 2.0 * math.pi, samples, endpoint=False)
        if theta is None
        else np.asarray(theta, dtype=float).ravel()
    )
    samples = theta.size
    crank = design.theta_f_rad + CRANK_RATIO * theta

    p1 = design.q_1 * np.stack([np.cos(theta), np.sin(theta)], axis=1)
    p2 = design.crank_centre + design.q_2 * np.stack([np.cos(crank), np.sin(crank)], axis=1)

    span = p2 - p1
    distance = np.hypot(span[:, 0], span[:, 1])
    if np.any(distance <= 0.0):
        msg = "the two crank pins coincide at some angle"
        raise AssemblyError(msg)

    along = (design.L_1**2 - design.L_2**2 + distance**2) / (2.0 * distance)
    height_squared = design.L_1**2 - along**2
    if np.any(height_squared < 0.0):
        worst = float(np.min(height_squared))
        msg = f"the two rods do not reach each other (worst gap {worst:.3f} mm^2)"
        raise AssemblyError(msg)
    height = np.sqrt(height_squared)

    unit = span / distance[:, None]
    normal = np.stack([-unit[:, 1], unit[:, 0]], axis=1)
    floating = p1 + along[:, None] * unit + design.branch * height[:, None] * normal

    offset = design.x_1 - floating[:, 0]
    reach_squared = design.e**2 - offset**2
    if np.any(reach_squared < 0.0):
        worst = float(np.min(reach_squared))
        msg = f"the connecting rod does not reach the cylinder axis (worst {worst:.3f} mm^2)"
        raise AssemblyError(msg)
    piston_height = floating[:, 1] + design.piston_branch * np.sqrt(reach_squared)
    piston = np.stack([np.full(samples, design.x_1), piston_height], axis=1)

    return FiveBarKinematics(theta=theta, P1=p1, P2=p2, F=floating, S=piston, lam=piston_height)


def from_topology(x: ArrayLike, layout: object | None = None) -> tuple[FiveBar, float]:
    """Read a five-bar out of a §5.9 synthesis result.

    The synthesis works in a ground structure whose coordinates are node
    positions and pin radii, not link lengths; this turns one of its design
    vectors into dimensions.  The element lengths are distances in the pose at
    zero input angle, which is exactly how :mod:`exlink.topology` defines them,
    so nothing is refitted here -- the mechanism is only re-expressed.

    Args:
        x: A design vector on :func:`exlink.topology.geared_ground_structure`,
            carrying the five-bar's elements.
        layout: The layout it is drawn on; the geared one by default.

    Returns:
        The dimensions, and the datum offset in radians -- the phase of the
        half-speed crank pin, which this module takes as zero and the ground
        structure does not.
    """
    from .topology import geared_ground_structure

    values = np.asarray(x, dtype=float).ravel()
    if layout is None:
        _structure, layout = geared_ground_structure()

    q_1, phi_1 = float(values[4]), float(values[5])
    q_2, phi_2 = float(values[6]), float(values[7])
    centre = values[0:2]
    span = float(np.hypot(*centre))

    # The crankshaft's angle at zero input is *not* the value seeded in the
    # design vector.  That slot is only a starting guess for a free coordinate,
    # and the seeded pose is not an equilibrium: the rods are all at their rest
    # lengths there by construction, but nothing has made the mesh close, so the
    # solver turns the geared shaft until it does.  Reading the seed instead of
    # the mesh puts the whole linkage on the wrong branch -- 17.6 mm of piston
    # error, which is how this was found.
    presences = np.asarray(values[layout.presence_offset :], dtype=float)  # type: ignore[attr-defined]
    structure = layout.structure  # type: ignore[attr-defined]
    fitted = [
        g for g, _gear in enumerate(structure.gears) if presences[structure.n_members + g] > 0.5
    ]
    if len(fitted) != 1:
        msg = f"expected exactly one fitted gear pair, found {len(fitted)}"
        raise ValueError(msg)
    gear = structure.gears[fitted[0]]
    ratio = float(gear.ratio)
    phase = float(values[gear.phase_slot])
    crank_angle = phase * (1.0 + ratio) / span

    # Two different poses are in play and they must not be mixed.  The link
    # *lengths* are distances in the seeded pose, because that is how
    # exlink.topology.member_lengths defines a rest length -- so the crank pin
    # is taken at the seeded shaft angle here.  The crank *angle at zero input*
    # is the one the mesh imposes, computed above.  Using the mesh angle for the
    # lengths too stretches L_2 from 99.16 mm to 135.71 mm and the linkage then
    # will not close at all, which is how this was found.
    seeded = float(values[14])
    floating = values[10:12]
    piston = np.array([float(values[12]), float(values[13])], dtype=float)
    pin_1 = q_1 * np.array([math.cos(phi_1), math.sin(phi_1)])
    pin_2 = centre + q_2 * np.array([math.cos(seeded + phi_2), math.sin(seeded + phi_2)])

    # The ground structure puts the half-speed pin at phi_1 when its input angle
    # is zero; this module puts it at zero.  That is a relabelling of the *angle*
    # and not a rotation of the frame, so the shaft axis and the cylinder keep
    # their bearings and only the crank dephasing moves -- by twice phi_1,
    # because the crankshaft turns twice for every turn of the half-speed shaft.
    design = FiveBar(
        q_1=q_1,
        q_2=q_2,
        I=float(np.hypot(*centre)),
        theta_r=math.degrees(math.atan2(centre[1], centre[0])),
        theta_f=math.degrees(crank_angle + phi_2 + ratio * phi_1),
        L_1=float(np.hypot(*(floating - pin_1))),
        L_2=float(np.hypot(*(floating - pin_2))),
        e=float(np.hypot(*(floating - piston))),
        x_1=float(values[12]),
    )
    return design, phi_1


def in_the_domain(design: FiveBar, samples: int = 720) -> tuple[object, object, FloatArray]:
    """Write a five-bar into the synthesis domain, as a design vector.

    The inverse of :func:`from_topology`, and what makes the closed form
    checkable without storing a design vector: the mechanism goes into
    :func:`exlink.topology.geared_ground_structure`, the spring model solves it
    from first principles, and the two motions are compared.  It is the same
    check :func:`exlink.topology.exlink_in_the_domain` provides for the EX-link.

    The datum is chosen so that the half-speed crank pin sits at angle zero, and
    the gear phase is set to close the mesh at that pose, so the seeded
    configuration *is* the equilibrium and the two models share a datum.  Both
    matter: a seed that does not satisfy the mesh is not an equilibrium, and the
    solver will move to a different branch of the assembly.

    Args:
        design: The dimensions.
        samples: Angles used to read the pose at zero input angle.

    Returns:
        The ground structure, its layout, and the design vector.

    Raises:
        AssemblyError: If the linkage does not close.
    """
    from .topology import geared_ground_structure

    structure, layout = geared_ground_structure()
    pose = solve(design, theta=np.zeros(1))

    x = np.zeros(layout.size, dtype=float)
    x[0:2] = design.crank_centre
    x[2:4] = [-250.0, -250.0]
    x[4], x[5] = design.q_1, 0.0
    x[6], x[7] = design.q_2, 0.0
    floating = pose.F[0]
    piston = pose.S[0]
    # The dead corner has to make a non-degenerate triangle with the other two
    # and is otherwise free -- that is what "dead" means.
    midpoint = 0.5 * (floating + piston)
    along = piston - floating
    x[8:10] = midpoint + 0.25 * np.array([-along[1], along[0]])
    x[10:12] = floating
    x[12], x[13] = design.x_1, float(piston[1])
    x[14] = design.theta_f_rad

    ratio = -CRANK_RATIO
    pair = next(g for g, gear in enumerate(structure.gears) if gear.ratio == ratio)
    # Close the mesh at the seeded pose: with the input shaft at zero the
    # residual is the second radius times the crankshaft angle, less the phase.
    x[structure.gears[pair].phase_slot] = design.theta_f_rad * design.I / (1.0 + ratio)

    names = [tuple(structure.names[n] for n in e.nodes) for e in structure.elements]
    rho = np.zeros(structure.n_presences, dtype=float)
    for element in (("P1", "F2"), ("P2", "F2"), ("F1", "F2", "S")):
        rho[names.index(element)] = 1.0
    rho[structure.n_members + pair] = 1.0
    x[layout.presence_offset :] = rho
    return structure, layout, x


SYNTHESISED = FiveBar(
    q_1=14.220180997332037,
    q_2=25.75582395872892,
    I=77.21651849347572,
    theta_r=69.51465926467662,
    theta_f=8.152678478806534,
    L_1=90.46002498640651,
    L_2=99.16459687028852,
    e=322.07325109283,
    x_1=-175.96192099022505,
    branch=-1,
    piston_branch=1,
)
"""The five-bar exactly as §5.9's precision-point fit produced it.

Kept as it came out, unimproved, because it is the *input* to the
multidisciplinary study and not its result.  It meets the kinematic
specification -- 0.484 mm of requirement error, an expansion stroke of 74.12 mm
against 74.000, a compression ratio of 15.98 against 16 -- and it is mechanically
indefensible: its connecting rod is 322 mm long and runs up to 60.6 degrees from
the cylinder axis, so the piston presses on the liner with 1.77 times the gas
force where a conventional engine stays below 0.3.

That contrast is the point of carrying it further.  A synthesis asked for a
motion will buy the motion with whatever it is not charged for, and what it is
not charged for here is the entire mechanical design.  §5.10 puts a price on it.
"""


#: The five moving links, and the design variable each takes its length from.
MEMBERS: tuple[tuple[str, str, str, str], ...] = (
    ("crank_1", "q_1", "O1", "P1"),
    ("rod_1", "L_1", "P1", "F"),
    ("crank_2", "q_2", "O2", "P2"),
    ("rod_2", "L_2", "P2", "F"),
    ("con_rod", "e", "F", "S"),
)
"""``(name, length attribute, start joint, end joint)`` for each link.

Every body is one member here, where the EX-link's trigonal link is three.  That
is the whole structural difference between the two mechanisms, and it is why this
one is worth pricing: five moving members against seven, four journals against
seven, and the same gear pair.
"""

#: Order of the unknowns in the equilibrium system.
_UNKNOWNS: tuple[str, ...] = (
    "O1x", "O1y", "P1x", "P1y", "F1x", "F1y",
    "O2x", "O2y", "P2x", "P2y", "F2x", "F2y",
    "Sx", "Sy", "N", "M_liner", "W", "T",
)  # fmt: skip
_INDEX = {name: i for i, name in enumerate(_UNKNOWNS)}
_N_UNKNOWNS = len(_UNKNOWNS)


@dataclass(frozen=True)
class FiveBarMassProperties:
    """Mass, centre of mass, inertia and orientation of every body."""

    member_mass: dict[str, float]
    """[tonne]"""
    body_com: dict[str, FloatArray]
    """Centre of mass over the revolution, ``(n, 2)`` [mm]."""
    body_inertia: dict[str, float]
    """About the body's own centre of mass [tonne.mm^2]."""
    body_angle: dict[str, FloatArray]
    """Orientation [rad], for the angular acceleration."""


def mass_properties(
    kinematics: FiveBarKinematics,
    design: FiveBar,
    diameters: dict[str, float],
    density: float,
    piston_mass: float,
    piston_length: float,
    bore_ratio: float = 0.0,
) -> FiveBarMassProperties:
    """Mass properties of the five links and the piston.

    Identical in form to :func:`exlink.dynamics.mass_properties`: a link is a
    uniform bar of its own diameter and length, its centre of mass at its
    midpoint, its inertia that of a cylinder about a transverse axis through its
    centroid, and the piston translates so its rotational inertia never enters.
    Using the same formulae is what makes the mass comparison a comparison of
    mechanisms rather than of modelling conventions.
    """
    joints = {
        "O1": np.zeros_like(kinematics.P1),
        "O2": np.zeros_like(kinematics.P1) + design.crank_centre,
        "P1": kinematics.P1,
        "P2": kinematics.P2,
        "F": kinematics.F,
        "S": kinematics.S,
    }

    member_mass: dict[str, float] = {}
    body_com: dict[str, FloatArray] = {}
    body_inertia: dict[str, float] = {}
    for name, attribute, start, end in MEMBERS:
        length = abs(float(getattr(design, attribute)))
        diameter = float(diameters[name])
        area = math.pi * diameter**2 / 4.0 * (1.0 - bore_ratio**2)
        mass = density * area * length
        member_mass[name] = mass
        body_com[name] = 0.5 * (joints[start] + joints[end])
        body_inertia[name] = (
            mass * (0.75 * diameter**2 * (1.0 + bore_ratio**2) + length**2) / 12.0
        )

    body_com["piston"] = kinematics.S + np.array([0.0, 0.5 * piston_length])
    body_inertia["piston"] = 0.0

    zeros = np.zeros_like(kinematics.theta)
    crank = design.theta_f_rad + CRANK_RATIO * kinematics.theta
    body_angle = {
        "crank_1": kinematics.theta,
        "rod_1": np.arctan2(
            kinematics.F[:, 1] - kinematics.P1[:, 1], kinematics.F[:, 0] - kinematics.P1[:, 0]
        ),
        "crank_2": crank,
        "rod_2": np.arctan2(
            kinematics.F[:, 1] - kinematics.P2[:, 1], kinematics.F[:, 0] - kinematics.P2[:, 0]
        ),
        "con_rod": np.arctan2(
            kinematics.S[:, 1] - kinematics.F[:, 1], kinematics.S[:, 0] - kinematics.F[:, 0]
        ),
        "piston": zeros,
    }
    masses = dict(member_mass)
    masses["piston"] = piston_mass
    return FiveBarMassProperties(
        member_mass=masses,
        body_com=body_com,
        body_inertia=body_inertia,
        body_angle=body_angle,
    )


@dataclass(frozen=True)
class FiveBarLoads:
    """Joint reactions, gear load and output torque over one cycle."""

    kinematics: FiveBarKinematics
    speed: float
    """Crankshaft speed [rad/s]."""
    properties: FiveBarMassProperties
    reaction: dict[str, FloatArray]
    """Force transmitted forward along the chain at each joint, ``(n, 2)`` [N].

    Keys ``O1``, ``P1``, ``F1`` follow the half-speed branch; ``O2``, ``P2``,
    ``F2`` the crankshaft branch; ``S`` the wrist pin.
    """
    liner_force: FloatArray
    """Side force from the cylinder liner [N]."""
    liner_moment: FloatArray
    """Reaction moment on the piston from the guide [N.mm]."""
    gear_force: FloatArray
    """Tooth load along the line of action [N]."""
    torque: FloatArray
    """Output torque on the crankshaft [N.mm]."""
    gas_force: FloatArray
    """Applied gas force on the crown [N], as solved with."""
    conditioning: float
    """Worst condition number of the equilibrium matrix."""

    @property
    def peak_bearing_load(self) -> float:
        """Largest load on either main journal [N]."""
        return float(
            max(
                np.max(np.linalg.norm(self.reaction["O1"], axis=1)),
                np.max(np.linalg.norm(self.reaction["O2"], axis=1)),
            )
        )

    @property
    def indicated_work(self) -> float:
        """Work the gas does on the piston over the cycle [N.mm].

        Read off the gas force and the piston motion, independently of the
        torque, so that the two can be compared -- which is the check that says
        whether the equilibrium was assembled correctly.
        """
        step = 2.0 * math.pi / self.kinematics.theta.size
        velocity = np.gradient(self.kinematics.lam, step, edge_order=2)
        return float(-np.sum(self.gas_force * velocity) * step)

    @property
    def shaft_work(self) -> float:
        """Work delivered at the crankshaft over the cycle [N.mm].

        The crankshaft turns :data:`CRANK_RATIO` times for each turn of the
        half-speed shaft, so its angle advances by that factor per step.
        """
        step = 2.0 * math.pi / self.kinematics.theta.size
        return float(np.sum(self.torque) * abs(CRANK_RATIO) * step)


def solve_loads(
    kinematics: FiveBarKinematics,
    design: FiveBar,
    gas_force: FloatArray,
    properties: FiveBarMassProperties,
    speed: float,
    pressure_angle: float,
    piston_length: float,
) -> FiveBarLoads:
    """Equilibrium of all six bodies at once, inertia included.

    Eighteen unknowns and eighteen equations: seven joint force pairs, the liner
    normal force and its reaction moment, the gear tooth load and the output
    torque.  Assembled exactly as :func:`exlink.dynamics.solve` assembles the
    EX-link's, and solved simultaneously for the same reason -- with mass in the
    rods nothing is a two-force member and sequential elimination does not close.

    The floating pin is carried by the connecting rod, so both rods push on that
    body and no separate massless-pin equation is needed.

    Args:
        kinematics: A solved five-bar.
        design: Its dimensions.
        gas_force: Gas force on the crown at each angle [N], positive downwards.
        properties: Masses, centres of mass and inertias.
        speed: Crankshaft speed [rad/s].  Zero recovers the quasi-static result.
        pressure_angle: Gear pressure angle [rad].
        piston_length: Piston length [mm], for the crown offset.

    Returns:
        Reactions, tooth load and torque over one revolution of the half-speed
        shaft, which is one four-stroke cycle.
    """
    from .derivatives import ramp_derivative, spectral_derivative

    n = kinematics.theta.size
    half_speed = speed / abs(CRANK_RATIO)

    def second(points: FloatArray) -> FloatArray:
        return np.stack(
            [spectral_derivative(points[:, 0], 2), spectral_derivative(points[:, 1], 2)],
            axis=-1,
        )

    # Derivatives are taken with respect to the half-speed shaft's angle, since
    # that is the variable the trajectories are sampled in, so the scale is its
    # speed and not the crankshaft's.
    scale = half_speed**2
    body_acceleration = {body: scale * second(com) for body, com in properties.body_com.items()}
    body_angular_acceleration = {
        body: scale * ramp_derivative(angle, 2) for body, angle in properties.body_angle.items()
    }
    # Both shafts turn at a constant rate; force the exact zero rather than
    # leaving spectral round-off in the inertia couples.
    for body in ("crank_1", "crank_2", "piston"):
        body_angular_acceleration[body] = np.zeros(n)

    joints = {
        "O1": np.zeros((n, 2)),
        "O2": np.zeros((n, 2)) + design.crank_centre,
        "P1": kinematics.P1,
        "P2": kinematics.P2,
        "F": kinematics.F,
        "S": kinematics.S,
    }
    crown = np.stack([np.full(n, design.x_1), kinematics.lam], axis=-1)

    # -- gear mesh geometry, on the EX-link's convention --------------------------
    axis = np.array([math.cos(design.theta_r_rad), math.sin(design.theta_r_rad)])
    line_of_action = np.array(
        [
            math.cos(design.theta_r_rad - math.pi / 2.0 + pressure_angle),
            math.sin(design.theta_r_rad - math.pi / 2.0 + pressure_angle),
        ]
    )
    radius_1 = design.I * abs(CRANK_RATIO) / (1.0 + abs(CRANK_RATIO))
    radius_2 = design.I / (1.0 + abs(CRANK_RATIO))
    contact_1 = np.zeros((n, 2)) + radius_1 * axis
    contact_2 = joints["O2"] - radius_2 * axis

    matrix = np.zeros((n, _N_UNKNOWNS, _N_UNKNOWNS))
    rhs = np.zeros((n, _N_UNKNOWNS))

    def add_joint(row: int, body: str, joint: str, key: str, sign: float) -> None:
        """Enter an unknown joint force into a body's three equations."""
        ix, iy = _INDEX[key + "x"], _INDEX[key + "y"]
        matrix[:, row, ix] += sign
        matrix[:, row + 1, iy] += sign
        arm = joints[joint] - properties.body_com[body]
        matrix[:, row + 2, ix] += -sign * arm[:, 1]
        matrix[:, row + 2, iy] += sign * arm[:, 0]

    def set_inertia(row: int, body: str) -> None:
        mass = properties.member_mass[body]
        rhs[:, row] += mass * body_acceleration[body][:, 0]
        rhs[:, row + 1] += mass * body_acceleration[body][:, 1]
        rhs[:, row + 2] += properties.body_inertia[body] * body_angular_acceleration[body]

    def add_gear(row: int, body: str, contact: FloatArray, sign: float) -> None:
        matrix[:, row, _INDEX["W"]] += sign * line_of_action[0]
        matrix[:, row + 1, _INDEX["W"]] += sign * line_of_action[1]
        arm = contact - properties.body_com[body]
        matrix[:, row + 2, _INDEX["W"]] += sign * (
            arm[:, 0] * line_of_action[1] - arm[:, 1] * line_of_action[0]
        )

    # -- the half-speed shaft's crank throw ---------------------------------------
    add_joint(0, "crank_1", "O1", "O1", +1.0)
    add_joint(0, "crank_1", "P1", "P1", -1.0)
    add_gear(0, "crank_1", contact_1, -1.0)
    set_inertia(0, "crank_1")

    # -- the rod from that crank to the floating pin ------------------------------
    add_joint(3, "rod_1", "P1", "P1", +1.0)
    add_joint(3, "rod_1", "F", "F1", -1.0)
    set_inertia(3, "rod_1")

    # -- the crankshaft's crank throw, which the work is taken from ---------------
    add_joint(6, "crank_2", "O2", "O2", +1.0)
    add_joint(6, "crank_2", "P2", "P2", -1.0)
    add_gear(6, "crank_2", contact_2, +1.0)
    # Signed so that a positive unknown is torque delivered *out* of the shaft.
    # The check that fixes the sign is that its integral must equal the p-V loop
    # area, not a convention -- and it did fix it: the opposite sign balanced the
    # power to one part in twelve thousand with the ratio at minus one.
    matrix[:, 8, _INDEX["T"]] += +1.0
    set_inertia(6, "crank_2")

    # -- the rod from the crankshaft to the floating pin --------------------------
    add_joint(9, "rod_2", "P2", "P2", +1.0)
    add_joint(9, "rod_2", "F", "F2", -1.0)
    set_inertia(9, "rod_2")

    # -- the connecting rod, which carries the floating pin -----------------------
    add_joint(12, "con_rod", "F", "F1", +1.0)
    add_joint(12, "con_rod", "F", "F2", +1.0)
    add_joint(12, "con_rod", "S", "S", -1.0)
    set_inertia(12, "con_rod")

    # -- the piston ---------------------------------------------------------------
    add_joint(15, "piston", "S", "S", +1.0)
    matrix[:, 15, _INDEX["N"]] += 1.0
    matrix[:, 17, _INDEX["M_liner"]] += 1.0
    set_inertia(15, "piston")
    arm = crown - properties.body_com["piston"]
    # The gas resultant is (0, -gas_force) on the crown, so it crosses to the
    # right-hand side as +gas_force.
    rhs[:, 16] += gas_force
    rhs[:, 17] += arm[:, 0] * gas_force

    solution: FloatArray = np.linalg.solve(matrix, rhs[..., None])[..., 0].astype(
        np.float64, copy=False
    )
    reaction = {
        "O1": solution[:, 0:2],
        "P1": solution[:, 2:4],
        "F1": solution[:, 4:6],
        "O2": solution[:, 6:8],
        "P2": solution[:, 8:10],
        "F2": solution[:, 10:12],
        "S": solution[:, 12:14],
    }
    return FiveBarLoads(
        kinematics=kinematics,
        speed=speed,
        properties=properties,
        reaction=reaction,
        liner_force=solution[:, _INDEX["N"]],
        liner_moment=solution[:, _INDEX["M_liner"]],
        gear_force=solution[:, _INDEX["W"]],
        torque=solution[:, _INDEX["T"]],
        gas_force=np.asarray(gas_force, dtype=float),
        conditioning=float(np.max(np.linalg.cond(matrix))),
    )
