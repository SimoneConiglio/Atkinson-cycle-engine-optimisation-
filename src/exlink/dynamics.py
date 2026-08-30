"""Dynamic load analysis: the quasi-static chain with inertia restored.

The 2015 report stops short of this deliberately -- "to have the masses of the
pieces we have to know their shape so we should have a first design, so those
passages are for another iteration".  This module is that iteration.

Why the report's solution method cannot simply be extended
----------------------------------------------------------
Without inertia every rod is a *two-force member*: the forces at its two ends
are equal, opposite and collinear with the rod.  That is what lets the report
eliminate unknowns one body at a time, from the piston down to the crankshaft.

Add inertia and that collapses.  A rod with mass has a distributed d'Alembert
load along it, so its end forces are neither collinear nor equal, and no body
can be solved before its neighbours.  The whole mechanism has to be solved
simultaneously.

Assembling the whole system
---------------------------
Counting with Gruebler over 7 links (6 moving plus ground), 8 lower pairs
(7 revolutes and the piston's prismatic guide) and 1 higher pair (the gear
mesh)::

    M = 3(7 - 1) - 2(8) - 1 = 1

One degree of freedom, so the load problem is statically determinate: exactly
18 scalar unknowns against 6 bodies x 3 equilibrium equations.  Assembling and
solving that 18x18 system at every crank angle is what :func:`solve` does.

The determinant of that matrix is the same quantity the report's condition (4a)
protects: at a critical configuration the mechanism gains a degree of freedom,
the matrix goes singular, and the internal forces blow up.  The conditioning is
reported in :attr:`DynamicLoads.conditioning` so the connection is visible.

Constant crankshaft speed
-------------------------
The report analyses the mechanism at constant ``Omega``, which is inherited
here: ``d/dt = Omega d/dtheta_1`` exactly, with no angular acceleration of the
crankshaft.  Two consequences worth stating, because they simplify the
bookkeeping and are easy to get wrong:

* Both shaft assemblies turn at constant rate (``theta_2 = -2 theta_1 +
  theta_f``), so neither has an angular acceleration and neither contributes an
  inertia couple.
* The shafts and gears are concentric with their own axes, so their centres of
  mass do not move and they exert no inertia force at all.  Only the *offset*
  crank throws do, which is why the bodies below are the crank arms rather than
  whole shaft assemblies.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec
from .derivatives import ramp_derivative, spectral_derivative
from .kinematics import Kinematics

FloatArray = NDArray[np.float64]

DEFAULT_SPEED_RPM = 1500.0
"""Crankshaft speed used by default [rev/min].

Chosen because it is where the problem is interesting.  A Shell Eco-marathon
engine runs slowly by design, and for this linkage the structural mass grows
roughly as the *cube* of the acceleration level (see :mod:`exlink.coupled`):
the published geometry needs 0.25 kg of links at rest, 1.0 kg at 1000 rpm,
8.4 kg at 1500 rpm, and cannot be built at all much above 2000 rpm.  1500 rpm
therefore exercises the coupling hard while leaving sensible designs reachable.
"""


def rpm_to_rad_per_s(rpm: float) -> float:
    """Convert rev/min to rad/s."""
    return float(rpm) * 2.0 * np.pi / 60.0


@dataclass(frozen=True)
class Member:
    """One sized structural member, spanning two joints of the linkage.

    Attributes:
        name: Identifier, also the sizing output key.
        start: Joint at the member's first end.
        end: Joint at its second end.
        kind: ``"link"`` for a free two-node member (swing rod, piston rod),
            ``"cantilever"`` for a crank throw built into its shaft, or
            ``"truss"`` for a side of the trigonal link.
        length_attribute: Name of the :class:`~exlink.design.Design` property
            giving its length.
    """

    name: str
    start: str
    end: str
    kind: str
    length_attribute: str


MEMBERS: tuple[Member, ...] = (
    Member("crank_1", "Q", "R1", "cantilever", "q_1"),
    Member("swing_rod", "Q", "A", "link", "a"),
    Member("trigonal_ad", "A", "D", "truss", "c"),
    Member("trigonal_ae", "A", "E", "truss", "b"),
    Member("trigonal_de", "D", "E", "truss", "d"),
    Member("crank_2", "D", "R2", "cantilever", "q_2"),
    Member("piston_rod", "E", "P", "link", "e"),
)
"""The members whose cross-sections are sized, in a stable order.

Cranks are listed loaded-end first: a crank throw is a cantilever whose root is
the shaft, so the pin force at its tip is what sizes it.
"""

MEMBER_NAMES: tuple[str, ...] = tuple(member.name for member in MEMBERS)

#: Which members make up each moving body of the mechanism.
BODY_MEMBERS: dict[str, tuple[str, ...]] = {
    "crank_1": ("crank_1",),
    "swing_rod": ("swing_rod",),
    "trigonal": ("trigonal_ad", "trigonal_ae", "trigonal_de"),
    "crank_2": ("crank_2",),
    "piston_rod": ("piston_rod",),
}

BODY_NAMES: tuple[str, ...] = (
    "crank_1",
    "swing_rod",
    "trigonal",
    "crank_2",
    "piston_rod",
    "piston",
)
"""The six moving bodies, in the order their equations appear in the system."""


@dataclass(frozen=True)
class MassProperties:
    """Inertial description of every moving body, at every crank angle.

    All in the consistent unit system (mm, N, MPa, tonne, s).
    """

    member_mass: dict[str, float]
    """Mass of each sized member [tonne]."""

    body_mass: dict[str, float]
    """Total mass of each moving body [tonne]."""

    body_com: dict[str, FloatArray]
    """Centre-of-mass trajectory of each body, shaped ``(n_angles, 2)`` [mm]."""

    body_inertia: dict[str, float]
    """Rotational inertia about the body's own centre of mass [tonne mm^2]."""

    body_angle: dict[str, FloatArray]
    """Orientation history of each body [rad], for its angular acceleration."""

    @property
    def total_mass_kg(self) -> float:
        """Total moving mass [kg]."""
        return 1000.0 * float(sum(self.body_mass.values()))


@dataclass(frozen=True)
class DynamicLoads:
    """Joint reactions and output torque over one revolution, with inertia."""

    kinematics: Kinematics
    speed: float
    """Crankshaft speed ``Omega`` [rad/s]."""

    mass_properties: MassProperties

    reaction: dict[str, FloatArray]
    """Joint force histories, each ``(n_angles, 2)`` [N].

    Keys are ``R1``, ``Q``, ``A``, ``D``, ``R2``, ``E``, ``P``; each is the
    force transmitted *forward* along the chain, from the body nearer the
    crankshaft to the body nearer the piston.
    """

    liner_force: FloatArray
    """Side force the cylinder liner applies to the piston [N]."""

    liner_moment: FloatArray
    """Reaction moment on the piston from the liner [N.mm]."""

    gear_force: FloatArray
    """Gear tooth load ``T`` along the line of action [N]."""

    torque: FloatArray
    """Output torque on the crankshaft [N.mm]."""

    joint_acceleration: dict[str, FloatArray]
    """Acceleration of each joint, ``(n_angles, 2)`` [mm/s^2]."""

    body_acceleration: dict[str, FloatArray]
    """Centre-of-mass acceleration of each body [mm/s^2]."""

    body_angular_acceleration: dict[str, FloatArray]
    """Angular acceleration of each body [rad/s^2]."""

    gas_force: FloatArray
    """The applied gas force on the piston crown [N], as solved with."""

    matrix: FloatArray
    """The assembled equilibrium matrix, ``(n_angles, 18, 18)``.

    Kept so that :mod:`exlink.dynamics_jacobian` can differentiate the solve
    without reassembling or re-factorising it.
    """

    conditioning: float
    """Worst 2-norm condition number of the 18x18 equilibrium matrix.

    Large values mean the mechanism is near a configuration where it gains a
    degree of freedom -- exactly what the report's ``W`` constraint keeps it
    away from.
    """

    @property
    def mean_torque(self) -> float:
        """``M_r,ave`` over one revolution [N.mm]."""
        return float(np.mean(self.torque))


# -- indices into the 18-unknown vector -------------------------------------------
_JOINT_ORDER = ("R1", "Q", "A", "D", "R2", "E", "P")
_IDX = {name: (2 * i, 2 * i + 1) for i, name in enumerate(_JOINT_ORDER)}
_I_LINER_FORCE = 14
_I_LINER_MOMENT = 15
_I_GEAR = 16
_I_TORQUE = 17
_N_UNKNOWNS = 18


def mass_properties(
    kinematics: Kinematics,
    diameters: dict[str, float],
    density: float,
    piston_mass: float,
    spec: EngineSpec = DEFAULT_SPEC,
) -> MassProperties:
    """Assemble body masses, centres of mass and inertias from member sections.

    Every member is a solid round bar of the given diameter, so its mass and
    inertia follow from its length; a body's centre of mass is the mass-weighted
    mean of its members' midpoints, and its inertia the parallel-axis sum.

    Args:
        kinematics: A solved mechanism, supplying the joint trajectories.
        diameters: Section diameter of each member of :data:`MEMBERS` [mm].
        density: Material density [tonne/mm^3].
        piston_mass: Mass of the piston assembly [tonne].
        spec: Fixed engine data.

    Returns:
        The mass properties of all six moving bodies.
    """
    design = kinematics.design
    joints = {name: getattr(kinematics, name) for name in _JOINT_ORDER}

    member_mass: dict[str, float] = {}
    member_midpoint: dict[str, FloatArray] = {}
    member_own_inertia: dict[str, float] = {}
    for member in MEMBERS:
        length = abs(float(getattr(design, member.length_attribute)))
        diameter = float(diameters[member.name])
        area = np.pi * diameter**2 / 4.0
        mass = density * area * length
        member_mass[member.name] = mass
        member_midpoint[member.name] = 0.5 * (joints[member.start] + joints[member.end])
        # Solid cylinder about a transverse axis through its own centroid.
        member_own_inertia[member.name] = mass * (0.75 * diameter**2 + length**2) / 12.0

    body_mass: dict[str, float] = {}
    body_com: dict[str, FloatArray] = {}
    body_inertia: dict[str, float] = {}
    for body, names in BODY_MEMBERS.items():
        masses = np.array([member_mass[n] for n in names])
        total = float(masses.sum())
        midpoints = np.stack([member_midpoint[n] for n in names])  # (k, n, 2)
        com = np.tensordot(masses, midpoints, axes=(0, 0)) / total
        offsets = midpoints - com[None]
        inertia = float(
            sum(
                member_own_inertia[n] + member_mass[n] * float(np.mean(np.sum(o**2, axis=-1)))
                for n, o in zip(names, offsets, strict=True)
            )
        )
        body_mass[body] = total
        body_com[body] = com
        body_inertia[body] = inertia

    # The piston translates: its centre of mass sits half a piston length above
    # the wrist pin, and its rotational inertia never enters.
    body_mass["piston"] = piston_mass
    body_com["piston"] = joints["P"] + np.array([0.0, 0.5 * spec.piston_length])
    body_inertia["piston"] = 0.0

    zeros = np.zeros_like(kinematics.theta_1)
    body_angle = {
        "crank_1": kinematics.theta_1,
        "swing_rod": kinematics.theta_a,
        "trigonal": kinematics.theta_T,
        "crank_2": kinematics.theta_2,
        "piston_rod": kinematics.theta_e,
        "piston": zeros,
    }
    return MassProperties(
        member_mass=member_mass,
        body_mass=body_mass,
        body_com=body_com,
        body_inertia=body_inertia,
        body_angle=body_angle,
    )


def _second_derivative_2d(points: FloatArray) -> FloatArray:
    """Spectral second derivative of a periodic ``(n, 2)`` trajectory."""
    return np.stack(
        [spectral_derivative(points[:, 0], 2), spectral_derivative(points[:, 1], 2)],
        axis=-1,
    )


def solve(
    kinematics: Kinematics,
    gas_force: FloatArray,
    properties: MassProperties,
    speed: float,
    spec: EngineSpec = DEFAULT_SPEC,
) -> DynamicLoads:
    """Solve the mechanism's equilibrium including inertia.

    Args:
        kinematics: A solved mechanism.
        gas_force: Gas force on the piston crown at each crank angle [N],
            positive downwards.
        properties: Body masses, centres of mass and inertias.
        speed: Crankshaft speed ``Omega`` [rad/s].  Zero recovers the report's
            quasi-static result exactly.
        spec: Fixed engine data.

    Returns:
        Joint reactions, gear load and output torque over one revolution.
    """
    design = kinematics.design
    n = kinematics.theta_1.size
    joints = {name: getattr(kinematics, name) for name in _JOINT_ORDER}

    # -- accelerations ------------------------------------------------------------
    scale = speed**2
    joint_acceleration = {
        name: scale * _second_derivative_2d(points) for name, points in joints.items()
    }
    body_acceleration = {
        body: scale * _second_derivative_2d(com) for body, com in properties.body_com.items()
    }
    body_angular_acceleration = {
        body: scale * ramp_derivative(angle, 2) for body, angle in properties.body_angle.items()
    }
    # Both shafts turn at constant rate; force the exact zero rather than
    # leaving spectral round-off in the inertia couples.
    for body in ("crank_1", "crank_2", "piston"):
        body_angular_acceleration[body] = np.zeros(n)

    # -- gear mesh geometry -------------------------------------------------------
    theta_r = design.theta_r_rad
    alpha = spec.pressure_angle
    axis = np.array([np.cos(theta_r), np.sin(theta_r)])
    line_of_action = np.array(
        [np.cos(theta_r - np.pi / 2.0 + alpha), np.sin(theta_r - np.pi / 2.0 + alpha)]
    )
    contact_1 = np.zeros((n, 2)) + design.r_1 * axis
    contact_2 = joints["R2"] - design.r_2 * axis

    crown = np.stack(
        [np.full(n, design.x_1), kinematics.lam], axis=-1
    )  # centre of the piston crown, where the gas resultant acts

    matrix = np.zeros((n, _N_UNKNOWNS, _N_UNKNOWNS))
    rhs = np.zeros((n, _N_UNKNOWNS))

    def add_force(row: int, joint_index: tuple[int, int], sign: float) -> None:
        """Enter an unknown joint force into a body's two force equations."""
        matrix[:, row, joint_index[0]] += sign
        matrix[:, row + 1, joint_index[1]] += sign

    def add_moment(
        row: int,
        body: str,
        joint_index: tuple[int, int],
        sign: float,
        point: FloatArray,
    ) -> None:
        """Enter an unknown joint force into a body's moment equation."""
        arm = point - properties.body_com[body]
        matrix[:, row + 2, joint_index[0]] += -sign * arm[:, 1]
        matrix[:, row + 2, joint_index[1]] += sign * arm[:, 0]

    def add_joint(
        row: int, body: str, joint: str, sign: float, point: FloatArray | None = None
    ) -> None:
        """Enter an unknown joint force into all three of a body's equations."""
        index = _IDX[joint]
        location = joints[joint] if point is None else point
        add_force(row, index, sign)
        add_moment(row, body, index, sign, location)

    def set_inertia(row: int, body: str) -> None:
        """Put ``m a`` and ``I alpha`` on the right-hand side."""
        mass = properties.body_mass[body]
        rhs[:, row] += mass * body_acceleration[body][:, 0]
        rhs[:, row + 1] += mass * body_acceleration[body][:, 1]
        rhs[:, row + 2] += properties.body_inertia[body] * body_angular_acceleration[body]

    # -- body 1: the crankshaft's crank throw ------------------------------------
    add_joint(0, "crank_1", "R1", +1.0)
    add_joint(0, "crank_1", "Q", -1.0)
    matrix[:, 0, _I_GEAR] += -line_of_action[0]
    matrix[:, 1, _I_GEAR] += -line_of_action[1]
    arm = contact_1 - properties.body_com["crank_1"]
    matrix[:, 2, _I_GEAR] += -(arm[:, 0] * line_of_action[1] - arm[:, 1] * line_of_action[0])
    # Signed so that the unknown is the torque delivered *out* of the
    # crankshaft, matching the report's convention: the moment equation reads
    # ``sum(M) - M_r = 0``, so a positive M_r is useful work out.
    matrix[:, 2, _I_TORQUE] += -1.0
    set_inertia(0, "crank_1")

    # -- body 2: the swing rod ----------------------------------------------------
    add_joint(3, "swing_rod", "Q", +1.0)
    add_joint(3, "swing_rod", "A", -1.0)
    set_inertia(3, "swing_rod")

    # -- body 3: the trigonal link ------------------------------------------------
    add_joint(6, "trigonal", "A", +1.0)
    add_joint(6, "trigonal", "D", +1.0)
    add_joint(6, "trigonal", "E", -1.0)
    set_inertia(6, "trigonal")

    # -- body 4: the eccentric shaft's crank throw --------------------------------
    add_joint(9, "crank_2", "R2", +1.0)
    add_joint(9, "crank_2", "D", -1.0)
    matrix[:, 9, _I_GEAR] += line_of_action[0]
    matrix[:, 10, _I_GEAR] += line_of_action[1]
    arm = contact_2 - properties.body_com["crank_2"]
    matrix[:, 11, _I_GEAR] += arm[:, 0] * line_of_action[1] - arm[:, 1] * line_of_action[0]
    set_inertia(9, "crank_2")

    # -- body 5: the piston rod ---------------------------------------------------
    add_joint(12, "piston_rod", "E", +1.0)
    add_joint(12, "piston_rod", "P", -1.0)
    set_inertia(12, "piston_rod")

    # -- body 6: the piston -------------------------------------------------------
    add_joint(15, "piston", "P", +1.0)
    matrix[:, 15, _I_LINER_FORCE] += 1.0
    arm = crown - properties.body_com["piston"]
    # The liner's normal force acts along x at the guide; its line of action is
    # not resolved, so the guide also carries a reaction moment.
    matrix[:, 17, _I_LINER_MOMENT] += 1.0
    set_inertia(15, "piston")
    # Known applied loads move to the right-hand side with a sign change: the
    # gas resultant is (0, -gas_force), so it enters as +gas_force.
    rhs[:, 16] += gas_force
    rhs[:, 17] += -(arm[:, 0] * (-gas_force) - arm[:, 1] * 0.0)

    # ``rhs`` needs a trailing axis: with a stack of matrices, numpy reads a
    # 2-D right-hand side as a batch of columns rather than a batch of vectors.
    solution: FloatArray = np.linalg.solve(matrix, rhs[..., None])[..., 0].astype(
        np.float64, copy=False
    )
    conditioning = float(np.max(np.linalg.cond(matrix)))

    reaction = {name: solution[:, _IDX[name][0] : _IDX[name][1] + 1] for name in _JOINT_ORDER}
    return DynamicLoads(
        kinematics=kinematics,
        speed=speed,
        mass_properties=properties,
        reaction=reaction,
        liner_force=solution[:, _I_LINER_FORCE],
        liner_moment=solution[:, _I_LINER_MOMENT],
        gear_force=solution[:, _I_GEAR],
        torque=solution[:, _I_TORQUE],
        gas_force=np.asarray(gas_force, dtype=float),
        matrix=matrix,
        joint_acceleration=joint_acceleration,
        body_acceleration=body_acceleration,
        body_angular_acceleration=body_angular_acceleration,
        conditioning=conditioning,
    )
