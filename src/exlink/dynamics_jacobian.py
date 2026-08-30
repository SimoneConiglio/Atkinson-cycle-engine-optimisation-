"""Exact derivatives through the sizing / dynamics coupling.

Without these, the coupled problem has to be differentiated by differencing the
whole MDA: eleven converged fixed points per gradient, each some fifty sweeps,
which measured at ten to twenty minutes for a short optimizer run.  With them,
GEMSEO assembles the coupled derivative from each discipline's *local*
Jacobian and the cost collapses to a handful of single evaluations plus one
small linear solve.

Three things have to be differentiated, and each has a clean route:

**The spectral operator.**  Accelerations are ``Omega^2`` times a second
derivative with respect to crank angle, and that derivative is a *linear*
operator.  So the derivative of an acceleration is the same operator applied to
the derivative of the trajectory -- no new mathematics, just
``spectral_derivative`` applied to the arrays :mod:`exlink.jacobian` already
produces.

**The 18x18 equilibrium solve.**  For ``A x = b``,

.. math:: \\frac{\\partial x}{\\partial p} = A^{-1}\\left(
          \\frac{\\partial b}{\\partial p} - \\frac{\\partial A}{\\partial p} x\\right)

and the factorisation of ``A`` is already available from the solve itself.  The
entries of ``A`` are moment arms -- differences of joint and centre-of-mass
positions -- so ``dA/dp`` follows directly from the position derivatives.

**The sizing bisection.**  The diameter is defined implicitly by driving the
worst utilisation to one, ``F(d, N, M) = 0``.  The implicit function theorem
gives ``dd/dN = -(dF/dN)/(dF/dd)`` without differentiating the bisection at all.

Parameter layout
----------------
Derivative arrays carry a trailing axis of length ``N_PARAMETERS``: the eleven
design variables first, then the seven section diameters.  :data:`DESIGN_SLICE`
and :data:`DIAMETER_SLICE` name the two halves.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec
from .derivatives import spectral_derivative
from .design import VARIABLE_NAMES, Design
from .dynamics import BODY_MEMBERS, MEMBER_NAMES, MEMBERS, DynamicLoads, MassProperties
from .jacobian import INDEX, KinematicJacobian
from .kinematics import Kinematics
from .materials import Material, SafetyFactors

FloatArray = NDArray[np.float64]

N_DESIGN = len(VARIABLE_NAMES)
N_DIAMETERS = len(MEMBER_NAMES)
N_PARAMETERS = N_DESIGN + N_DIAMETERS
DESIGN_SLICE = slice(0, N_DESIGN)
DIAMETER_SLICE = slice(N_DESIGN, N_PARAMETERS)

_JOINT_ORDER = ("R1", "Q", "A", "D", "R2", "E", "P")


def member_length_jacobian(design: Design) -> FloatArray:
    """``d L_k / dX`` for every sized member, shaped ``(n_members, 11)``.

    Most members take their length straight from a design variable.  The two
    sides of the trigonal link that are *derived* -- ``b = hypot(x_b, y_b)`` and
    ``d = hypot(x_b - c, y_b)`` -- carry the chain rule of that
    reparametrisation.
    """
    jacobian = np.zeros((N_DIAMETERS, N_DESIGN))
    b, side = design.b, design.d
    for row, member in enumerate(MEMBERS):
        attribute = member.length_attribute
        if attribute in VARIABLE_NAMES:
            jacobian[row, INDEX[attribute]] = np.sign(getattr(design, attribute)) or 1.0
        elif attribute == "b":
            jacobian[row, INDEX["x_b"]] = design.x_b / b
            jacobian[row, INDEX["y_b"]] = design.y_b / b
        elif attribute == "d":
            jacobian[row, INDEX["x_b"]] = (design.x_b - design.c) / side
            jacobian[row, INDEX["c"]] = -(design.x_b - design.c) / side
            jacobian[row, INDEX["y_b"]] = design.y_b / side
        else:  # pragma: no cover - guards a future member with a new length
            msg = f"no length derivative for {attribute!r}"
            raise ValueError(msg)
    return jacobian


@dataclass(frozen=True)
class MassJacobian:
    """Derivatives of the inertial description with respect to ``(X, d)``."""

    member_mass: FloatArray
    """``d m_k / dp``, shaped ``(n_members, N_PARAMETERS)``."""

    body_mass: dict[str, FloatArray]
    """``d M_i / dp``, each shaped ``(N_PARAMETERS,)``."""

    body_com: dict[str, FloatArray]
    """``d G_i / dp``, each shaped ``(n_angles, 2, N_PARAMETERS)``."""

    body_inertia: dict[str, FloatArray]
    """``d I_i / dp``, each shaped ``(N_PARAMETERS,)``."""


def mass_property_jacobian(
    design: Design,
    kinematics: Kinematics,
    kinematic: KinematicJacobian,
    properties: MassProperties,
    diameters: dict[str, float],
    density: float,
) -> MassJacobian:
    """Differentiate the body masses, centres of mass and inertias.

    Args:
        design: The mechanism dimensions.
        kinematics: The solved kinematics.
        kinematic: Its derivatives from :func:`exlink.jacobian.kinematic_jacobian`.
        properties: The mass properties being differentiated.
        diameters: Section diameter of each member [mm].
        density: Material density [tonne/mm^3].

    Returns:
        The derivatives of every inertial quantity.
    """
    joints = {name: getattr(kinematics, name) for name in _JOINT_ORDER}
    n = kinematics.theta_1.size
    length_jacobian = member_length_jacobian(design)

    lengths = np.array([abs(float(getattr(design, m.length_attribute))) for m in MEMBERS])
    widths = np.array([float(diameters[m.name]) for m in MEMBERS])

    # m_k = rho * (pi d_k^2 / 4) * L_k
    d_member_mass = np.zeros((N_DIAMETERS, N_PARAMETERS))
    area = np.pi * widths**2 / 4.0
    d_member_mass[:, DESIGN_SLICE] = density * area[:, None] * length_jacobian
    for row in range(N_DIAMETERS):
        d_member_mass[row, N_DESIGN + row] = density * np.pi * widths[row] / 2.0 * lengths[row]

    # Midpoints move with the linkage but not with the sections.
    d_midpoint: dict[str, FloatArray] = {}
    for row, member in enumerate(MEMBERS):
        stacked = np.zeros((n, 2, N_PARAMETERS))
        stacked[..., DESIGN_SLICE] = 0.5 * (
            kinematic.joints[member.start] + kinematic.joints[member.end]
        )
        d_midpoint[member.name] = stacked
        del row

    midpoint = {m.name: 0.5 * (joints[m.start] + joints[m.end]) for m in MEMBERS}

    d_body_mass: dict[str, FloatArray] = {}
    d_body_com: dict[str, FloatArray] = {}
    d_body_inertia: dict[str, FloatArray] = {}
    index_of = {m.name: i for i, m in enumerate(MEMBERS)}

    for body, names in BODY_MEMBERS.items():
        rows = [index_of[name] for name in names]
        masses = np.array([properties.member_mass[name] for name in names])
        total = float(masses.sum())
        d_total = d_member_mass[rows].sum(axis=0)
        d_body_mass[body] = d_total

        com = properties.body_com[body]
        weighted = sum(
            d_member_mass[row][None, None, :] * midpoint[name][..., None]
            + masses[k] * d_midpoint[name]
            for k, (row, name) in enumerate(zip(rows, names, strict=True))
        )
        d_com = weighted / total - com[..., None] * d_total[None, None, :] / total
        d_body_com[body] = d_com

        # I = sum_k [ m_k (3 d_k^2/4 + L_k^2)/12 + m_k |mid_k - G|^2 ]
        derivative = np.zeros(N_PARAMETERS)
        for k, (row, name) in enumerate(zip(rows, names, strict=True)):
            own = (0.75 * widths[row] ** 2 + lengths[row] ** 2) / 12.0
            d_own = np.zeros(N_PARAMETERS)
            d_own[DESIGN_SLICE] = 2.0 * lengths[row] * length_jacobian[row] / 12.0
            d_own[N_DESIGN + row] = 1.5 * widths[row] / 12.0

            offset = midpoint[name] - com
            d_offset = d_midpoint[name] - d_com
            squared = float(np.mean(np.sum(offset**2, axis=-1)))
            d_squared = 2.0 * np.mean(np.sum(offset[..., None] * d_offset, axis=1), axis=0)

            derivative += d_member_mass[row] * (own + squared) + masses[k] * (d_own + d_squared)
        d_body_inertia[body] = derivative

    # The piston: its mass is fixed by the gas pressure, not by any section, and
    # it translates, so only its centre of mass moves.
    piston_com = np.zeros((n, 2, N_PARAMETERS))
    piston_com[..., DESIGN_SLICE] = kinematic.joints["P"]
    d_body_com["piston"] = piston_com
    d_body_mass["piston"] = np.zeros(N_PARAMETERS)
    d_body_inertia["piston"] = np.zeros(N_PARAMETERS)

    return MassJacobian(
        member_mass=d_member_mass,
        body_mass=d_body_mass,
        body_com=d_body_com,
        body_inertia=d_body_inertia,
    )


def _spectral_2d(values: FloatArray) -> FloatArray:
    """Second crank-angle derivative of a ``(n, 2, ...)`` array, along axis 0."""
    moved = np.moveaxis(values, 0, -1)
    return np.moveaxis(spectral_derivative(moved, 2), -1, 0)


@dataclass(frozen=True)
class AccelerationJacobian:
    """Derivatives of the acceleration field with respect to ``(X, d)``."""

    joint: dict[str, FloatArray]
    """``d a_J / dp``, each ``(n_angles, 2, N_PARAMETERS)``."""

    body: dict[str, FloatArray]
    """``d a_Gi / dp``, each ``(n_angles, 2, N_PARAMETERS)``."""

    angular: dict[str, FloatArray]
    """``d alpha_i / dp``, each ``(n_angles, N_PARAMETERS)``."""


def acceleration_jacobian(
    kinematic: KinematicJacobian,
    mass_jacobian: MassJacobian,
    speed: float,
    n_angles: int,
) -> AccelerationJacobian:
    """Differentiate the accelerations, by reusing the spectral operator.

    ``a = Omega^2 D^2 r`` with ``D^2`` linear, so ``da/dp = Omega^2 D^2 (dr/dp)``.

    Args:
        kinematic: Derivatives of the kinematics.
        mass_jacobian: Derivatives of the centres of mass.
        speed: Crankshaft speed ``Omega`` [rad/s].
        n_angles: Crank angles per revolution.

    Returns:
        The derivatives of every acceleration.
    """
    scale = speed**2
    joint = {}
    for name, derivative in kinematic.joints.items():
        padded = np.zeros((n_angles, 2, N_PARAMETERS))
        padded[..., DESIGN_SLICE] = derivative
        joint[name] = scale * _spectral_2d(padded)

    body = {
        name: scale * _spectral_2d(derivative)
        for name, derivative in mass_jacobian.body_com.items()
    }

    # Both shafts and the piston turn or translate at constant rate: their
    # angular accelerations are identically zero, derivative included.
    angular: dict[str, FloatArray] = {}
    for name in ("crank_1", "crank_2", "piston"):
        angular[name] = np.zeros((n_angles, N_PARAMETERS))
    for name, derivative in (
        ("swing_rod", kinematic.theta_a),
        ("trigonal", kinematic.theta_T),
        ("piston_rod", kinematic.theta_e),
    ):
        padded = np.zeros((n_angles, N_PARAMETERS))
        padded[:, DESIGN_SLICE] = derivative
        angular[name] = scale * spectral_derivative(padded.T, 2).T

    return AccelerationJacobian(joint=joint, body=body, angular=angular)


# -- indices into the 18-unknown vector, mirroring exlink.dynamics ----------------
_IDX = {name: (2 * i, 2 * i + 1) for i, name in enumerate(_JOINT_ORDER)}
_I_LINER_FORCE = 14
_I_LINER_MOMENT = 15
_I_GEAR = 16
_I_TORQUE = 17
_N_UNKNOWNS = 18

#: ``(row, body, joint, sign)`` for every unknown joint force in the system,
#: in the same order the forward assembly enters them.
_JOINT_TERMS: tuple[tuple[int, str, str, float], ...] = (
    (0, "crank_1", "R1", +1.0),
    (0, "crank_1", "Q", -1.0),
    (3, "swing_rod", "Q", +1.0),
    (3, "swing_rod", "A", -1.0),
    (6, "trigonal", "A", +1.0),
    (6, "trigonal", "D", +1.0),
    (6, "trigonal", "E", -1.0),
    (9, "crank_2", "R2", +1.0),
    (9, "crank_2", "D", -1.0),
    (12, "piston_rod", "E", +1.0),
    (12, "piston_rod", "P", -1.0),
    (15, "piston", "P", +1.0),
)

#: ``(row, body)`` for each body's three equilibrium equations.
_BODY_ROWS: tuple[tuple[int, str], ...] = (
    (0, "crank_1"),
    (3, "swing_rod"),
    (6, "trigonal"),
    (9, "crank_2"),
    (12, "piston_rod"),
    (15, "piston"),
)


def equilibrium_jacobian(
    design: Design,
    kinematics: Kinematics,
    kinematic: KinematicJacobian,
    properties: MassProperties,
    mass_jacobian: MassJacobian,
    acceleration: AccelerationJacobian,
    loads: DynamicLoads,
    gas_force_derivative: FloatArray,
    spec: EngineSpec = DEFAULT_SPEC,
) -> FloatArray:
    """Differentiate the 18x18 equilibrium solve.

    For ``A x = b`` the derivative is ``dx/dp = A^-1 (db/dp - dA/dp x)``, so the
    same matrix the forward solve already assembled is reused -- no second
    factorisation, and no differentiation of the solve itself.

    The entries of ``A`` are moment arms, differences of joint and
    centre-of-mass positions, so ``dA/dp`` follows straight from the position
    derivatives.  ``b`` holds ``m a`` and ``I alpha`` and the applied gas load,
    so ``db/dp`` follows from the mass, acceleration and gas-force derivatives.

    Args:
        design: The mechanism dimensions.
        kinematics: The solved kinematics.
        kinematic: Its derivatives.
        properties: The mass properties.
        mass_jacobian: Their derivatives.
        acceleration: The acceleration derivatives.
        loads: The solved load case, supplying ``x`` and the matrix.
        gas_force_derivative: ``d P_gas / dX``, shaped ``(n_angles, 11)``.
        spec: Fixed engine data.

    Returns:
        ``dx/dp``, shaped ``(n_angles, 18, N_PARAMETERS)``, in the unknown order
        of :mod:`exlink.dynamics`.
    """
    n = kinematics.theta_1.size
    # Joint position derivatives padded onto the full parameter axis.
    d_joint = {}
    for name in _JOINT_ORDER:
        padded = np.zeros((n, 2, N_PARAMETERS))
        padded[..., DESIGN_SLICE] = kinematic.joints[name]
        d_joint[name] = padded

    d_matrix = np.zeros((n, _N_UNKNOWNS, _N_UNKNOWNS, N_PARAMETERS))
    d_rhs = np.zeros((n, _N_UNKNOWNS, N_PARAMETERS))

    # -- moment arms of the unknown joint forces ---------------------------------
    # Only the moment rows depend on the design; the two force rows are +-1.
    for row, body, joint, sign in _JOINT_TERMS:
        first, second = _IDX[joint]
        arm = d_joint[joint] - mass_jacobian.body_com[body]
        d_matrix[:, row + 2, first, :] = -sign * arm[:, 1, :]
        d_matrix[:, row + 2, second, :] = sign * arm[:, 0, :]

    # -- the gear mesh ------------------------------------------------------------
    theta_r = design.theta_r_rad
    alpha = spec.pressure_angle
    angle = theta_r - np.pi / 2.0 + alpha
    line = np.array([np.cos(angle), np.sin(angle)])
    d_line = np.zeros((2, N_PARAMETERS))
    d_line[0, INDEX["theta_r"]] = -np.sin(angle) * np.pi / 180.0
    d_line[1, INDEX["theta_r"]] = np.cos(angle) * np.pi / 180.0

    # Both gears touch at the same point, (2 I / 3) along the shaft axis.
    axis = np.array([np.cos(theta_r), np.sin(theta_r)])
    d_axis = np.zeros((2, N_PARAMETERS))
    d_axis[0, INDEX["theta_r"]] = -np.sin(theta_r) * np.pi / 180.0
    d_axis[1, INDEX["theta_r"]] = np.cos(theta_r) * np.pi / 180.0
    radius = design.r_1
    d_radius = np.zeros(N_PARAMETERS)
    d_radius[INDEX["I"]] = 2.0 / 3.0
    contact = radius * axis
    d_contact = radius * d_axis + np.outer(axis, d_radius)

    for row, body, sign in ((0, "crank_1", -1.0), (9, "crank_2", +1.0)):
        d_matrix[:, row, _I_GEAR, :] = sign * d_line[0][None, :]
        d_matrix[:, row + 1, _I_GEAR, :] = sign * d_line[1][None, :]
        arm_value = contact[None, :] - properties.body_com[body]
        d_arm = d_contact[None, :, :] - mass_jacobian.body_com[body]
        cross = (
            d_arm[:, 0, :] * line[1]
            + arm_value[:, 0, None] * d_line[1][None, :]
            - d_arm[:, 1, :] * line[0]
            - arm_value[:, 1, None] * d_line[0][None, :]
        )
        d_matrix[:, row + 2, _I_GEAR, :] = sign * cross

    # -- inertia on the right-hand side -------------------------------------------
    for row, body in _BODY_ROWS:
        mass = properties.body_mass[body]
        d_mass = mass_jacobian.body_mass[body]
        body_acceleration = loads.body_acceleration[body]
        d_rhs[:, row, :] = (
            d_mass[None, :] * body_acceleration[:, 0, None]
            + mass * acceleration.body[body][:, 0, :]
        )
        d_rhs[:, row + 1, :] = (
            d_mass[None, :] * body_acceleration[:, 1, None]
            + mass * acceleration.body[body][:, 1, :]
        )
        angular = loads.body_angular_acceleration[body]
        d_rhs[:, row + 2, :] = (
            mass_jacobian.body_inertia[body][None, :] * angular[:, None]
            + properties.body_inertia[body] * acceleration.angular[body]
        )

    # -- the applied gas load on the piston ---------------------------------------
    d_gas = np.zeros((n, N_PARAMETERS))
    d_gas[:, DESIGN_SLICE] = gas_force_derivative
    force = loads.gas_force
    crown = np.stack([np.full(n, design.x_1), kinematics.lam], axis=-1)
    d_crown = np.zeros((n, 2, N_PARAMETERS))
    d_crown[:, 0, INDEX["x_1"]] = 1.0
    d_crown[:, 1, DESIGN_SLICE] = kinematic.lam
    arm_value = crown - properties.body_com["piston"]
    d_arm = d_crown - mass_jacobian.body_com["piston"]

    d_rhs[:, 16, :] += d_gas
    d_rhs[:, 17, :] += d_arm[:, 0, :] * force[:, None] + arm_value[:, 0, None] * d_gas

    solution = np.stack(
        [loads.reaction[name][:, component] for name in _JOINT_ORDER for component in (0, 1)]
        + [loads.liner_force, loads.liner_moment, loads.gear_force, loads.torque],
        axis=-1,
    )
    residual = d_rhs - np.einsum("nijp,nj->nip", d_matrix, solution)
    return np.linalg.solve(loads.matrix, residual)


def _cross(first: FloatArray, second: FloatArray) -> FloatArray:
    """z-component of a 2-D cross product, over a trailing parameter axis."""
    return first[..., 0, :] * second[..., 1] - first[..., 1, :] * second[..., 0]


def member_load_jacobian(
    design: Design,
    kinematics: Kinematics,
    kinematic: KinematicJacobian,
    properties: MassProperties,
    mass_jacobian: MassJacobian,
    acceleration: AccelerationJacobian,
    loads: DynamicLoads,
    reaction_jacobian: FloatArray,
    stations: int,
) -> tuple[FloatArray, FloatArray]:
    """Differentiate the internal axial force and bending moment of every member.

    The closed forms of :func:`exlink.sizing.internal_loads` are differentiated
    directly.  For the trigonal link's three sides the pin-jointed idealisation
    is differentiated too: the joint-equilibrium least squares gives
    ``dz = (M^T M)^-1 M^T (dr - dM z)``, exact because the system is consistent.

    Args:
        design: The mechanism dimensions.
        kinematics: The solved kinematics.
        kinematic: Its derivatives.
        properties: The mass properties.
        mass_jacobian: Their derivatives.
        acceleration: The acceleration derivatives.
        loads: The solved load case.
        reaction_jacobian: ``dx/dp`` from :func:`equilibrium_jacobian`.
        stations: Sections evaluated along each member.

    Returns:
        ``(d_axial, d_bending)``, each shaped
        ``(n_members, n_angles, stations, N_PARAMETERS)``.
    """
    from .sizing import MEMBER_IS_SLENDER  # noqa: F401  (kept for symmetry)

    n = kinematics.theta_1.size
    joints = {name: getattr(kinematics, name) for name in _JOINT_ORDER}
    d_joint = {}
    for name in _JOINT_ORDER:
        padded = np.zeros((n, 2, N_PARAMETERS))
        padded[..., DESIGN_SLICE] = kinematic.joints[name]
        d_joint[name] = padded

    def reaction_slice(name: str) -> FloatArray:
        first, second = _IDX[name]
        return reaction_jacobian[:, first : second + 1, :]

    index_of = {m.name: i for i, m in enumerate(MEMBERS)}

    # -- the trigonal link, as a pin-jointed triangle ------------------------------
    sides = [m for m in MEMBERS if m.kind == "truss"]
    vertices = ("A", "D", "E")
    external = {
        "A": (loads.reaction["A"], reaction_slice("A")),
        "D": (loads.reaction["D"], reaction_slice("D")),
        "E": (-loads.reaction["E"], -reaction_slice("E")),
    }
    matrix = np.zeros((n, 6, 3))
    d_matrix = np.zeros((n, 6, 3, N_PARAMETERS))
    rhs = np.zeros((n, 6))
    d_rhs = np.zeros((n, 6, N_PARAMETERS))
    for row, vertex in enumerate(vertices):
        attached = sum(
            properties.member_mass[m.name] for m in sides if vertex in (m.start, m.end)
        )
        d_attached = sum(
            mass_jacobian.member_mass[index_of[m.name]]
            for m in sides
            if vertex in (m.start, m.end)
        )
        lumped = -0.5 * attached * loads.joint_acceleration[vertex]
        d_lumped = -0.5 * (
            d_attached[None, None, :] * loads.joint_acceleration[vertex][..., None]
            + attached * acceleration.joint[vertex]
        )
        for column, side in enumerate(sides):
            if vertex not in (side.start, side.end):
                continue
            other = side.end if vertex == side.start else side.start
            delta = joints[other] - joints[vertex]
            d_delta = d_joint[other] - d_joint[vertex]
            length = np.linalg.norm(delta, axis=-1, keepdims=True)
            unit = delta / length
            d_length = np.sum(unit[..., None] * d_delta, axis=1, keepdims=True)
            d_unit = d_delta / length[..., None] - (
                delta[..., None] * d_length / length[..., None] ** 2
            )
            matrix[:, 2 * row, column] = unit[:, 0]
            matrix[:, 2 * row + 1, column] = unit[:, 1]
            d_matrix[:, 2 * row, column, :] = d_unit[:, 0, :]
            d_matrix[:, 2 * row + 1, column, :] = d_unit[:, 1, :]
        value, derivative = external[vertex]
        rhs[:, 2 * row] = -(value[:, 0] + lumped[:, 0])
        rhs[:, 2 * row + 1] = -(value[:, 1] + lumped[:, 1])
        d_rhs[:, 2 * row, :] = -(derivative[:, 0, :] + d_lumped[:, 0, :])
        d_rhs[:, 2 * row + 1, :] = -(derivative[:, 1, :] + d_lumped[:, 1, :])

    normal = np.einsum("nij,nik->njk", matrix, matrix)
    projected = np.einsum("nij,ni->nj", matrix, rhs)
    axial_truss = np.linalg.solve(normal, projected[..., None])[..., 0]
    # z = (M^T M)^-1 M^T r, so dz = (M^T M)^-1 [dM^T (r - M z) + M^T (dr - dM z)].
    # The joint-equilibrium system is consistent -- the vertex loads already
    # balance the link's own inertia -- so the first term is zero to round-off;
    # it is kept anyway, since it costs nothing and guards the assumption.
    residual = d_rhs - np.einsum("nijp,nj->nip", d_matrix, axial_truss)
    d_projected = np.einsum("nij,nip->njp", matrix, residual) + np.einsum(
        "nijp,ni->njp", d_matrix, rhs - np.einsum("nij,nj->ni", matrix, axial_truss)
    )
    d_axial_truss = np.linalg.solve(normal, d_projected)

    # -- the applied end force of every member -------------------------------------
    applied = {
        "crank_1": (-loads.reaction["Q"], -reaction_slice("Q")),
        "swing_rod": (loads.reaction["Q"], reaction_slice("Q")),
        "crank_2": (-loads.reaction["D"], -reaction_slice("D")),
        "piston_rod": (loads.reaction["E"], reaction_slice("E")),
    }

    d_axial = np.zeros((len(MEMBERS), n, stations, N_PARAMETERS))
    d_bending = np.zeros((len(MEMBERS), n, stations, N_PARAMETERS))
    grid = np.linspace(0.0, 1.0, stations)[None, :]

    for row, member in enumerate(MEMBERS):
        start, end = joints[member.start], joints[member.end]
        d_start, d_end = d_joint[member.start], d_joint[member.end]
        delta = end - start
        d_delta = d_end - d_start
        length = np.linalg.norm(delta, axis=-1, keepdims=True)
        unit = delta / length
        d_length = np.sum(unit[..., None] * d_delta, axis=1, keepdims=True)
        d_unit = d_delta / length[..., None] - (
            delta[..., None] * d_length / length[..., None] ** 2
        )

        first = loads.joint_acceleration[member.start]
        second = loads.joint_acceleration[member.end]
        d_first = acceleration.joint[member.start]
        d_second = acceleration.joint[member.end]
        mass = properties.member_mass[member.name]
        d_mass = mass_jacobian.member_mass[row]

        if member.kind == "truss":
            column = [m.name for m in sides].index(member.name)
            normal_direction = np.stack([-unit[:, 1], unit[:, 0]], axis=-1)
            d_normal = np.stack([-d_unit[:, 1, :], d_unit[:, 0, :]], axis=1)
            first_n = np.sum(first * normal_direction, axis=-1)
            second_n = np.sum(second * normal_direction, axis=-1)
            d_first_n = np.sum(
                d_first * normal_direction[..., None] + first[..., None] * d_normal,
                axis=1,
            )
            d_second_n = np.sum(
                d_second * normal_direction[..., None] + second[..., None] * d_normal,
                axis=1,
            )
            shear = -mass * (first_n / 2.0 + (second_n - first_n) / 6.0)
            d_shear = -(
                d_mass[None, :] * (first_n / 2.0 + (second_n - first_n) / 6.0)[:, None]
                + mass * (d_first_n / 2.0 + (d_second_n - d_first_n) / 6.0)
            )
            value = axial_truss[:, column]
            derivative = d_axial_truss[:, column, :]
            force = -value[:, None] * unit + shear[:, None] * normal_direction
            d_force = (
                -derivative[:, None, :] * unit[..., None]
                - value[:, None, None] * d_unit
                + d_shear[:, None, :] * normal_direction[..., None]
                + shear[:, None, None] * d_normal
            )
        else:
            force, d_force = applied[member.name]

        delta_a = second - first
        d_delta_a = d_second - d_first

        # F(s) = m [a1 s + da s^2/2] - F1
        internal = (
            d_mass[None, None, None, :]
            * (
                first[:, None, :, None] * grid[..., None, None]
                + delta_a[:, None, :, None] * grid[..., None, None] ** 2 / 2.0
            )
            + mass
            * (
                d_first[:, None, :, :] * grid[..., None, None]
                + d_delta_a[:, None, :, :] * grid[..., None, None] ** 2 / 2.0
            )
            - d_force[:, None, :, :]
        )
        value_internal = (
            mass
            * (
                first[:, None, :] * grid[..., None]
                + delta_a[:, None, :] * grid[..., None] ** 2 / 2.0
            )
            - force[:, None, :]
        )
        d_axial[row] = np.sum(
            internal * unit[:, None, :, None] + value_internal[..., None] * d_unit[:, None],
            axis=2,
        )

        # M(s) = -m [(dr x a1) s^2/2 + (dr x da) s^3/6] + s (dr x F1)
        def cross_pair(u: FloatArray, v: FloatArray) -> FloatArray:
            return u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]

        def d_cross_pair(
            u: FloatArray, du: FloatArray, v: FloatArray, dv: FloatArray
        ) -> FloatArray:
            return (
                du[:, 0, :] * v[:, 1, None]
                + u[:, 0, None] * dv[:, 1, :]
                - du[:, 1, :] * v[:, 0, None]
                - u[:, 1, None] * dv[:, 0, :]
            )

        cross_a1 = cross_pair(delta, first)
        cross_da = cross_pair(delta, delta_a)
        d_cross_a1 = d_cross_pair(delta, d_delta, first, d_first)
        d_cross_da = d_cross_pair(delta, d_delta, delta_a, d_delta_a)
        d_cross_f = d_cross_pair(delta, d_delta, force, d_force)

        d_bending[row] = (
            -d_mass[None, None, :]
            * (
                cross_a1[:, None, None] * grid[..., None] ** 2 / 2.0
                + cross_da[:, None, None] * grid[..., None] ** 3 / 6.0
            )
            - mass
            * (
                d_cross_a1[:, None, :] * grid[..., None] ** 2 / 2.0
                + d_cross_da[:, None, :] * grid[..., None] ** 3 / 6.0
            )
            + d_cross_f[:, None, :] * grid[..., None]
        )

    del design
    return d_axial, d_bending


@dataclass(frozen=True)
class SizingJacobian:
    """Derivatives of the required diameters, with respect to what sets them."""

    d_axial: FloatArray
    """``d(diameter_k) / d(axial)``, shaped ``(n_members, n_angles, stations)``.

    Sparse in practice: only the sections and crank angles that *attain* the
    binding utilisation carry a non-zero entry.
    """

    d_bending: FloatArray
    """``d(diameter_k) / d(bending)``, same shape."""

    d_length: FloatArray
    """``d(diameter_k) / d(L_k)``, shaped ``(n_members,)``.

    Non-zero only where buckling binds; yield and fatigue do not see the span.
    """


def sizing_jacobian(
    axial: FloatArray,
    bending: FloatArray,
    diameters: FloatArray,
    lengths: FloatArray,
    material: Material,
    safety: SafetyFactors,
) -> SizingJacobian:
    """Differentiate the sizing solve, without differentiating the bisection.

    The diameter is defined implicitly: it is the value that drives the worst
    utilisation to one, ``U(d, N, M, L) = 1``.  The implicit function theorem
    then gives

    .. math:: \\frac{\\partial d}{\\partial q} =
              -\\frac{\\partial U/\\partial q}{\\partial U/\\partial d}

    so the bisection itself never has to be differentiated -- only the closed
    form of whichever failure mode binds.  Each mode is an extremum over fibre,
    crank angle and station, so the envelope theorem applies again: the
    derivative is taken at the entry that attains it.

    Args:
        axial: ``(n_members, n_angles, n_stations)`` [N], tension positive.
        bending: Same shape [N.mm].
        diameters: The solved diameters ``(n_members,)`` [mm].
        lengths: Member lengths ``(n_members,)`` [mm].
        material: The :class:`~exlink.materials.Material`.
        safety: The :class:`~exlink.materials.SafetyFactors`.

    Returns:
        The derivative of each diameter with respect to its loads and length.
    """
    from .sizing import MEMBER_FIXITY, _utilisations

    n_members = diameters.size
    area = np.pi * diameters**2 / 4.0
    modulus = np.pi * diameters**3 / 32.0

    direct = axial / area[:, None, None]
    flexural = bending / modulus[:, None, None]
    fibres = np.stack([direct + flexural, direct - flexural], axis=0)
    signs = np.array([1.0, -1.0])

    static, fatigue, buckling = _utilisations(
        diameters, axial, bending, lengths, MEMBER_FIXITY, material, safety
    )

    d_axial = np.zeros_like(axial)
    d_bending = np.zeros_like(bending)
    d_length = np.zeros(n_members)

    # dU/dd is common to all three modes through 1/A ~ d^-2 and 1/Z ~ d^-3.
    for member in range(n_members):
        diameter = diameters[member]
        mode = int(np.argmax([static[member], fatigue[member], buckling[member]]))
        grad_axial = np.zeros(axial.shape[1:])
        grad_bending = np.zeros(bending.shape[1:])
        grad_length = 0.0

        if mode == 0:  # static yield
            flat = int(np.argmax(np.abs(fibres[:, member])))
            fibre, angle, station = (
                int(v) for v in np.unravel_index(flat, fibres[:, member].shape)
            )
            stress = fibres[fibre, member, angle, station]
            sign = np.sign(stress)
            factor = safety.static / material.yield_strength
            grad_axial[angle, station] = sign * factor / area[member]
            grad_bending[angle, station] = sign * factor * signs[fibre] / modulus[member]
            slope = (
                -sign
                * factor
                * (
                    2.0 * axial[member, angle, station] / area[member]
                    + 3.0 * signs[fibre] * bending[member, angle, station] / modulus[member]
                )
                / diameter
            )

        elif mode == 1:  # fatigue
            highest = np.max(fibres[:, member], axis=1)
            lowest = np.min(fibres[:, member], axis=1)
            alternating = 0.5 * (highest - lowest)
            mean = 0.5 * (highest + lowest)
            endurance = float(material.endurance_limit(diameter)[0])
            utilisation = (
                alternating / endurance + np.maximum(mean, 0.0) / material.ultimate_strength
            )
            fibre, station = (
                int(v) for v in np.unravel_index(int(np.argmax(utilisation)), utilisation.shape)
            )
            high_angle: int = int(np.argmax(fibres[fibre, member, :, station]))
            low_angle: int = int(np.argmin(fibres[fibre, member, :, station]))
            positive_mean = mean[fibre, station] > 0.0
            weight_high = safety.fatigue * (
                0.5 / endurance + (0.5 / material.ultimate_strength if positive_mean else 0.0)
            )
            weight_low = safety.fatigue * (
                -0.5 / endurance + (0.5 / material.ultimate_strength if positive_mean else 0.0)
            )
            for angle, weight in ((high_angle, weight_high), (low_angle, weight_low)):
                grad_axial[angle, station] += weight / area[member]
                grad_bending[angle, station] += weight * signs[fibre] / modulus[member]

            # Bound as defaults so the closure captures this member's indices,
            # not whatever the loop holds when it is finally called.
            def stress_slope(
                angle: int,
                _member: int = member,
                _station: int = int(station),
                _fibre: int = int(fibre),
                _diameter: float = float(diameter),
            ) -> float:
                """``d sigma / dd``, from ``1/A ~ d^-2`` and ``1/Z ~ d^-3``."""
                return (
                    -(
                        2.0 * axial[_member, angle, _station] / area[_member]
                        + 3.0
                        * signs[_fibre]
                        * bending[_member, angle, _station]
                        / modulus[_member]
                    )
                    / _diameter
                )

            slope_high = stress_slope(high_angle)
            slope_low = stress_slope(low_angle)
            d_alternating = 0.5 * (slope_high - slope_low)
            d_mean = 0.5 * (slope_high + slope_low)
            exponent = -0.107 if diameter <= 51.0 else -0.157
            d_endurance = endurance * exponent / diameter
            slope = safety.fatigue * (
                d_alternating / endurance
                - alternating[fibre, station] * d_endurance / endurance**2
                + (d_mean / material.ultimate_strength if positive_mean else 0.0)
            )

        else:  # Euler buckling
            flat = int(np.argmin(axial[member]))
            angle, station = (int(v) for v in np.unravel_index(flat, axial[member].shape))
            second_moment = np.pi * diameter**4 / 64.0
            critical = (
                np.pi**2
                * material.youngs_modulus
                * second_moment
                / (MEMBER_FIXITY[member] * lengths[member]) ** 2
            )
            grad_axial[angle, station] = -safety.buckling / critical
            slope = -4.0 * buckling[member] / diameter
            grad_length = 2.0 * buckling[member] / lengths[member]

        if slope == 0.0:  # pragma: no cover - a mode with no sensitivity to d
            continue
        d_axial[member] = -grad_axial / slope
        d_bending[member] = -grad_bending / slope
        d_length[member] = -grad_length / slope

    return SizingJacobian(d_axial=d_axial, d_bending=d_bending, d_length=d_length)


@dataclass(frozen=True)
class CoupledJacobian:
    """Everything the two coupled disciplines need to report a local Jacobian."""

    axial: FloatArray
    """``d(axial)/dp``, ``(n_members, n_angles, stations, N_PARAMETERS)``."""

    bending: FloatArray
    """``d(bending)/dp``, same shape."""

    mean_torque: FloatArray
    """``d(mean M_r)/dp``, ``(N_PARAMETERS,)``."""

    peak_bearing_load: FloatArray
    """``d(max |R_1|)/dp``, ``(N_PARAMETERS,)``."""


def coupled_jacobian(
    design: Design,
    analysis: object,
    diameters: dict[str, float],
    properties: MassProperties,
    loads: DynamicLoads,
    stations: int,
    material: Material,
    spec: EngineSpec = DEFAULT_SPEC,
) -> CoupledJacobian:
    """Assemble the dynamics discipline's local Jacobian, end to end.

    Chains :func:`exlink.jacobian.kinematic_jacobian` through the mass
    properties, the accelerations, the 18x18 solve and the internal loads.

    Args:
        design: The mechanism dimensions.
        analysis: A valid :class:`~exlink.model.Analysis`.
        diameters: Section diameter of each member [mm].
        properties: The mass properties used for the load case.
        loads: The solved load case.
        stations: Sections evaluated along each member.
        material: The material, for its density.
        spec: Fixed engine data.

    Returns:
        The derivatives the discipline reports.
    """
    from .jacobian import (
        _at_extremum,
        _refined_extremum,
        gas_force_jacobian,
        kinematic_jacobian,
    )

    solved = analysis.require_solved()  # type: ignore[attr-defined]
    kinematics = solved.kinematics
    phases = solved.thermodynamics.phases

    kinematic = kinematic_jacobian(design, kinematics, spec)
    mass = mass_property_jacobian(
        design, kinematics, kinematic, properties, diameters, material.density
    )
    acceleration = acceleration_jacobian(kinematic, mass, loads.speed, kinematics.theta_1.size)

    tops = phases.maxima_indices
    deep, shallow = phases.minima_indices
    offsets = {i: _refined_extremum(kinematics.lam, i)[0] for i in (*tops, deep, shallow)}
    values = {i: _refined_extremum(kinematics.lam, i)[1] for i in (*tops, deep, shallow)}
    top = tops[0] if values[tops[0]] >= values[tops[1]] else tops[1]
    d_lam_top = _at_extremum(kinematic.lam, top, offsets[top])
    d_compression = d_lam_top - _at_extremum(kinematic.lam, shallow, offsets[shallow])
    gas = gas_force_jacobian(analysis, kinematic, d_lam_top, d_compression, spec)

    reaction = equilibrium_jacobian(
        design, kinematics, kinematic, properties, mass, acceleration, loads, gas, spec
    )
    axial, bending = member_load_jacobian(
        design,
        kinematics,
        kinematic,
        properties,
        mass,
        acceleration,
        loads,
        reaction,
        stations,
    )

    mean_torque = reaction[:, _I_TORQUE, :].mean(axis=0)

    first, second = _IDX["R1"]
    bearing = loads.reaction["R1"]
    magnitude = np.linalg.norm(bearing, axis=1)
    peak = int(np.argmax(magnitude))
    peak_bearing = (
        bearing[peak, 0] * reaction[peak, first, :]
        + bearing[peak, 1] * reaction[peak, second, :]
    ) / magnitude[peak]

    return CoupledJacobian(
        axial=axial,
        bending=bending,
        mean_torque=mean_torque,
        peak_bearing_load=peak_bearing,
    )
