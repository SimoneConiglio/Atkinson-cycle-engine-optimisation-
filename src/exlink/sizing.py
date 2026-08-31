"""Sizing every member against static yield, fatigue and buckling.

Each link is taken as a solid round bar and its diameter is solved for: the
smallest section that satisfies all three failure modes over the whole
revolution.

Internal loads
--------------
For a member spanning two joints, the internal force and moment at a section a
fraction ``s`` along it follow from the free body of the piece ``[0, s]``: the
force applied at its first end, plus its own distributed d'Alembert load.
Because a rigid body's acceleration varies *linearly* along any straight line
through it, that load is linear in ``s`` and both integrals are closed form ---
see :func:`internal_loads`.

The trigonal link is the one place an idealisation is needed.  It is a single
rigid part, so as a frame it is three times statically indeterminate.  It is
treated here as a pin-jointed triangle: each side carries an axial force from
joint equilibrium (which *is* determinate), and bends only under its own
distributed inertia as a simply supported beam.  That is the standard first-cut
decomposition; it captures the axial load exactly and under-estimates bending
at the corners, which is noted rather than hidden.

Failure modes
-------------
``static``
    Peak absolute fibre stress against ``S_y / n_y``.  The state is uniaxial
    (axial plus bending), so von Mises reduces to ``|sigma|``.
``fatigue``
    Goodman, evaluated per extreme fibre so that the alternating and mean
    components are taken at a fixed material point rather than at whichever
    fibre happens to be worst at each instant.
``buckling``
    Euler, for the peak compressive axial load, with an end-fixity factor of 1
    for pin-ended links and 2 for a crank throw cantilevered off its shaft.

All three utilisations fall as the diameter grows, so the required diameter is
found by bisection --- robust, derivative-free, and needed anyway because the
fatigue size factor ``k_b`` itself depends on the diameter being solved for.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec
from .cycle import Thermodynamics
from .dynamics import MEMBER_NAMES, MEMBERS, DynamicLoads
from .materials import (
    DEFAULT_MATERIAL,
    DEFAULT_SAFETY,
    Material,
    SafetyFactors,
    goodman_utilisation,
)

FloatArray = NDArray[np.float64]

STATIONS = 9
"""Sections evaluated along each member, ends included."""

MIN_DIAMETER = 2.0
"""Smallest diameter considered [mm]; below this a link is not manufacturable."""

MAX_DIAMETER = 250.0
"""Largest diameter considered [mm]; hitting it means the member cannot be sized."""

BISECTION_STEPS = 60
"""Bisection iterations; 60 halvings of the bracket reach machine precision."""

#: End-fixity factor ``K`` in the Euler load, by member kind.
END_FIXITY: dict[str, float] = {"link": 1.0, "truss": 1.0, "cantilever": 2.0}


def internal_loads(
    start: FloatArray,
    end: FloatArray,
    end_force: FloatArray,
    mass: float,
    start_acceleration: FloatArray,
    end_acceleration: FloatArray,
    stations: int = STATIONS,
) -> tuple[FloatArray, FloatArray]:
    """Axial force and bending moment along a uniform straight member.

    Equilibrium of the piece between the first end and a section at fraction
    ``s``, carrying the applied end force and its own share of the distributed
    inertia load.  With ``a(s) = a_1 + s (a_2 - a_1)`` and ``r(s)`` linear, the
    two integrals evaluate in closed form to

    .. code-block:: text

        F(s) = m [a_1 s + (a_2 - a_1) s^2 / 2] - F_1
        M(s) = -m [(dr x a_1)_z s^2 / 2 + (dr x da)_z s^3 / 6] + s (dr x F_1)_z

    Args:
        start: First-end positions, ``(n_angles, 2)`` [mm].
        end: Second-end positions, ``(n_angles, 2)`` [mm].
        end_force: Force applied to the member at its first end, ``(n_angles, 2)`` [N].
        mass: Member mass [tonne].
        start_acceleration: Acceleration of the first end, ``(n_angles, 2)`` [mm/s^2].
        end_acceleration: Acceleration of the second end [mm/s^2].
        stations: Number of sections, spread evenly including both ends.

    Returns:
        ``(axial, bending)``, each ``(n_angles, stations)``, in N and N.mm.
        Axial is positive in tension.
    """
    s = np.linspace(0.0, 1.0, stations)[None, :]  # (1, stations)
    delta_r = (end - start)[:, None, :]  # (n, 1, 2)
    length = np.linalg.norm(delta_r, axis=-1)
    direction = delta_r / np.where(length[..., None] == 0.0, 1.0, length[..., None])

    a_1 = start_acceleration[:, None, :]
    delta_a = (end_acceleration - start_acceleration)[:, None, :]
    force_1 = end_force[:, None, :]

    inertia_force = mass * (a_1 * s[..., None] + delta_a * s[..., None] ** 2 / 2.0)
    internal_force = inertia_force - force_1

    def cross(u: FloatArray, v: FloatArray) -> FloatArray:
        return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]

    bending = -mass * (
        cross(delta_r, a_1) * s**2 / 2.0 + cross(delta_r, delta_a) * s**3 / 6.0
    ) + s * cross(delta_r, force_1)
    axial = np.sum(internal_force * direction, axis=-1)
    return axial, bending


def _simply_supported_reaction(
    normal: FloatArray,
    mass: float,
    start_acceleration: FloatArray,
    end_acceleration: FloatArray,
) -> FloatArray:
    """Transverse end reaction of a pin-ended member under its own inertia.

    For a distributed load linear along the span, the reaction at the first end
    is ``-m (a_1n / 2 + (a_2n - a_1n) / 6)`` in the transverse direction.
    """
    a_1n = np.sum(start_acceleration * normal, axis=-1)
    a_2n = np.sum(end_acceleration * normal, axis=-1)
    return -mass * (a_1n / 2.0 + (a_2n - a_1n) / 6.0)


def _truss_axial_forces(loads: DynamicLoads) -> dict[str, FloatArray]:
    """Axial forces in the three sides of the trigonal link.

    Joint equilibrium of a pin-jointed triangle: six scalar equations in three
    unknowns, consistent because the vertex loads already balance the link's
    own inertia.  Solved in the least-squares sense, which is exact here.
    """
    kinematics = loads.kinematics
    properties = loads.mass_properties
    vertices = {"A": kinematics.A, "D": kinematics.D, "E": kinematics.E}
    sides = [m for m in MEMBERS if m.kind == "truss"]

    external = {
        "A": loads.reaction["A"],
        "D": loads.reaction["D"],
        "E": -loads.reaction["E"],
    }
    # Half of each side's inertia force is lumped at each of its two ends; the
    # halves sum to the side's true resultant because acceleration varies
    # linearly along it.
    lumped = {name: np.zeros(0) for name in vertices}
    for name, point in vertices.items():
        del point
        attached = sum(
            properties.member_mass[m.name] for m in sides if name in (m.start, m.end)
        )
        lumped[name] = -0.5 * attached * loads.joint_acceleration[name]

    n = kinematics.theta_1.size
    matrix = np.zeros((n, 6, 3))
    rhs = np.zeros((n, 6))
    for row, name in enumerate(vertices):
        for column, side in enumerate(sides):
            if name not in (side.start, side.end):
                continue
            other = side.end if name == side.start else side.start
            delta = vertices[other] - vertices[name]
            unit = delta / np.linalg.norm(delta, axis=-1, keepdims=True)
            matrix[:, 2 * row, column] = unit[:, 0]
            matrix[:, 2 * row + 1, column] = unit[:, 1]
        rhs[:, 2 * row] = -(external[name][:, 0] + lumped[name][:, 0])
        rhs[:, 2 * row + 1] = -(external[name][:, 1] + lumped[name][:, 1])

    normal = np.einsum("nij,nik->njk", matrix, matrix)
    projected = np.einsum("nij,ni->nj", matrix, rhs)
    solution = np.linalg.solve(normal, projected[..., None])[..., 0]
    return {side.name: solution[:, column] for column, side in enumerate(sides)}


def member_loads(
    loads: DynamicLoads, stations: int = STATIONS
) -> dict[str, tuple[FloatArray, FloatArray]]:
    """Internal axial force and bending moment for every sized member.

    Args:
        loads: A solved dynamic load case.
        stations: Sections evaluated along each member.

    Returns:
        ``{member name: (axial, bending)}``, each array ``(n_angles, stations)``.
    """
    kinematics = loads.kinematics
    properties = loads.mass_properties
    joints = {name: getattr(kinematics, name) for name in ("R1", "Q", "A", "D", "R2", "E", "P")}
    truss_axial = _truss_axial_forces(loads)

    # Force applied to each member at its first end.  Reactions are stored
    # pointing forward along the chain, so a body on the crankshaft side of a
    # joint feels the negative.
    applied: dict[str, FloatArray] = {
        "crank_1": -loads.reaction["Q"],
        "swing_rod": loads.reaction["Q"],
        "crank_2": -loads.reaction["D"],
        "piston_rod": loads.reaction["E"],
    }

    result: dict[str, tuple[FloatArray, FloatArray]] = {}
    for member in MEMBERS:
        start, end = joints[member.start], joints[member.end]
        a_start = loads.joint_acceleration[member.start]
        a_end = loads.joint_acceleration[member.end]
        mass = properties.member_mass[member.name]

        if member.kind == "truss":
            delta = end - start
            unit = delta / np.linalg.norm(delta, axis=-1, keepdims=True)
            normal = np.stack([-unit[:, 1], unit[:, 0]], axis=-1)
            axial_force = truss_axial[member.name]
            shear = _simply_supported_reaction(normal, mass, a_start, a_end)
            end_force = -axial_force[:, None] * unit + shear[:, None] * normal
        else:
            end_force = applied[member.name]

        result[member.name] = internal_loads(
            start, end, end_force, mass, a_start, a_end, stations=stations
        )
    return result


@dataclass(frozen=True)
class MemberSizing:
    """Sizing outcome for one member."""

    name: str
    diameter: float
    """Required section diameter [mm]."""

    static_utilisation: float
    """Peak stress over ``S_y / n_y``; ``<= 1`` is safe."""

    fatigue_utilisation: float
    """Goodman utilisation times ``n_f``; ``<= 1`` is safe."""

    buckling_utilisation: float
    """Peak compression over the Euler load divided by ``n_b``; ``<= 1`` is safe."""

    peak_stress: float
    """Largest absolute fibre stress at the sized diameter [MPa]."""

    mass: float
    """Member mass at the sized diameter [tonne]."""

    critical_mode: str
    """Which of the three modes set the diameter."""

    @property
    def mass_kg(self) -> float:
        """Member mass [kg]."""
        return 1000.0 * self.mass


def _utilisations(
    diameter: FloatArray,
    axial: FloatArray,
    bending: FloatArray,
    lengths: FloatArray,
    fixity: FloatArray,
    material: Material,
    safety: SafetyFactors,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Static, fatigue and buckling utilisations for a vector of diameters.

    Args:
        diameter: One diameter per member, ``(n_members,)`` [mm].
        axial: ``(n_members, n_angles, n_stations)`` [N], tension positive.
        bending: Same shape [N.mm].
        lengths: Member lengths ``(n_members,)`` [mm].
        fixity: Euler end-fixity factor per member ``(n_members,)``.
        material: The material.
        safety: The design factors.

    Returns:
        ``(static, fatigue, buckling)``, each ``(n_members,)`` and ``<= 1`` when safe.
    """
    area = np.pi * diameter**2 / 4.0
    modulus = np.pi * diameter**3 / 32.0
    second_moment = np.pi * diameter**4 / 64.0

    direct = axial / area[:, None, None]
    flexural = bending / modulus[:, None, None]
    # Both extreme fibres, so that the alternating and mean components below are
    # taken at a fixed material point.
    fibres = np.stack([direct + flexural, direct - flexural], axis=0)

    peak = np.max(np.abs(fibres), axis=(0, 2, 3))
    static = peak * safety.static / material.yield_strength

    highest = np.max(fibres, axis=2)
    lowest = np.min(fibres, axis=2)
    alternating = 0.5 * (highest - lowest)
    mean = 0.5 * (highest + lowest)
    # ``alternating`` is (fibre, member, station); the endurance limit varies
    # only with the member's diameter.
    endurance = material.endurance_limit(diameter)[None, :, None]
    utilisation = goodman_utilisation(alternating, mean, endurance, material.ultimate_strength)
    fatigue = np.max(utilisation, axis=(0, 2)) * safety.fatigue

    compression = np.maximum(-np.min(axial, axis=(1, 2)), 0.0)
    euler = np.pi**2 * material.youngs_modulus * second_moment / (fixity * lengths) ** 2
    buckling = compression * safety.buckling / euler
    return static, fatigue, buckling


def member_lengths(design: object) -> FloatArray:
    """Length of each sized member, in the order of :data:`MEMBERS` [mm]."""
    return np.array([abs(float(getattr(design, m.length_attribute))) for m in MEMBERS])


#: Euler end-fixity factor of each member, in the order of :data:`MEMBERS`.
MEMBER_FIXITY: FloatArray = np.array([END_FIXITY[m.kind] for m in MEMBERS])

#: Whether a member's slenderness is worth policing, in the order of :data:`MEMBERS`.
#:
#: Only the connecting links. A crank *throw* is routinely shorter than the pin
#: diameter it carries -- 8 mm of throw on a 19 mm crankpin is an ordinary
#: engine, not a modelling failure -- so applying a slenderness limit to it
#: would reject sound designs. The beam model over-predicts stress for such a
#: stubby part, which errs the safe way.
MEMBER_IS_SLENDER: NDArray[np.bool_] = np.array([m.kind != "cantilever" for m in MEMBERS])


def size_from_arrays(
    axial: FloatArray,
    bending: FloatArray,
    lengths: FloatArray,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    fixity: FloatArray | None = None,
    names: Sequence[str] | None = None,
) -> dict[str, MemberSizing]:
    """Solve for the smallest safe diameter, from raw internal-load arrays.

    Split out from :func:`size_members` so that the GEMSEO structural
    discipline can be driven by the coupling variables directly, without having
    to reconstruct a :class:`~exlink.dynamics.DynamicLoads` -- and so that
    :mod:`exlink.slidercrank` can size a *different* mechanism through exactly
    the same yield, fatigue and buckling checks.  Comparing two mechanisms is
    only meaningful if the structural model is literally the same code, so the
    member list is a parameter rather than a constant.

    Args:
        axial: ``(n_members, n_angles, n_stations)`` [N], tension positive.
        bending: Same shape [N.mm].
        lengths: Member lengths ``(n_members,)`` [mm].
        material: The material.
        safety: The design factors.
        fixity: Euler end-fixity factor per member; the EX-link's if omitted.
        names: Member names; the EX-link's if omitted.

    Returns:
        ``{member name: MemberSizing}``.
    """
    fixity = MEMBER_FIXITY if fixity is None else np.asarray(fixity, dtype=float)
    labels = tuple(MEMBER_NAMES) if names is None else tuple(names)
    count = len(labels)

    low = np.full(count, MIN_DIAMETER)
    high = np.full(count, MAX_DIAMETER)
    for _ in range(BISECTION_STEPS):
        middle = 0.5 * (low + high)
        static, fatigue, buckling = _utilisations(
            middle, axial, bending, lengths, fixity, material, safety
        )
        safe = np.maximum(np.maximum(static, fatigue), buckling) <= 1.0
        high = np.where(safe, middle, high)
        low = np.where(safe, low, middle)

    diameter = high
    static, fatigue, buckling = _utilisations(
        diameter, axial, bending, lengths, fixity, material, safety
    )
    area = np.pi * diameter**2 / 4.0
    modulus = np.pi * diameter**3 / 32.0
    peak = np.max(
        np.abs(axial / area[:, None, None]) + np.abs(bending / modulus[:, None, None]),
        axis=(1, 2),
    )
    mass = material.density * area * lengths

    modes = np.stack([static, fatigue, buckling])
    mode_names = ("static", "fatigue", "buckling")
    return {
        label: MemberSizing(
            name=label,
            diameter=float(diameter[i]),
            static_utilisation=float(static[i]),
            fatigue_utilisation=float(fatigue[i]),
            buckling_utilisation=float(buckling[i]),
            peak_stress=float(peak[i]),
            mass=float(mass[i]),
            critical_mode=mode_names[int(np.argmax(modes[:, i]))],
        )
        for i, label in enumerate(labels)
    }


def size_members(
    loads: DynamicLoads,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    stations: int = STATIONS,
) -> dict[str, MemberSizing]:
    """Solve for the smallest safe diameter of every member.

    Args:
        loads: A solved dynamic load case.
        material: The material.
        safety: The design factors.
        stations: Sections evaluated along each member.

    Returns:
        ``{member name: MemberSizing}``.
    """
    per_member = member_loads(loads, stations=stations)
    return size_from_arrays(
        np.stack([per_member[name][0] for name in MEMBER_NAMES]),
        np.stack([per_member[name][1] for name in MEMBER_NAMES]),
        member_lengths(loads.kinematics.design),
        material,
        safety,
    )


def piston_mass(
    thermodynamics: Thermodynamics,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    spec: EngineSpec = DEFAULT_SPEC,
    skirt_thickness: float = 2.5,
) -> tuple[float, float]:
    """Size the piston crown against peak gas pressure and return its mass.

    The crown is treated as a clamped circular plate of radius ``phi / 2`` under
    uniform pressure, whose peak stress is ``3 p r^2 / (4 t^2)``; the skirt is a
    thin cylinder of the given wall thickness.

    The piston sits *outside* the sizing loop: its thickness follows from the
    gas pressure alone, which no structural choice changes.  Its mass still
    matters, because it is the largest reciprocating mass in the mechanism.

    Args:
        thermodynamics: The solved cycle, for the peak gauge pressure.
        material: The material.
        safety: The design factors.
        spec: Fixed engine data.
        skirt_thickness: Wall thickness of the skirt [mm].

    Returns:
        ``(crown_thickness_mm, mass_tonne)``.
    """
    return piston_mass_from_pressure(
        float(np.max(thermodynamics.gauge_pressure)), material, safety, spec, skirt_thickness
    )


def piston_mass_from_pressure(
    peak_pressure: float,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    spec: EngineSpec = DEFAULT_SPEC,
    skirt_thickness: float = 2.5,
) -> tuple[float, float]:
    """Size the piston from a peak gauge pressure alone.

    The cycle-model-agnostic half of :func:`piston_mass`, so that an Otto
    slider-crank and an Atkinson linkage get identically-sized pistons for the
    same peak pressure.

    Args:
        peak_pressure: Peak in-cylinder gauge pressure [MPa].
        material: The material.
        safety: The design factors.
        spec: Fixed engine data.
        skirt_thickness: Wall thickness of the skirt [mm].

    Returns:
        ``(crown_thickness_mm, mass_tonne)``.
    """
    allowable = material.yield_strength / safety.static
    radius = 0.5 * spec.bore
    crown = radius * np.sqrt(3.0 * peak_pressure / (4.0 * allowable))
    crown = max(float(crown), 2.0)

    crown_volume = np.pi * radius**2 * crown
    skirt_volume = (
        np.pi * ((radius) ** 2 - (radius - skirt_thickness) ** 2) * spec.piston_length
    )
    return crown, float(material.density * (crown_volume + skirt_volume))
