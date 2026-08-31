"""The gear pair, and the discrete variable hiding inside the inter-axle distance.

The linkage carries one 2:1 external gear pair between the crankshaft and the
eccentric shaft.  In the geometric problem it appears only as a pair of
primitive radii, ``r_1 = 2I/3`` and ``r_2 = I/3``, which are continuous
functions of the design variable ``I``.  That is a fiction, and it is the
interesting kind.

Why ``I`` is not really continuous
-----------------------------------
A gear has an integer number of teeth cut with a standard-module hob.  For a
pair on a fixed centre distance,

.. math:: r = \\frac{m z}{2}, \\qquad I = r_1 + r_2 = \\frac{m (z_1 + z_2)}{2}

and the 2:1 ratio forces ``z_1 = 2 z_2``.  Together:

.. math:: I = \\tfrac{3}{2} m z_2, \\qquad m \\in \\text{ISO 54}, \\; z_2 \\in \\mathbb{Z}

So the inter-axle distance lives on a **lattice**, not on an interval.  At
``m = 1.5`` the spacing is 2.25 mm; at ``m = 3`` it is 4.5 mm.  A continuous
optimizer that returns ``I = 61.37 mm`` has returned a mechanism that cannot be
geared, and the nearest buildable neighbours are 60.75 and 63.00.

This matters more than it sounds, because ``I`` is one of the variables the
equality constraints ``STE = 74`` and ``epsilon = 16`` are satisfied *with*.
Snapping ``I`` to the lattice moves the design off both equalities, and the
remaining continuous variables have to be re-solved to restore them.  That is
the classic mixed-integer structure: choose the discrete variable, then repair
the continuous ones.  :func:`lattice_neighbours` enumerates the candidates and
:mod:`exlink.scenarios` does the repair.

Undercutting sets the floor.  A 20 deg involute pinion needs
``z >= 2 / sin^2(alpha) = 17`` teeth to avoid an undercut root, so the smaller
gear pins ``z_2 >= 17`` and hence ``I >= 25.5 m``.

Sizing the teeth
----------------
Two failure modes, both standard:

* **Root bending** (Lewis): ``sigma = F_t / (b m Y)``, with ``Y`` the form
  factor of the tooth count.
* **Flank contact** (Hertz line contact, the AGMA pitting check):
  ``sigma_H = Z_E sqrt(F_t (u + 1) / (b d_1 u))``.

Both are inversely proportional to face width, so ``b`` is what carries the
load once ``m`` is chosen.  The face width is bounded above at ``12 m``: wider
than that and the tooth load stops being uniform across the face, because no
real shaft is stiff enough to keep the mesh aligned.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .manufacturing import (
    MIN_FACE_WIDTH,
    MIN_TEETH,
    STANDARD_MODULES,
    round_to_module,
)
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors

ELASTIC_COEFFICIENT = 191.0
"""``Z_E`` for a steel-on-steel pair [sqrt(MPa)].

``Z_E = sqrt(1 / (pi ((1 - nu_1^2)/E_1 + (1 - nu_2^2)/E_2)))`` at
``E = 210 GPa``, ``nu = 0.3``.
"""

CONTACT_ALLOWABLE = 1100.0
"""Allowable flank contact stress for a case-hardened steel gear [MPa]."""

MAX_WIDTH_FACTOR = 12.0
"""Largest face width as a multiple of the module, for mesh alignment."""

WEB_FRACTION = 0.35
"""Web thickness as a fraction of face width, for the blank mass model."""


def lewis_form_factor(teeth: float) -> float:
    """Lewis form factor ``Y`` for a 20 deg full-depth involute tooth.

    A smooth fit to the standard table over ``z in [17, 200]``, accurate to
    about 2 % across that range, which is finer than the uncertainty in the
    allowable stress it multiplies.

    Args:
        teeth: Number of teeth.

    Returns:
        ``Y``, dimensionless.
    """
    z = max(float(teeth), float(MIN_TEETH))
    return 0.484 - 2.87 / z


def tooth_count(inter_axle: float, module: float) -> int:
    """Teeth on the smaller gear for a given centre distance and module.

    Args:
        inter_axle: ``I``, the centre distance [mm].
        module: Gear module [mm].

    Returns:
        ``z_2 = 2 I / (3 m)``, rounded to the nearest integer and floored at
        :data:`~exlink.manufacturing.MIN_TEETH`.
    """
    return max(round(2.0 * float(inter_axle) / (3.0 * float(module))), MIN_TEETH)


def lattice_inter_axle(module: float, teeth: int) -> float:
    """The centre distance a module and tooth count actually produce [mm].

    Args:
        module: Gear module [mm].
        teeth: Teeth on the smaller gear.

    Returns:
        ``I = 1.5 m z_2``.
    """
    return 1.5 * float(module) * float(teeth)


def lattice_neighbours(
    inter_axle: float,
    modules: np.ndarray = STANDARD_MODULES,
    count: int = 2,
) -> list[tuple[float, int, float]]:
    """Buildable centre distances bracketing a continuous one.

    For every standard module, the tooth counts whose centre distance sits
    closest to the requested one are returned.  This is the candidate list a
    mixed-integer search enumerates.

    Args:
        inter_axle: The continuous ``I`` an optimizer asked for [mm].
        modules: Standard modules to consider [mm].
        count: Tooth counts to keep on each side, per module.

    Returns:
        ``[(module, teeth, inter_axle), ...]`` sorted by distance from the
        requested value.  Entries below the undercut limit are dropped.
    """
    candidates: list[tuple[float, int, float]] = []
    for module in np.atleast_1d(np.asarray(modules, dtype=float)):
        centre = tooth_count(inter_axle, float(module))
        for teeth in range(centre - count, centre + count + 1):
            if teeth < MIN_TEETH:
                continue
            candidates.append(
                (float(module), int(teeth), lattice_inter_axle(float(module), teeth))
            )
    unique = {(m, z): value for m, z, value in candidates}
    ordered = sorted(
        ((m, z, value) for (m, z), value in unique.items()),
        key=lambda item: abs(item[2] - float(inter_axle)),
    )
    return ordered


@dataclass(frozen=True)
class GearPair:
    """A sized 2:1 external gear pair."""

    module: float
    """Standard module ``m`` [mm]."""

    teeth_small: int
    """``z_2``, on the eccentric shaft."""

    teeth_large: int
    """``z_1 = 2 z_2``, on the crankshaft."""

    face_width: float
    """``b``, common to both gears [mm]."""

    inter_axle: float
    """``I = 1.5 m z_2`` [mm]."""

    bending_utilisation: float
    """Lewis root stress over its allowable; ``<= 1`` is safe."""

    contact_utilisation: float
    """Hertz flank stress over its allowable; ``<= 1`` is safe."""

    mass: float
    """Mass of both blanks [tonne]."""

    inertia_small: float
    """Polar inertia of the small gear about its own axis [tonne mm^2]."""

    inertia_large: float
    """Polar inertia of the large gear [tonne mm^2]."""

    @property
    def mass_kg(self) -> float:
        """Mass of both blanks [kg]."""
        return 1000.0 * self.mass

    @property
    def radius_small(self) -> float:
        """Primitive radius of the small gear [mm]."""
        return 0.5 * self.module * self.teeth_small

    @property
    def radius_large(self) -> float:
        """Primitive radius of the large gear [mm]."""
        return 0.5 * self.module * self.teeth_large

    @property
    def width_factor(self) -> float:
        """Face width as a multiple of the module."""
        return self.face_width / self.module

    @property
    def feasible(self) -> bool:
        """Whether the pair is both safe and within the face-width limit."""
        return (
            max(self.bending_utilisation, self.contact_utilisation) <= 1.0
            and self.width_factor <= MAX_WIDTH_FACTOR + 1.0e-9
            and self.teeth_small >= MIN_TEETH
        )


def _blank_mass_and_inertia(
    primitive: float,
    module: float,
    face: float,
    bore: float,
    density: float,
) -> tuple[float, float]:
    """Mass and polar inertia of one gear blank.

    A rim of full face width carrying the teeth, a thinner web, and a hub
    around the shaft.  Crude, but it captures the two things that matter: mass
    grows with ``r^2 b``, and inertia with ``r^4 b``, so a large-module gear on
    a wide face is expensive twice over.

    Args:
        primitive: Primitive radius [mm].
        module: Module [mm].
        face: Face width [mm].
        bore: Shaft bore diameter [mm].
        density: Material density [tonne/mm^3].

    Returns:
        ``(mass_tonne, polar_inertia)``.
    """
    outer = primitive + module
    root = max(primitive - 1.25 * module, 0.55 * outer)
    hub_outer = max(0.9 * bore, 0.5 * bore + 4.0)
    hub_inner = 0.5 * bore

    def annulus(r_outer: float, r_inner: float, width: float) -> tuple[float, float]:
        r_outer = max(r_outer, r_inner)
        volume = math.pi * (r_outer**2 - r_inner**2) * width
        mass = density * volume
        # Polar inertia of a hollow cylinder: m (r_o^2 + r_i^2) / 2.
        return mass, 0.5 * mass * (r_outer**2 + r_inner**2)

    rim = annulus(outer, root, face)
    web = annulus(root, hub_outer, WEB_FRACTION * face)
    hub = annulus(hub_outer, hub_inner, face)
    mass = rim[0] + web[0] + hub[0]
    inertia = rim[1] + web[1] + hub[1]
    return mass, inertia


def size_pair(
    inter_axle: float,
    tangential_force: float,
    module: float | None = None,
    teeth: int | None = None,
    shaft_bore: float = 12.0,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
) -> GearPair:
    """Size a 2:1 pair for a peak tangential tooth load.

    The module and tooth count may be given (a mixed-integer search fixes them
    and asks what it costs) or left out, in which case the smallest standard
    module that can carry the load within the face-width limit is chosen.

    Args:
        inter_axle: Requested centre distance ``I`` [mm].  Used to pick the
            tooth count when ``teeth`` is not given; the *realised* centre
            distance is on the lattice and is returned in the result.
        tangential_force: Peak tangential tooth load ``F_t`` [N].
        module: Standard module to use [mm]; chosen automatically if omitted.
        teeth: Teeth on the small gear; derived from ``inter_axle`` if omitted.
        shaft_bore: Shaft diameter the gears are bored for [mm].
        material: Gear material.
        safety: Design factors; the static factor is applied to root bending.

    Returns:
        The sized pair.  Check :attr:`GearPair.feasible`.
    """
    force = max(abs(float(tangential_force)), 1.0e-9)
    bending_allowable = material.yield_strength / safety.static

    def build(module_value: float, teeth_value: int) -> GearPair:
        z_small = max(int(teeth_value), MIN_TEETH)
        z_large = 2 * z_small
        # Both checks are on the small gear: it sees the same tooth load, has
        # the weaker form factor and turns more often.
        form = lewis_form_factor(z_small)
        face_bending = force / (module_value * form * bending_allowable)

        pitch_diameter = module_value * z_small
        ratio = 2.0
        face_contact = (
            force
            * (ratio + 1.0)
            * ELASTIC_COEFFICIENT**2
            / (pitch_diameter * ratio * CONTACT_ALLOWABLE**2)
        )
        face = max(face_bending, face_contact, MIN_FACE_WIDTH)

        sigma_bending = force / (face * module_value * form)
        sigma_contact = ELASTIC_COEFFICIENT * math.sqrt(
            force * (ratio + 1.0) / (face * pitch_diameter * ratio)
        )

        small = _blank_mass_and_inertia(
            0.5 * module_value * z_small, module_value, face, shaft_bore, material.density
        )
        large = _blank_mass_and_inertia(
            0.5 * module_value * z_large, module_value, face, shaft_bore, material.density
        )
        return GearPair(
            module=module_value,
            teeth_small=z_small,
            teeth_large=z_large,
            face_width=face,
            inter_axle=lattice_inter_axle(module_value, z_small),
            bending_utilisation=sigma_bending / bending_allowable,
            contact_utilisation=sigma_contact / CONTACT_ALLOWABLE,
            mass=small[0] + large[0],
            inertia_small=small[1],
            inertia_large=large[1],
        )

    if module is not None:
        chosen = round_to_module(module)
        return build(chosen, teeth if teeth is not None else tooth_count(inter_axle, chosen))

    # No module given: take the lightest standard one that stays inside the
    # face-width limit.  Small modules need less blank material but a wider
    # face, so the limit is what decides.
    best: GearPair | None = None
    for candidate in STANDARD_MODULES:
        pair = build(float(candidate), teeth or tooth_count(inter_axle, float(candidate)))
        if not pair.feasible:
            continue
        if best is None or pair.mass < best.mass:
            best = pair
    if best is not None:
        return best
    return build(
        float(STANDARD_MODULES[-1]),
        teeth or tooth_count(inter_axle, float(STANDARD_MODULES[-1])),
    )


def buildable_neighbours(
    inter_axle: float,
    tangential_force: float,
    shaft_bore: float = 12.0,
    count: int = 3,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
) -> list[tuple[float, int, float, GearPair]]:
    """Lattice points that can actually carry the tooth load, nearest first.

    :func:`lattice_neighbours` ranks purely by distance from the requested
    centre distance, and that is the wrong order to optimise in.  The nearest
    lattice points are reached with the *smallest* modules, and a small module
    needs a wide face to carry a given tooth load -- so the geometrically
    closest candidates are routinely the structurally worst, and can be
    unbuildable outright.

    Ranking instead by whether the pair is within its face-width limit, then by
    distance, is what an enumeration over the discrete variable should do.

    Args:
        inter_axle: The continuous ``I`` an optimizer asked for [mm].
        tangential_force: Peak tangential tooth load ``F_t`` [N].
        shaft_bore: Shaft diameter the gears are bored for [mm].
        count: Tooth counts to keep on each side, per module.
        material: Gear material.
        safety: Design factors.

    Returns:
        ``[(module, teeth, inter_axle, pair), ...]``, feasible pairs first and
        each group ordered by distance from the request.  The list is never
        empty: when nothing is feasible the infeasible candidates are still
        returned, least-overloaded first, so a caller can start from the best
        available and report that none was buildable.
    """
    candidates = []
    for module, teeth, value in lattice_neighbours(inter_axle, count=count):
        pair = size_pair(
            value,
            tangential_force,
            module=module,
            teeth=teeth,
            shaft_bore=shaft_bore,
            material=material,
            safety=safety,
        )
        candidates.append((module, teeth, value, pair))
    candidates.sort(
        key=lambda item: (
            not item[3].feasible,
            item[3].width_factor if not item[3].feasible else 0.0,
            abs(item[2] - float(inter_axle)),
        )
    )
    return candidates
