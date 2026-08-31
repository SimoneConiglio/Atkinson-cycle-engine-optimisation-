"""What the engine actually weighs, rather than what its seven links weigh.

The coupled sizing loop returns the mass of the sized members: about 0.25 kg.
No engine weighs 0.25 kg.  Reporting that number as "the mass" and then
optimising it is a category error -- it optimises a quarter-kilogram tail while
a four-kilogram dog goes unmodelled, and worse, it hides the fact that the two
*envelope* objectives ``H`` and ``B`` are themselves mass in disguise.

This module assembles the whole reciprocating assembly.  Eight contributions,
each sized from something the analysis already knows:

===================  =====================================================
contribution         what sets it
===================  =====================================================
linkage members      the coupled sizing fixed point
piston assembly      peak gas pressure
gear pair            peak tooth load, and the chosen module
shafts               combined bending and torsion at the journals
bearings             journal diameter, via a catalogue regression
crankcase            **the envelope H x B**, at minimum castable wall
cylinder and head    bore and peak pressure
flywheel             the cyclic torque fluctuation the linkage produces
===================  =====================================================

Two of these deserve emphasis, because they are what change the shape of the
optimization problem rather than just its scale.

**The crankcase makes H and B physical.**  A box has to enclose the mechanism,
and its walls scale with the envelope area.  In the geometric problem ``H`` and
``B`` are two of three competing objectives with no stated exchange rate
against efficiency -- which is precisely why that problem could only produce a
Pareto front and never a design.  Here they convert to kilograms at a rate the
physics fixes, and kilograms convert to range.  The three-objective problem
collapses to one objective, and the exchange rate is no longer the designer's
to choose.

**The flywheel makes torque smoothness physical.**  A single-cylinder engine
needs enough rotating inertia to carry it through the compression stroke.  The
requirement follows from the *fluctuation* of the torque the linkage produces,
so a design whose torque curve is flat needs less flywheel and therefore less
mass -- a design driver that no constraint in the geometric problem expresses,
and one that pushes in the opposite direction from a long lever arm.  It is
also frequently the single heaviest item in the budget.

Fidelity
--------
Every item is a first-order sizing model of the kind used in conceptual design:
right scaling laws, right sensitivities, absolute values good to perhaps 20 %.
The budget is not a weight statement for a built engine.  It is a consistent
way to convert design changes into kilograms, which is what an objective needs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import DEFAULT_SPEC, EngineSpec
from .derivatives import spectral_derivative
from .dynamics import DynamicLoads
from .gears import GearPair, size_pair
from .manufacturing import (
    MIN_MACHINED_THICKNESS,
    MIN_WALL_THICKNESS,
)
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors

ALUMINIUM_DENSITY = 2.70e-9
"""Density of a cast aluminium crankcase alloy [tonne/mm^3]."""

ALUMINIUM_YIELD = 200.0
"""Yield strength of the same alloy [MPa]."""

SPEED_FLUCTUATION = 0.10
"""Coefficient of cyclic speed fluctuation ``delta`` the flywheel must hold [-].

``delta = (omega_max - omega_min) / omega_mean``.  Ten per cent suits a road
vehicle, where the wheels and the driveline are themselves part of the
smoothing and nothing is driving a generator; a machine tool would want half
that.  The flywheel is frequently the heaviest item in the budget and its mass
goes as ``1 / delta``, so this number is swept rather than trusted.
"""

CASE_SLENDERNESS = 60.0
"""Span-to-wall ratio a cast housing is stiffened to [-].

A cast case is not sized by stress -- at these loads a 2.5 mm wall is
overwhelmingly strong -- but by stiffness, castability and the need to survive
handling.  Housings are conventionally proportioned at a span-to-thickness
ratio of 50 to 80.
"""

BEARING_FILL = 0.55
"""Fraction of a deep-groove bearing's swept volume that is steel [-]."""

CASE_CLEARANCE = 8.0
"""Clearance between the mechanism envelope and the inside of the case [mm]."""

FLYWHEEL_WEB_FACTOR = 1.25
"""Multiplier turning a bare rim mass into a rim with web and hub [-]."""


@dataclass(frozen=True)
class MassBudget:
    """Mass of the whole engine, itemised."""

    items: dict[str, float] = field(default_factory=dict)
    """Each contribution [tonne]."""

    gears: GearPair | None = None
    shaft_diameter: float = 0.0
    """Journal diameter both shafts are sized to [mm]."""

    flywheel_inertia: float = 0.0
    """Flywheel polar inertia actually fitted [tonne mm^2]."""

    flywheel_radius: float = 0.0
    """Radius the flywheel rim sits at [mm]."""

    required_inertia: float = 0.0
    """Rotating inertia the speed-fluctuation limit demands [tonne mm^2]."""

    inherent_inertia: float = 0.0
    """Rotating inertia already present in the cranks and gears [tonne mm^2]."""

    @property
    def total(self) -> float:
        """Total engine mass [tonne]."""
        return float(sum(self.items.values()))

    @property
    def total_kg(self) -> float:
        """Total engine mass [kg]."""
        return 1000.0 * self.total

    def kilograms(self) -> dict[str, float]:
        """The itemised budget in kilograms, heaviest first."""
        return dict(
            sorted(
                ((name, 1000.0 * value) for name, value in self.items.items()),
                key=lambda item: -item[1],
            )
        )

    def shares(self) -> dict[str, float]:
        """Each item as a fraction of the total."""
        total = self.total
        if total <= 0.0:
            return {}
        return {name: value / total for name, value in self.items.items()}


def gas_torque(loads: DynamicLoads) -> np.ndarray:
    """The crank torque the gas force alone would produce [N.mm].

    By virtual work, ``M = -F_gas d(lambda)/d(theta_1)``.  This is the classic
    turning-moment diagram, and it is what the flywheel calculation needs: the
    *inertia* part of the full torque history represents energy traded with the
    mechanism's own moving masses, which is accounted for separately as
    inherent rotating inertia.  Feeding the full torque into a flywheel sizing
    double-counts that exchange and, at low speed, overstates the flywheel by
    an order of magnitude.

    Args:
        loads: A solved dynamic load case.

    Returns:
        Gas torque at each crank angle [N.mm].
    """
    slope = spectral_derivative(loads.kinematics.lam, 1)
    return -np.asarray(loads.gas_force, dtype=float) * slope


def flywheel_requirement(
    torque: np.ndarray,
    speed: float,
    fluctuation: float = SPEED_FLUCTUATION,
    span: float = 2.0 * math.pi,
) -> tuple[float, float]:
    """Rotating inertia needed to hold the cyclic speed fluctuation.

    The standard construction: integrate the excess of instantaneous torque
    over its mean around the cycle, take the peak-to-peak of that accumulated
    energy, and divide by ``delta omega^2``.

    .. math:: J = \\frac{\\Delta E}{\\delta \\, \\omega^2},
        \\qquad \\Delta E = \\max_\\theta E(\\theta) - \\min_\\theta E(\\theta),
        \\qquad E(\\theta) = \\int_0^\\theta (M_r - \\bar{M_r}) \\, d\\theta_1

    Args:
        torque: Output torque at uniformly spaced crank angles [N.mm].
        speed: Mean crankshaft speed ``omega`` [rad/s].
        fluctuation: Allowed ``delta``.
        span: Crank angle the samples cover [rad].  ``2 pi`` for a cycle that
            completes in one crankshaft revolution, ``4 pi`` for a
            conventional four-stroke, whose flywheel has to carry it through
            two revolutions on one firing and is correspondingly larger.

    Returns:
        ``(required_inertia, energy_swing)`` in tonne mm^2 and N.mm.  The
        inertia is zero at zero speed, where the concept does not apply.
    """
    values = np.asarray(torque, dtype=float)
    n = values.size
    step = span / n
    excess = values - float(np.mean(values))
    accumulated = np.cumsum(excess) * step
    swing = float(np.max(accumulated) - np.min(accumulated))
    if speed <= 0.0 or fluctuation <= 0.0:
        return 0.0, swing
    return swing / (fluctuation * speed**2), swing


def shaft_diameter(
    peak_reaction: float,
    peak_torque: float,
    overhang: float,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
) -> float:
    """Journal diameter from combined bending and torsion.

    The ASME shaft equation with fatigue factors ``k_b = 1.5``, ``k_t = 1.0``:

    .. math:: d^3 = \\frac{16}{\\pi \\tau_a}
        \\sqrt{(k_b M)^2 + (k_t T)^2}

    Args:
        peak_reaction: Largest journal reaction over the cycle [N].
        peak_torque: Largest shaft torque [N.mm].
        overhang: Distance from the journal to the load [mm].
        material: Shaft material.
        safety: Design factors.

    Returns:
        Required journal diameter [mm], floored at 6 mm.
    """
    allowable_shear = 0.30 * material.yield_strength / safety.static
    bending = 1.5 * abs(peak_reaction) * max(overhang, 1.0)
    torsion = abs(peak_torque)
    cubed = 16.0 * math.hypot(bending, torsion) / (math.pi * allowable_shear)
    return max(float(np.cbrt(cubed)), 6.0)


def bearing_mass(bore: float, count: int, density: float) -> float:
    """Mass of a set of deep-groove ball bearings [tonne].

    Outer diameter and width follow the usual catalogue proportions for a
    medium series, ``D ~ 2.2 d + 8`` and ``B ~ 0.55 d + 3``; the swept annulus
    is then filled to :data:`BEARING_FILL` with steel to account for the
    raceways, balls and cage.

    Args:
        bore: Bearing bore, i.e. journal diameter [mm].
        count: How many bearings.
        density: Steel density [tonne/mm^3].

    Returns:
        Total mass [tonne].
    """
    outer = 2.2 * bore + 8.0
    width = 0.55 * bore + 3.0
    swept = math.pi * (outer**2 - bore**2) / 4.0 * width
    return density * BEARING_FILL * swept * count


def crankcase_mass(
    height: float,
    width: float,
    depth: float,
    peak_reaction: float,
    thickness: float | None = None,
) -> tuple[float, float]:
    """Mass of a thin-walled cast crankcase enclosing the mechanism.

    This is where the envelope objectives become kilograms.  The case is a
    closed box one clearance larger than the mechanism on every side.  Its wall
    is the greater of the minimum castable thickness and what a plate bending
    check needs to carry the main-bearing reaction into the structure, treating
    the bearing boss as a load on a plate of span ``min(H, B)``:

    .. math:: t = \\sqrt{\\frac{3 F \\, s}{2 \\, \\sigma_a \\, w}}

    with ``s`` the span, ``w`` an effective load-spreading width, and
    ``sigma_a`` the allowable.

    Args:
        height: Mechanism envelope ``H`` [mm].
        width: Mechanism envelope ``B`` [mm].
        depth: Case depth, across the crankshaft axis [mm].
        peak_reaction: Largest main-bearing reaction [N].
        thickness: Wall thickness [mm]; derived if omitted.

    Returns:
        ``(mass_tonne, wall_thickness_mm)``.
    """
    outer_h = height + 2.0 * CASE_CLEARANCE
    outer_b = width + 2.0 * CASE_CLEARANCE
    span = max(min(outer_h, outer_b), 1.0)
    stiffness = span / CASE_SLENDERNESS
    wall = max(stiffness, MIN_WALL_THICKNESS) if thickness is None else float(thickness)

    area = 2.0 * (outer_h * outer_b + outer_h * depth + outer_b * depth)
    walls = ALUMINIUM_DENSITY * area * wall

    # The main-bearing reactions do not go through the wall in bending: they go
    # into a local boss that spreads the load over a few journal diameters.
    # The boss is sized as a ring in bearing against the case, which is what
    # sets its thickness at these loads.
    allowable = ALUMINIUM_YIELD / 2.0
    boss_thickness = max(abs(peak_reaction) / (allowable * max(span, 1.0)), MIN_WALL_THICKNESS)
    boss_radius = 1.6 * CASE_CLEARANCE + 0.5 * span / CASE_SLENDERNESS
    bosses = ALUMINIUM_DENSITY * 4.0 * math.pi * boss_radius**2 * boss_thickness

    return walls + bosses, wall


def cylinder_mass(
    stroke: float,
    peak_pressure: float,
    spec: EngineSpec = DEFAULT_SPEC,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
) -> float:
    """Mass of the liner and head [tonne].

    The liner is a thick-walled tube under internal pressure (Lame, thin-wall
    limit); the head is a clamped circular plate over the bore.

    Args:
        stroke: Total piston travel the liner must guide [mm].
        peak_pressure: Peak in-cylinder gauge pressure [MPa].
        spec: Fixed engine data.
        material: Liner and head material.
        safety: Design factors.

    Returns:
        Mass of liner plus head [tonne].
    """
    allowable = material.yield_strength / safety.static
    radius = 0.5 * spec.bore
    liner_wall = max(peak_pressure * radius / allowable, MIN_MACHINED_THICKNESS)
    liner_length = stroke + spec.piston_length + 2.0 * CASE_CLEARANCE
    liner_volume = math.pi * ((radius + liner_wall) ** 2 - radius**2) * liner_length

    head_thickness = max(
        radius * math.sqrt(3.0 * peak_pressure / (4.0 * allowable)),
        3.0 * MIN_MACHINED_THICKNESS,
    )
    head_volume = math.pi * (radius + liner_wall + 6.0) ** 2 * head_thickness
    return material.density * (liner_volume + head_volume)


def assemble(
    loads: DynamicLoads,
    diameters: dict[str, float],
    member_mass: dict[str, float],
    piston_mass: float,
    height: float,
    width: float,
    stroke: float,
    peak_pressure: float,
    module: float | None = None,
    teeth: int | None = None,
    fluctuation: float = SPEED_FLUCTUATION,
    spec: EngineSpec = DEFAULT_SPEC,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
) -> MassBudget:
    """Assemble the full engine mass budget.

    Args:
        loads: A solved dynamic load case.
        diameters: Sized member diameters [mm].
        member_mass: Mass of each sized member [tonne].
        piston_mass: Piston assembly mass [tonne].
        height: Mechanism envelope ``H`` [mm].
        width: Mechanism envelope ``B`` [mm].
        stroke: Total piston travel [mm].
        peak_pressure: Peak in-cylinder gauge pressure [MPa].
        module: Gear module to use [mm]; chosen automatically if omitted.
        teeth: Teeth on the small gear; derived from ``I`` if omitted.
        fluctuation: Allowed cyclic speed fluctuation.
        spec: Fixed engine data.
        material: Steel parts' material.
        safety: Design factors.

    Returns:
        The itemised budget.
    """
    design = loads.kinematics.design

    peak_reaction = float(np.max(np.linalg.norm(loads.reaction["R1"], axis=1)))
    peak_reaction_2 = float(np.max(np.linalg.norm(loads.reaction["R2"], axis=1)))
    peak_torque = float(np.max(np.abs(loads.torque)))
    peak_tangential = float(np.max(np.abs(loads.gear_force))) * math.cos(spec.pressure_angle)

    # -- shafts and bearings ------------------------------------------------------
    journal = shaft_diameter(
        max(peak_reaction, peak_reaction_2),
        peak_torque,
        max(design.q_1, design.q_2),
        material,
        safety,
    )
    case_depth = spec.bore + 2.0 * (0.55 * journal + 3.0) + 4.0 * CASE_CLEARANCE
    shaft_length = case_depth + 2.0 * (0.55 * journal + 3.0)
    shaft = material.density * math.pi * journal**2 / 4.0 * shaft_length * 2.0
    bearings = bearing_mass(journal, 4, material.density)

    # -- gears --------------------------------------------------------------------
    pair = size_pair(
        design.I,
        peak_tangential,
        module=module,
        teeth=teeth,
        shaft_bore=journal,
        material=material,
        safety=safety,
    )

    # -- flywheel -----------------------------------------------------------------
    required, _swing = flywheel_requirement(gas_torque(loads), loads.speed, fluctuation)
    # What already turns with the crankshaft: the large gear, and the crank
    # throw treated as a bar rotating about the shaft axis.
    crank_inertia = member_mass.get("crank_1", 0.0) * design.q_1**2 / 3.0
    inherent = pair.inertia_large + crank_inertia
    deficit = max(required - inherent, 0.0)
    # The flywheel lives on the crankshaft nose, so the case width is what
    # bounds it; a bigger wheel is always lighter for the same inertia, so the
    # bound is active whenever a flywheel is needed at all.
    radius = min(max(0.45 * width, 30.0), 150.0)
    flywheel = FLYWHEEL_WEB_FACTOR * deficit / radius**2 if radius > 0.0 else 0.0

    # -- case, cylinder -----------------------------------------------------------
    case, _wall = crankcase_mass(height, width, case_depth, peak_reaction)
    cylinder = cylinder_mass(stroke, peak_pressure, spec, material, safety)

    items = {
        "linkage": float(sum(member_mass.values())),
        "piston": float(piston_mass),
        "gears": pair.mass,
        "shafts": float(shaft),
        "bearings": float(bearings),
        "crankcase": float(case),
        "cylinder_head": float(cylinder),
        "flywheel": float(flywheel),
    }
    return MassBudget(
        items=items,
        gears=pair,
        shaft_diameter=journal,
        flywheel_inertia=deficit,
        flywheel_radius=radius,
        required_inertia=required,
        inherent_inertia=inherent,
    )
