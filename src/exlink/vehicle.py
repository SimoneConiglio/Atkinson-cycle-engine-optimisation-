"""Range: the one objective the whole problem was always about.

The geometric problem has three objectives -- maximise ``eta``, minimise ``H``,
minimise ``B`` -- and no exchange rate between them.  That is why it can only
ever produce a Pareto front, and why picking a point on that front is left to
the designer's taste.  But the engine is for a Shell Eco-marathon car, and that
competition has exactly one figure of merit: **how far you get on a given
quantity of fuel**.  Once that is written down, the exchange rates are no
longer anyone's taste.  They are physics.

.. code-block:: text

    eta_mech  --.
                +--> brake efficiency --.
    cycle work -'                        +--> fuel per metre --> RANGE
                                        /
    H, B  --> crankcase --.            /
    torque ripple --> flywheel --> engine mass --> rolling resistance

Efficiency and mass both land on range, in commensurable units, and they pull
in opposite directions: a long lever arm flatters ``eta`` but builds a big
crankcase and a heavy flywheel, and the heavier car then drags that back out
through rolling resistance.

How the car is driven
---------------------
Not at constant speed.  Every serious Eco-marathon team drives **burn and
coast**: run the engine hard from ``v_lo`` up to ``v_hi``, declutch, coast back
down to ``v_lo``, repeat.  This is not a detail -- it is the reason the engine
can be optimised at one operating point at all, and it is where the
accelerations and decelerations enter.

The energy bookkeeping is worth stating, because the naive expectation is
wrong.  Over one full burn-and-coast cycle the car starts and ends at the same
speed, so the kinetic energy nets to zero and

.. math:: W_{\\text{burn}} = \\int_{\\text{burn} + \\text{coast}} F_{\\text{res}}(v) \\, dx

Accelerating hard costs **nothing** in resistance work.  What burn-and-coast
actually buys is that the engine spends its whole running time at high load,
where brake efficiency is far better than at the part load a constant-speed
drive would need.  What it costs is aerodynamic: drag goes as ``v^2``, so
swinging between ``v_lo`` and ``v_hi`` burns more than cruising at their mean.
Both effects are in the model below, and the optimum window is the balance
between them.

What is *not* modelled
----------------------
Transient combustion during the burn, gear-change losses, the driver's line
around the track, wind, and the engine's own warm-up.  Rolling resistance is a
constant coefficient and aerodynamic drag a constant ``C_d A``; neither
responds to speed or load beyond the explicit ``v^2``.  The absolute range is
therefore good to perhaps 20 %, and the model earns its place by ranking
designs consistently, not by predicting a competition result.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .materials import FloatArray

GRAVITY = 9.80665
"""Standard gravity [m/s^2]."""

AIR_DENSITY = 1.20
"""Air density at competition conditions [kg/m^3]."""

GASOLINE_LHV = 44.0e6
"""Lower heating value of gasoline [J/kg]."""

GASOLINE_DENSITY = 745.0
"""Density of gasoline [kg/m^3]."""

ENERGY_PER_LITRE = GASOLINE_LHV * GASOLINE_DENSITY / 1000.0
"""Energy in one litre of gasoline [J/L]."""


@dataclass(frozen=True)
class Vehicle:
    """The car the engine goes into, everything except the engine.

    Defaults describe a Shell Eco-marathon *Prototype* class entry: a
    three-wheeled, fully faired, single-seat lay-down car.
    """

    glider_mass: float = 35.0
    """Chassis, body, wheels and driveline, without the engine [kg]."""

    driver_mass: float = 50.0
    """Driver plus equipment [kg]; the competition sets a 50 kg minimum."""

    rolling_resistance: float = 0.0015
    """``C_rr`` for competition tyres on smooth asphalt [-]."""

    drag_area: float = 0.075
    """``C_d A`` [m^2]; roughly ``C_d = 0.15`` over a 0.5 m^2 frontal area."""

    transmission_efficiency: float = 0.90
    """Crankshaft to driven wheel, including the clutch [-]."""

    minimum_average_speed: float = 25.0
    """Average speed the competition requires [km/h]."""

    def mass(self, engine_mass_kg: float) -> float:
        """Total moving mass with a given engine [kg]."""
        return self.glider_mass + self.driver_mass + float(engine_mass_kg)

    def resistance(self, speed: FloatArray | float, total_mass: float) -> FloatArray:
        """Road-load force at a speed [N].

        Args:
            speed: Road speed [m/s].
            total_mass: Vehicle mass including engine and driver [kg].

        Returns:
            ``C_rr m g + 0.5 rho C_d A v^2`` [N].
        """
        v = np.asarray(speed, dtype=float)
        rolling = self.rolling_resistance * total_mass * GRAVITY
        aero = 0.5 * AIR_DENSITY * self.drag_area * v**2
        return rolling + aero


@dataclass(frozen=True)
class RangeResult:
    """Outcome of a burn-and-coast range calculation."""

    distance_per_joule: float
    """Distance covered per joule of fuel energy [m/J]."""

    km_per_litre: float
    """The competition's figure of merit [km/L]."""

    average_speed: float
    """Mean road speed over the burn-and-coast cycle [km/h]."""

    burn_fraction: float
    """Share of the cycle *distance* covered under power [-]."""

    burn_distance: float
    """Distance covered under power, one cycle [m]."""

    coast_distance: float
    """Distance covered coasting, one cycle [m]."""

    total_mass: float
    """Vehicle mass including engine and driver [kg]."""

    engine_mass: float
    """Engine mass [kg]."""

    brake_power: float
    """Brake power at the chosen operating point [W]."""

    brake_efficiency: float
    """Brake work over fuel energy [-]."""

    feasible: bool
    """Whether the engine can drive the car and meet the speed rule."""

    reason: str = ""
    """Why not, when ``feasible`` is false."""


def brake_efficiency(brake_work: float, heat_release: float) -> float:
    """Brake thermal efficiency from one cycle's works.

    Args:
        brake_work: Work leaving the crankshaft per cycle [N.mm].
        heat_release: Heat added per cycle [N.mm].

    Returns:
        ``W_brake / Q``, zero when either is non-positive.
    """
    if heat_release <= 0.0 or brake_work <= 0.0:
        return 0.0
    return brake_work / heat_release


def heat_release(
    dead_volume: float,
    p_compression_end: float,
    p_combustion: float,
    heat_capacity_ratio: float,
) -> float:
    """Heat added by the constant-volume combustion step [N.mm].

    .. math:: Q = \\frac{V_0 (P_3 - P_2)}{\\gamma - 1}

    This, not an assumed air-fuel ratio, is what fixes the fuel consumed:
    it is the only heat input the idealised cycle has, so deriving the fuel
    from it keeps the range model consistent with the thermodynamics rather
    than bolting a second, independent combustion model on the side.

    Args:
        dead_volume: Clearance volume ``V_0`` [mm^3].
        p_compression_end: ``P_2`` [MPa].
        p_combustion: ``P_3`` [MPa].
        heat_capacity_ratio: ``gamma``.

    Returns:
        ``Q`` per cycle [N.mm].
    """
    return dead_volume * (p_combustion - p_compression_end) / (heat_capacity_ratio - 1.0)


def _quadrature(
    low: float, high: float, total_mass: float, net_force: FloatArray, speeds: FloatArray
) -> tuple[float, float]:
    """Integrate distance and time for a speed sweep under a net force.

    ``M v dv/dx = F`` and ``M dv/dt = F`` give
    ``dx = M v dv / F`` and ``dt = M dv / F``.
    """
    if high <= low:
        return 0.0, 0.0
    safe = np.where(
        np.abs(net_force) < 1.0e-9, np.sign(net_force) * 1.0e-9 + 1.0e-12, net_force
    )
    distance = float(np.trapezoid(total_mass * speeds / safe, speeds))
    time = float(np.trapezoid(total_mass / safe, speeds))
    return distance, time


def burn_and_coast(
    vehicle: Vehicle,
    engine_mass_kg: float,
    brake_power: float,
    efficiency: float,
    speed_low: float,
    speed_high: float,
    samples: int = 201,
) -> RangeResult:
    """Range under a burn-and-coast strategy.

    Args:
        vehicle: The car.
        engine_mass_kg: Engine mass [kg].
        brake_power: Brake power at the operating point [W].
        efficiency: Brake thermal efficiency [-].
        speed_low: Speed the coast ends at [m/s].
        speed_high: Speed the burn ends at [m/s].
        samples: Quadrature points per phase.

    Returns:
        The range result.  Check :attr:`RangeResult.feasible`.
    """
    total = vehicle.mass(engine_mass_kg)
    blank = RangeResult(
        distance_per_joule=0.0,
        km_per_litre=0.0,
        average_speed=0.0,
        burn_fraction=0.0,
        burn_distance=0.0,
        coast_distance=0.0,
        total_mass=total,
        engine_mass=float(engine_mass_kg),
        brake_power=float(brake_power),
        brake_efficiency=float(efficiency),
        feasible=False,
    )
    if brake_power <= 0.0 or efficiency <= 0.0:
        return RangeResult(**{**blank.__dict__, "reason": "engine produces no net work"})
    if not 0.0 < speed_low < speed_high:
        return RangeResult(**{**blank.__dict__, "reason": "invalid speed window"})

    speeds = np.linspace(speed_low, speed_high, samples)
    resistance = vehicle.resistance(speeds, total)
    wheel_power = brake_power * vehicle.transmission_efficiency
    traction = wheel_power / speeds

    net_burn = traction - resistance
    if np.min(net_burn) <= 0.0:
        return RangeResult(
            **{**blank.__dict__, "reason": "engine cannot accelerate the car to v_high"}
        )

    burn_distance, burn_time = _quadrature(speed_low, speed_high, total, net_burn, speeds)
    coast_distance, coast_time = _quadrature(speed_low, speed_high, total, resistance, speeds)

    distance = burn_distance + coast_distance
    time = burn_time + coast_time
    if distance <= 0.0 or time <= 0.0:
        return RangeResult(**{**blank.__dict__, "reason": "degenerate burn-coast cycle"})

    fuel_energy = brake_power / efficiency * burn_time
    per_joule = distance / fuel_energy
    average = 3.6 * distance / time

    return RangeResult(
        distance_per_joule=per_joule,
        km_per_litre=per_joule * ENERGY_PER_LITRE / 1000.0,
        average_speed=average,
        burn_fraction=burn_distance / distance,
        burn_distance=burn_distance,
        coast_distance=coast_distance,
        total_mass=total,
        engine_mass=float(engine_mass_kg),
        brake_power=float(brake_power),
        brake_efficiency=float(efficiency),
        feasible=average >= vehicle.minimum_average_speed,
        reason=""
        if average >= vehicle.minimum_average_speed
        else "below minimum average speed",
    )


def best_strategy(
    vehicle: Vehicle,
    engine_mass_kg: float,
    brake_power: float,
    efficiency: float,
    window: tuple[float, float] = (4.0, 16.0),
    steps: int = 24,
) -> RangeResult:
    """Search the burn-and-coast window for the best range.

    A coarse grid over ``(v_low, v_high)``, refined once around the winner.
    The surface is smooth and single-peaked in practice -- aerodynamic drag
    pushes the window down, the minimum-average-speed rule pushes it up -- so a
    grid is both adequate and robust, and it never fails the way a gradient
    method would at the feasibility boundary.

    Args:
        vehicle: The car.
        engine_mass_kg: Engine mass [kg].
        brake_power: Brake power [W].
        efficiency: Brake thermal efficiency [-].
        window: Speeds to search between [m/s].
        steps: Grid points per axis.

    Returns:
        The best feasible result, or the least-bad infeasible one.
    """
    low_bound, high_bound = window
    best: RangeResult | None = None
    fallback: RangeResult | None = None

    def consider(v_low: float, v_high: float) -> None:
        nonlocal best, fallback
        result = burn_and_coast(vehicle, engine_mass_kg, brake_power, efficiency, v_low, v_high)
        if result.feasible:
            if best is None or result.km_per_litre > best.km_per_litre:
                best = result
        elif fallback is None or result.km_per_litre > fallback.km_per_litre:
            fallback = result

    grid = np.linspace(low_bound, high_bound, steps)
    for i, v_low in enumerate(grid[:-1]):
        for v_high in grid[i + 1 :]:
            consider(float(v_low), float(v_high))

    if best is not None:
        span = (high_bound - low_bound) / (steps - 1)
        anchor_low, anchor_high = None, None
        # Recover the winning window by re-deriving it from the stored result:
        # the average speed and burn fraction identify it uniquely enough for a
        # local refinement, so simply re-scan a neighbourhood of the best pair.
        for i, v_low in enumerate(grid[:-1]):
            for v_high in grid[i + 1 :]:
                candidate = burn_and_coast(
                    vehicle,
                    engine_mass_kg,
                    brake_power,
                    efficiency,
                    float(v_low),
                    float(v_high),
                )
                if candidate.feasible and math.isclose(
                    candidate.km_per_litre, best.km_per_litre, rel_tol=1e-12
                ):
                    anchor_low, anchor_high = float(v_low), float(v_high)
        if anchor_low is not None and anchor_high is not None:
            fine_low = np.linspace(max(anchor_low - span, 0.5), anchor_low + span, 7)
            fine_high = np.linspace(anchor_high - span, anchor_high + span, 7)
            for v_low in fine_low:
                for v_high in fine_high:
                    if v_high > v_low:
                        consider(float(v_low), float(v_high))
        return best
    return (
        fallback
        if fallback is not None
        else burn_and_coast(
            vehicle, engine_mass_kg, brake_power, efficiency, low_bound, high_bound
        )
    )
