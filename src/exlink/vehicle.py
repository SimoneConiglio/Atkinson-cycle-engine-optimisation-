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


def _average_speed(
    vehicle: Vehicle,
    engine_mass_kg: float,
    brake_power: float,
    efficiency: float,
    speed_low: float,
    speed_high: float,
) -> float:
    """Average speed of a burn-and-coast cycle, or ``-inf`` if it cannot run."""
    result = burn_and_coast(
        vehicle, engine_mass_kg, brake_power, efficiency, speed_low, speed_high
    )
    if result.average_speed <= 0.0:
        return -math.inf
    return result.average_speed


def top_speed(
    vehicle: Vehicle,
    engine_mass_kg: float,
    brake_power: float,
    ceiling: float = 40.0,
) -> float:
    """The speed at which traction and road load balance [m/s].

    A burn phase cannot be asked to reach beyond this, and the burn-and-coast
    window is therefore bounded by it rather than by an arbitrary constant.

    Args:
        vehicle: The car.
        engine_mass_kg: Engine mass [kg].
        brake_power: Brake power [W].
        ceiling: Upper bracket for the search [m/s].

    Returns:
        The balance speed, or 0 when the engine cannot move the car at all.
    """
    total = vehicle.mass(engine_mass_kg)
    wheel_power = brake_power * vehicle.transmission_efficiency
    if wheel_power <= 0.0:
        return 0.0

    def surplus(speed: float) -> float:
        return wheel_power / speed - float(vehicle.resistance(speed, total))

    low, high = 1.0e-3, ceiling
    if surplus(low) <= 0.0:
        return 0.0
    if surplus(high) > 0.0:
        return high
    for _ in range(80):
        middle = 0.5 * (low + high)
        if surplus(middle) > 0.0:
            low = middle
        else:
            high = middle
    return low


def _high_speed_for_rule(
    vehicle: Vehicle,
    engine_mass_kg: float,
    brake_power: float,
    efficiency: float,
    speed_low: float,
    ceiling: float,
    tolerance: float = 1.0e-6,
) -> float | None:
    """The ``v_high`` at which the average-speed rule is exactly met.

    Average speed rises monotonically with ``v_high`` at fixed ``v_low`` -- a
    wider window spends more of the cycle fast -- so a bisection is both valid
    and robust.  Solving the rule *exactly* rather than searching a grid is
    what keeps the objective smooth enough to differentiate.

    Returns:
        The bracketing speed, or ``None`` if the rule cannot be met at this
        ``v_low`` with this engine.
    """
    target = vehicle.minimum_average_speed
    low, high = speed_low * (1.0 + 1.0e-6), ceiling
    if (
        _average_speed(vehicle, engine_mass_kg, brake_power, efficiency, speed_low, high)
        < target
    ):
        return None
    if (
        _average_speed(vehicle, engine_mass_kg, brake_power, efficiency, speed_low, low)
        >= target
    ):
        return low
    for _ in range(60):
        middle = 0.5 * (low + high)
        if high - low < tolerance:
            break
        value = _average_speed(
            vehicle, engine_mass_kg, brake_power, efficiency, speed_low, middle
        )
        if value >= target:
            high = middle
        else:
            low = middle
    return high


def best_strategy(
    vehicle: Vehicle,
    engine_mass_kg: float,
    brake_power: float,
    efficiency: float,
    window: tuple[float, float] = (2.0, 25.0),
    steps: int = 24,
) -> RangeResult:
    """Choose the burn-and-coast window that maximises range.

    The structure of this sub-problem is worth stating, because it makes the
    search both fast and smooth.

    Range always improves as the whole speed window comes down: rolling
    resistance is speed-independent and aerodynamic drag is not, so the
    cheapest way to cover a metre is to cover it slowly.  What stops the
    optimizer driving the window to zero is the competition's
    **minimum average speed**, and that constraint is therefore active at
    every optimum where the engine has power to spare.

    So the two-dimensional search collapses to one dimension.  For each
    ``v_low``, the rule pins ``v_high`` exactly -- average speed is monotone in
    ``v_high``, so a bisection finds it -- and what remains is a scalar
    minimisation over ``v_low``, trading a narrow window at a higher mean speed
    against a wide one swinging around a lower mean.

    Solving the active constraint rather than sampling a grid matters for more
    than speed.  A grid makes the objective a step function of the design
    variables, and a step function has no useful finite difference; the
    gradient-based optimizer downstream would be differentiating quantisation
    noise.  Here the objective is smooth.

    Args:
        vehicle: The car.
        engine_mass_kg: Engine mass [kg].
        brake_power: Brake power [W].
        efficiency: Brake thermal efficiency [-].
        window: Speeds to search between [m/s].
        steps: Scalar-search resolution over ``v_low``.

    Returns:
        The best feasible result, or the least-bad infeasible one -- never an
        exception, so that an optimizer can walk through infeasible ground.
    """
    low_bound, high_bound = window
    # A burn cannot outrun the balance speed, so that -- not the nominal window
    # -- is what bounds v_high.  A margin keeps the acceleration away from zero,
    # where the burn distance integral diverges.
    reachable = min(high_bound, 0.98 * top_speed(vehicle, engine_mass_kg, brake_power))
    infeasible = burn_and_coast(
        vehicle,
        engine_mass_kg,
        brake_power,
        efficiency,
        low_bound,
        max(reachable, low_bound + 0.1),
    )
    if reachable <= low_bound:
        return infeasible

    def evaluate_low(v_low: float) -> RangeResult | None:
        v_high = _high_speed_for_rule(
            vehicle, engine_mass_kg, brake_power, efficiency, v_low, reachable
        )
        if v_high is None:
            return None
        result = burn_and_coast(vehicle, engine_mass_kg, brake_power, efficiency, v_low, v_high)
        return result if result.feasible else None

    ceiling = min(reachable, vehicle.minimum_average_speed / 3.6)
    candidates = np.linspace(low_bound, ceiling * 0.995, steps)
    best: RangeResult | None = None
    best_low = 0.0
    for v_low in candidates:
        result = evaluate_low(float(v_low))
        if result is not None and (best is None or result.km_per_litre > best.km_per_litre):
            best, best_low = result, float(v_low)

    if best is None:
        return infeasible

    # Golden-section refinement around the winner, on the same smooth objective.
    span = (ceiling * 0.995 - low_bound) / (steps - 1)
    left, right = max(best_low - span, low_bound), min(best_low + span, ceiling * 0.995)
    phi = 0.5 * (math.sqrt(5.0) - 1.0)
    for _ in range(40):
        if right - left < 1.0e-6:
            break
        a = right - phi * (right - left)
        b = left + phi * (right - left)
        value_a = evaluate_low(a)
        value_b = evaluate_low(b)
        score_a = value_a.km_per_litre if value_a else -math.inf
        score_b = value_b.km_per_litre if value_b else -math.inf
        for value in (value_a, value_b):
            if value is not None and value.km_per_litre > best.km_per_litre:
                best = value
        if score_a > score_b:
            right = b
        else:
            left = a
    return best
