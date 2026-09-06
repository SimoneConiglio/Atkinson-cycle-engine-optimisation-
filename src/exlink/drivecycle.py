"""Scoring a design over a schedule of operating points, not one.

Everything in sections 6.1 to 6.4 is reported at a single crankshaft speed with
a sweep around it, and section 7.2 lists that as a limitation: an optimum found
at one point may be an artefact of the point.  This module scores a design over
a weighted schedule instead, and the aggregation is where the interest is.

Why a weighted harmonic mean
----------------------------
The figure of merit is distance per unit fuel, so fuel per unit distance is what
adds.  Over a schedule whose points carry distance shares ``w_i``,

.. math:: \\frac{1}{R} = \\sum_i \\frac{w_i}{R_i}, \\qquad \\sum_i w_i = 1

which is the harmonic mean weighted by distance, not the arithmetic mean of the
``R_i``.  The difference is not cosmetic: the arithmetic mean flatters a design
whose range collapses at one point, and a schedule exists precisely to catch
that.

Why one engine, not one per point
---------------------------------
Scoring each point with :func:`exlink.performance.evaluate` and averaging is
wrong, and wrong in the optimistic direction, because it lets the *engine*
change from point to point: the sections shrink at the point where the inertia
is small and the flywheel shrinks at the point where the speed is high.  A real
engine is one engine.

The two ends of the schedule size different parts of it, and they pull opposite
ways:

* **the fastest point sizes the structure.**  Every inertia load grows as the
  square of speed and the sections follow, so the members, bearings and case
  must be those of the highest speed in the schedule.
* **the slowest point sizes the flywheel.**  The rotating inertia needed to
  hold a given speed fluctuation goes as the *inverse* square of speed, so the
  wheel is largest at the bottom of the schedule -- and on this engine it is
  frequently the heaviest single item.

:func:`score_cycle` therefore takes the structure from the fastest point and
the flywheel requirement from whichever point demands most, assembles one
engine from the two, and scores every point with that mass.  Sizing the
structure at the top of the schedule is exact in this model rather than
conservative: the flywheel is concentric with its shaft, so fitting a larger
one adds no inertia force and cannot feed back into the member loads
(:mod:`exlink.dynamics` states the same for shafts and gears).
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, replace
from typing import Any

from .constants import DEFAULT_SPEC, EngineSpec
from .design import Design
from .mass_budget import FLYWHEEL_WEB_FACTOR, MassBudget
from .performance import Performance, evaluate
from .slidercrank import SliderCrank
from .slidercrank import friction_work as slidercrank_friction_work
from .slidercrank import mass_budget as slidercrank_mass_budget
from .slidercrank import solve as solve_slidercrank
from .vehicle import Vehicle, best_strategy, brake_efficiency


@dataclass(frozen=True)
class CyclePoint:
    """One operating point of a schedule."""

    crankshaft_rpm: float
    """Speed of the shaft power is taken from [rev/min]."""

    share: float
    """Share of the run's *distance* covered at this point [-]."""

    @property
    def speed_rpm(self) -> float:
        """The analysis speed, half the crankshaft's [rev/min]."""
        return self.crankshaft_rpm / DEFAULT_SPEC.output_revolutions_per_cycle


@dataclass(frozen=True)
class DriveCycle:
    """A weighted schedule of operating points."""

    points: tuple[CyclePoint, ...]
    name: str = ""

    def normalised(self) -> DriveCycle:
        """The same schedule with the shares summing to one."""
        total = sum(point.share for point in self.points)
        if total <= 0.0:
            msg = "a drive cycle needs at least one point with a positive share"
            raise ValueError(msg)
        return replace(
            self,
            points=tuple(replace(point, share=point.share / total) for point in self.points),
        )

    @property
    def fastest(self) -> CyclePoint:
        """The point that sizes the structure."""
        return max(self.points, key=lambda point: point.crankshaft_rpm)


STANDARD_CYCLE = DriveCycle(
    points=(
        CyclePoint(1600.0, 0.20),
        CyclePoint(2000.0, 0.40),
        CyclePoint(2400.0, 0.25),
        CyclePoint(2800.0, 0.15),
    ),
    name="four-point schedule about the design speed",
)
"""The schedule the study reports against.

A competition run is not held at one operating point: gradient, wind and
traffic move the power the car needs, and the driver's burn-and-coast rhythm
moves with it.  Without a surveyed track there is no defensible way to derive
the weights, so they are stated as an assumption rather than computed -- a
spread from 0.8 to 1.4 times the design speed, weighted towards it.  What the
schedule is for is not a more accurate number but a *test*: a design whose
range holds across it was not an artefact of the point it was found at.
"""


@dataclass(frozen=True)
class CycleRow:
    """What one point of the schedule contributed."""

    point: CyclePoint
    km_per_litre: float
    brake_power: float
    brake_efficiency: float
    feasible: bool
    reason: str = ""


@dataclass(frozen=True)
class CycleResult:
    """A design scored over a schedule."""

    cycle: DriveCycle
    rows: tuple[CycleRow, ...]
    engine_mass_kg: float
    """Mass of the one engine every point was scored with [kg]."""

    flywheel_kg: float
    """The flywheel inside that mass [kg], sized by the slowest point."""

    sizing_rpm: float
    """Crankshaft speed the structure was sized at [rev/min]."""

    @property
    def km_per_litre(self) -> float:
        """Distance-weighted harmonic mean of the per-point ranges [km/L]."""
        fuel = 0.0
        for row in self.rows:
            if row.km_per_litre <= 0.0:
                return 0.0
            fuel += row.point.share / row.km_per_litre
        return 0.0 if fuel <= 0.0 else 1.0 / fuel

    @property
    def feasible(self) -> bool:
        """Whether the engine drives the car at every point of the schedule."""
        return all(row.feasible for row in self.rows)

    @property
    def spread(self) -> float:
        """Worst point's range over the best point's [-].

        One near unity says the design is not specific to the point it was found
        at; one well below says it is.
        """
        ranges = [row.km_per_litre for row in self.rows]
        best = max(ranges, default=0.0)
        return 0.0 if best <= 0.0 else min(ranges) / best


def score_cycle(
    design: Design,
    cycle: DriveCycle = STANDARD_CYCLE,
    vehicle: Vehicle | None = None,
    spec: EngineSpec = DEFAULT_SPEC,
    **kwargs: Any,
) -> CycleResult:
    """Score one design over a schedule, with one engine.

    Args:
        design: The mechanism.
        cycle: The schedule; shares are normalised on the way in.
        vehicle: The car; a default Prototype-class entry if omitted.
        spec: Fixed engine data.
        **kwargs: Forwarded to :func:`exlink.performance.evaluate`, so the gear
            pair can be pinned as every run in the study pins it.

    Returns:
        The per-point rows and the aggregate.
    """
    schedule = cycle.normalised()
    car = vehicle if vehicle is not None else Vehicle()

    scored: dict[float, Performance] = {
        point.crankshaft_rpm: evaluate(
            design, speed_rpm=point.speed_rpm, vehicle=car, spec=spec, **kwargs
        )
        for point in schedule.points
    }

    mass_kg, flywheel_kg = _one_engine(
        scored[schedule.fastest.crankshaft_rpm].budget,
        (item.budget for item in scored.values()),
    )

    rows = []
    for point in schedule.points:
        performance = scored[point.crankshaft_rpm]
        outcome = best_strategy(
            car,
            mass_kg,
            performance.brake_power,
            performance.brake_efficiency,
        )
        rows.append(
            CycleRow(
                point=point,
                km_per_litre=outcome.km_per_litre,
                brake_power=performance.brake_power,
                brake_efficiency=performance.brake_efficiency,
                feasible=outcome.feasible and performance.analysis.valid,
                reason=outcome.reason,
            )
        )

    return CycleResult(
        cycle=schedule,
        rows=tuple(rows),
        engine_mass_kg=mass_kg,
        flywheel_kg=flywheel_kg,
        sizing_rpm=schedule.fastest.crankshaft_rpm,
    )


def _one_engine(structure: MassBudget, scored: Iterable[MassBudget]) -> tuple[float, float]:
    """Mass of the single engine the whole schedule must share [kg].

    The structure of the fastest point, with its flywheel replaced by the
    largest any point of the schedule demands.  Returns the total and the
    flywheel inside it.

    Args:
        structure: Budget of the fastest point, which sizes everything but the
            flywheel.
        scored: Budgets of every point, which between them size the flywheel.

    Returns:
        Total engine mass and the flywheel inside it [kg].
    """
    fitted = 1000.0 * structure.items.get("flywheel", 0.0)
    if structure.flywheel_radius <= 0.0:
        return structure.total_kg, fitted

    required = max(budget.required_inertia for budget in scored)
    deficit = max(required - structure.inherent_inertia, 0.0)
    heaviest = 1000.0 * FLYWHEEL_WEB_FACTOR * deficit / structure.flywheel_radius**2
    return structure.total_kg - fitted + heaviest, heaviest


def score_slidercrank_cycle(
    mechanism: SliderCrank,
    cycle: DriveCycle = STANDARD_CYCLE,
    vehicle: Vehicle | None = None,
    samples: int = 360,
    **kwargs: Any,
) -> CycleResult:
    """Score a conventional engine over the same schedule, the same way.

    The baseline of section 6.3 is optimised at one speed exactly as the
    EX-link is, so it inherits the same objection and has to answer it the
    same way.  The mechanics are identical -- one engine, structure from the
    fastest point, flywheel from whichever point demands most -- with the one
    difference that the slider-crank's crankshaft *is* its analysis shaft, so
    the schedule's speeds are used unhalved.

    Args:
        mechanism: The slider-crank.
        cycle: The schedule.
        vehicle: The car; a default Prototype-class entry if omitted.
        samples: Crank angles per revolution.
        **kwargs: Forwarded to :func:`exlink.slidercrank.solve`.

    Returns:
        The per-point rows and the aggregate.
    """
    schedule = cycle.normalised()
    car = vehicle if vehicle is not None else Vehicle()

    solved = {
        point.crankshaft_rpm: solve_slidercrank(
            mechanism, point.crankshaft_rpm, samples=samples, **kwargs
        )
        for point in schedule.points
    }
    budgets = {rpm: slidercrank_mass_budget(item) for rpm, item in solved.items()}
    mass_kg, flywheel_kg = _one_engine(
        budgets[schedule.fastest.crankshaft_rpm], budgets.values()
    )

    rows = []
    for point in schedule.points:
        result = solved[point.crankshaft_rpm]
        brake = result.indicated_work - slidercrank_friction_work(result)
        efficiency = brake_efficiency(brake, result.heat_release)
        # Two revolutions per cycle, so half a cycle's work per revolution.
        power = brake / 1000.0 * (point.crankshaft_rpm / 60.0) / 2.0
        outcome = best_strategy(car, mass_kg, power, efficiency)
        rows.append(
            CycleRow(
                point=point,
                km_per_litre=outcome.km_per_litre,
                brake_power=power,
                brake_efficiency=efficiency,
                feasible=result.converged and outcome.feasible and brake > 0.0,
                reason=outcome.reason,
            )
        )

    return CycleResult(
        cycle=schedule,
        rows=tuple(rows),
        engine_mass_kg=mass_kg,
        flywheel_kg=flywheel_kg,
        sizing_rpm=schedule.fastest.crankshaft_rpm,
    )


def format_cycle(result: CycleResult) -> str:
    """Render a :class:`CycleResult` as an aligned table."""
    lines = [f"drive cycle: {result.cycle.name or 'unnamed'}", "=" * 46, ""]
    lines.append(f"  {'crank rpm':>10}{'share':>8}{'km/L':>10}{'W':>8}{'eta_b':>8}  ok")
    for row in result.rows:
        lines.append(
            f"  {row.point.crankshaft_rpm:>10.0f}{row.point.share:>8.2f}"
            f"{row.km_per_litre:>10.1f}{row.brake_power:>8.1f}"
            f"{row.brake_efficiency:>8.3f}  {'yes' if row.feasible else row.reason}"
        )
    lines.append("")
    lines.append(
        f"  one engine of {result.engine_mass_kg:.2f} kg, structure sized at "
        f"{result.sizing_rpm:.0f} rpm,"
    )
    lines.append(f"  carrying a {result.flywheel_kg:.2f} kg flywheel for the slowest point")
    lines.append(f"  cycle range   {result.km_per_litre:.1f} km/L")
    lines.append(f"  worst / best  {result.spread:.3f}")
    return "\n".join(lines)
