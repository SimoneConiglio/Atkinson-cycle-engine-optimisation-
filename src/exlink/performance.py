"""The whole chain, from eleven linkage dimensions to kilometres per litre.

This module does no physics of its own.  It composes the seven that do, in the
one order the couplings allow, and returns a single object holding every
intermediate so that a result can be interrogated rather than just believed.

.. code-block:: text

    Design (11 vars) + speed + gear module
        |
        v
    kinematics --> Atkinson cycle --> gas force
        |                                 |
        |          +----------------------+
        v          v
      [ dynamics <===> sizing ]        <-- MDA: neither runs first
        |
        +--> friction  -----> brake work --> brake efficiency --.
        |                                                        |
        +--> mass budget ---> engine mass ---------------------. |
                  ^                                            v v
        H, B, torque ripple ------------------------------>  RANGE

The one-way parts are ordinary function calls.  The two-way part -- sections
set masses, masses set inertia loads, loads set sections -- is the fixed point
in :mod:`exlink.coupled`, and everything downstream of it depends on it having
converged.

Why range is the right objective
---------------------------------
Every earlier formulation of this problem had to stop short of a single number.
Efficiency alone is maximised by an unboundedly large mechanism.  Adding ``H``
and ``B`` as competing objectives bounds it, but only produces a front, because
nothing in the problem says what a millimetre of height is worth in points of
efficiency.  Mass alone is minimised by a mechanism that produces no torque.

Range prices all of them, in units the competition actually scores, and the
prices are not adjustable:

* a point of brake efficiency is worth a fixed number of kilometres, through
  the fuel burnt per unit of work;
* a millimetre of ``H`` or ``B`` is worth a fixed number of grams of crankcase,
  and a gram is worth a fixed number of kilometres through rolling resistance;
* a newton-millimetre of torque ripple is worth a fixed number of grams of
  flywheel, on the same terms.

That is what makes this a multidisciplinary *design* problem rather than a
multi-objective one.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .constants import DEFAULT_SPEC, EngineSpec
from .coupled import CoupledResult, solve_coupled
from .design import Design
from .dynamics import DEFAULT_SPEED_RPM
from .friction import FrictionLosses
from .friction import losses as friction_losses
from .manufacturing import round_up_to_stock
from .mass_budget import SPEED_FLUCTUATION, MassBudget, assemble
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors
from .metrics import Metrics
from .model import Analysis, analyse
from .vehicle import RangeResult, Vehicle, best_strategy, brake_efficiency, heat_release


@dataclass(frozen=True)
class Performance:
    """Everything known about one design, at one operating point."""

    design: Design
    speed_rpm: float
    analysis: Analysis

    coupled: CoupledResult | None
    """The sizing fixed point, or ``None`` if the kinematics rejected the design."""

    friction: FrictionLosses | None
    """The loss breakdown, or ``None`` for a design that was never sized."""

    budget: MassBudget
    range: RangeResult

    heat_release: float
    """Heat added per cycle [N.mm]."""

    buildable_diameters: dict[str, float]
    """Member diameters rounded up to stock sizes [mm]."""

    stock_premium: float
    """Fractional mass added by that rounding [-]."""

    spec: EngineSpec = DEFAULT_SPEC
    """Fixed engine data, kept so the output shaft can be reported."""

    @property
    def metrics(self) -> Metrics:
        """The geometric metrics of the design."""
        return self.analysis.metrics

    @property
    def km_per_litre(self) -> float:
        """The objective [km/L]."""
        return self.range.km_per_litre

    @property
    def engine_mass_kg(self) -> float:
        """Total engine mass [kg]."""
        return self.budget.total_kg

    @property
    def output_speed_rpm(self) -> float:
        """Speed of the shaft power is taken from [rev/min].

        Twice :attr:`speed_rpm`, which is the half-speed shaft the kinematics
        are parametrised on.  The cycle spans 720 deg of the output shaft, as
        it does on a conventional four-stroke, so this is the speed to quote
        and the speed to compare at.
        """
        return self.spec.output_speed_rpm(self.speed_rpm)

    @property
    def cycles_per_minute(self) -> float:
        """Power strokes per minute [1/min], equal to :attr:`speed_rpm`."""
        return self.output_speed_rpm / self.spec.output_revolutions_per_cycle

    @property
    def output_torque(self) -> float:
        """Mean torque at the output shaft [N.mm]."""
        return self.spec.output_torque(self.analysis.metrics.mean_torque)

    @property
    def brake_power(self) -> float:
        """Brake power at this speed [W].

        One cycle per turn of ``theta_1``, hence one per two turns of the
        output shaft; reading the speed at either shaft gives the same power.
        """
        if self.friction is None:
            return 0.0
        return _brake_power(self.friction.brake_work, self.speed_rpm)

    @property
    def indicated_efficiency(self) -> float:
        """Indicated work over heat released [-]."""
        if self.friction is None or self.heat_release <= 0.0:
            return 0.0
        return self.friction.indicated_work / self.heat_release

    @property
    def brake_efficiency(self) -> float:
        """Brake work over heat released [-]."""
        if self.friction is None:
            return 0.0
        return brake_efficiency(self.friction.brake_work, self.heat_release)

    @property
    def geometrically_feasible(self) -> bool:
        """Whether the design meets the geometric constraint set.

        Kept separate from :attr:`feasible` because the two ask different
        questions and conflating them is how a violating design gets reported
        as a result: an engine can run perfectly well, and go a long way, while
        sitting outside the specification it was supposed to meet.
        """
        from .scenarios import is_feasible

        return is_feasible(self.analysis)

    @property
    def feasible(self) -> bool:
        """Whether this is a design that both works *and* meets its brief.

        Every discipline has to return something usable -- the kinematics
        closes, the sizing settles, the engine produces net work, the gears
        fit, the car meets the speed rule -- **and** the geometric constraints
        have to hold.  Leaving that last clause out reports designs that run
        but miss the stroke, the compression ratio, the top-dead-centre gap or
        the side-load limit, which is not a result.
        """
        return (
            self.analysis.valid
            and self.coupled is not None
            and self.coupled.feasible
            and self.friction is not None
            and self.friction.runs
            and self.range.feasible
            and (self.budget.gears is None or self.budget.gears.feasible)
            and self.geometrically_feasible
        )

    def reason(self) -> str:
        """Why the design is not feasible, or an empty string."""
        if not self.analysis.valid:
            return f"kinematics: {self.analysis.metrics.reason}"
        if self.coupled is None or not self.coupled.feasible:
            return "sizing did not converge on a buildable set of sections"
        if self.friction is None or not self.friction.runs:
            return "friction exceeds indicated work; the engine will not run"
        if self.budget.gears is not None and not self.budget.gears.feasible:
            return "gear pair outside its face-width or tooth-count limits"
        if not self.range.feasible:
            return f"vehicle: {self.range.reason}"
        if not self.geometrically_feasible:
            return "outside the geometric constraint set"
        return ""


def _brake_power(brake_work: float, speed_rpm: float) -> float:
    """Brake power [W] from work per revolution [N.mm] and speed [rev/min].

    One thermodynamic cycle per crankshaft revolution, and ``N.mm`` to ``J`` is
    a factor of 1000.
    """
    return brake_work / 1000.0 * (speed_rpm / 60.0)


def evaluate(
    design: Design,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    vehicle: Vehicle | None = None,
    module: float | None = None,
    teeth: int | None = None,
    samples: int = 360,
    fluctuation: float = SPEED_FLUCTUATION,
    spec: EngineSpec = DEFAULT_SPEC,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    **kwargs: Any,
) -> Performance:
    """Run the full chain for one design.

    Args:
        design: The mechanism dimensions.
        speed_rpm: Crankshaft speed [rev/min].
        vehicle: The car; a default Prototype-class entry if omitted.
        module: Gear module [mm]; the lightest workable standard one if omitted.
        teeth: Teeth on the small gear; derived from ``I`` if omitted.
        samples: Crank angles per revolution.
        fluctuation: Allowed cyclic speed fluctuation, for the flywheel.
        spec: Fixed engine data.
        material: Structural material.
        safety: Design factors.
        **kwargs: Forwarded to :func:`exlink.coupled.solve_coupled`.

    Returns:
        The full result.  Check :attr:`Performance.feasible` and
        :meth:`Performance.reason` -- an unanalysable design comes back as an
        infeasible :class:`Performance`, never as an exception, so that an
        optimizer can keep walking.
    """
    car = vehicle if vehicle is not None else Vehicle()
    analysis = analyse(design, samples=samples, spec=spec)
    if not analysis.valid:
        return _failed(design, speed_rpm, analysis, car)

    solved = analysis.require_solved()
    coupled = solve_coupled(
        solved.kinematics,
        solved.thermodynamics,
        speed_rpm=speed_rpm,
        material=material,
        safety=safety,
        spec=spec,
        **kwargs,
    )
    friction = friction_losses(coupled.loads, coupled.diameters)

    metrics = analysis.metrics
    thermo = solved.thermodynamics
    quantity = heat_release(
        spec.dead_volume,
        thermo.p_compression_end,
        thermo.p_combustion,
        spec.heat_capacity_ratio,
    )

    budget = assemble(
        coupled.loads,
        coupled.diameters,
        coupled.mass_properties.member_mass,
        coupled.piston_mass,
        metrics.height,
        metrics.width,
        metrics.expansion_stroke + metrics.compression_stroke,
        float(np.max(thermo.gauge_pressure)),
        module=module,
        teeth=teeth,
        fluctuation=fluctuation,
        spec=spec,
        material=material,
        safety=safety,
    )

    efficiency = brake_efficiency(friction.brake_work, quantity)
    power = _brake_power(friction.brake_work, speed_rpm)
    outcome = best_strategy(car, budget.total_kg, power, efficiency)

    names = list(coupled.diameters)
    continuous = np.array([coupled.diameters[n] for n in names])
    discrete = round_up_to_stock(continuous)
    from .manufacturing import stock_premium

    return Performance(
        design=design,
        speed_rpm=float(speed_rpm),
        analysis=analysis,
        coupled=coupled,
        friction=friction,
        budget=budget,
        range=outcome,
        heat_release=float(quantity),
        buildable_diameters={n: float(d) for n, d in zip(names, discrete, strict=True)},
        stock_premium=stock_premium(continuous, discrete),
        spec=spec,
    )


def _failed(
    design: Design, speed_rpm: float, analysis: Analysis, vehicle: Vehicle
) -> Performance:
    """A :class:`Performance` for a design the kinematics rejected."""
    empty_range = RangeResult(
        distance_per_joule=0.0,
        km_per_litre=0.0,
        average_speed=0.0,
        burn_fraction=0.0,
        burn_distance=0.0,
        coast_distance=0.0,
        total_mass=vehicle.mass(0.0),
        engine_mass=0.0,
        brake_power=0.0,
        brake_efficiency=0.0,
        feasible=False,
        reason=analysis.metrics.reason or "kinematically invalid",
    )
    return Performance(
        design=design,
        speed_rpm=float(speed_rpm),
        analysis=analysis,
        coupled=None,
        friction=None,
        budget=MassBudget(),
        range=empty_range,
        heat_release=0.0,
        buildable_diameters={},
        stock_premium=0.0,
    )


def speed_sweep(
    design: Design,
    speeds: tuple[float, ...] = (600.0, 800.0, 1000.0, 1250.0, 1500.0, 2000.0, 2500.0, 3000.0),
    **kwargs: Any,
) -> list[Performance]:
    """Evaluate one design across crankshaft speeds.

    Speed is the sharpest single variable in the coupled problem, and it is the
    one the geometric formulation could not see at all.  It pulls two ways:

    * **up**, because the flywheel inertia needed for a given speed
      fluctuation goes as ``1 / omega^2``, so a faster engine needs a
      dramatically lighter flywheel -- often the heaviest item in the budget;
    * **down**, because every inertia load goes as ``omega^2``, and the
      sections that carry them, and hence the structural mass, follow.

    The optimum sits where those two cross, and it is not near either end.

    Args:
        design: The mechanism dimensions.
        speeds: Crankshaft speeds to evaluate [rev/min].
        **kwargs: Forwarded to :func:`evaluate`.

    Returns:
        One :class:`Performance` per speed, in the order given.
    """
    return [evaluate(design, speed_rpm=speed, **kwargs) for speed in speeds]
