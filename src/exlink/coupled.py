"""The sizing / dynamics fixed point, and why it is a fixed point at all.

Adding inertia to the load analysis closes a loop that the quasi-static study
never had:

.. code-block:: text

    section diameters --> member masses --> inertia forces --> internal loads
            ^                                                        |
            +--------------- sizing against yield, ------------------+
                             fatigue and buckling

Neither half can be evaluated first.  The loads cannot be computed without the
masses, and the masses cannot be chosen without the loads.  That is a genuine
multidisciplinary coupling and it has to be *solved*, not sequenced.

Does it converge?
-----------------
Scaling tells you what to expect.  A bending-critical member needs
``d ~ F^(1/3)``, so its mass goes as ``m ~ d^2 ~ F^(2/3)``; the inertia force it
then generates is ``F ~ m a``.  Composing, ``m ~ (C a) m^(2/3)``: the loop gain
is sub-linear, so a fixed point exists at ``m = (C a)^3`` and plain iteration
reaches it -- but slowly, and with a mass that grows as the *cube* of the
acceleration level.  Under-relaxation is provided for the stiff cases, and the
cubic sensitivity is why the answer is so sensitive to engine speed.

The result is that speed, which does not enter the geometric problem at all,
becomes one of the strongest drivers here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .constants import DEFAULT_SPEC, EngineSpec
from .cycle import Thermodynamics
from .design import Design
from .dynamics import (
    DEFAULT_SPEED_RPM,
    MEMBER_NAMES,
    DynamicLoads,
    MassProperties,
    mass_properties,
    rpm_to_rad_per_s,
)
from .dynamics import solve as solve_dynamics
from .kinematics import Kinematics
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors
from .sizing import MAX_DIAMETER, MemberSizing, piston_mass, size_members

INITIAL_DIAMETER = 8.0
"""Diameter every member starts from [mm]."""

DEFAULT_TOLERANCE = 1.0e-6
"""Convergence tolerance on the largest diameter change between sweeps [mm]."""

DEFAULT_MAX_ITERATIONS = 200
DEFAULT_RELAXATION = 1.0
"""Under-relaxation factor; below 1 damps the loop at the cost of more sweeps."""


@dataclass(frozen=True)
class CoupledResult:
    """Outcome of the sizing / dynamics fixed point."""

    diameters: dict[str, float]
    """Converged section diameter of each member [mm]."""

    sizing: dict[str, MemberSizing]
    loads: DynamicLoads
    mass_properties: MassProperties
    piston_crown_thickness: float
    piston_mass: float
    """[tonne]"""

    speed: float
    """Crankshaft speed ``Omega`` [rad/s]."""

    iterations: int
    residual: float
    """Largest diameter change in the final sweep [mm]."""

    converged: bool
    history: list[float]
    """Residual after each sweep, for diagnosing a stiff or runaway case."""

    @property
    def total_mass_kg(self) -> float:
        """Total moving mass, piston included [kg]."""
        return self.mass_properties.total_mass_kg + 1000.0 * self.piston_mass

    @property
    def saturated(self) -> bool:
        """Whether any member was driven to the diameter ceiling.

        A saturated member means the mechanism cannot be built to survive its
        own inertia at this speed -- the loop has run away rather than settled.
        """
        return any(d >= MAX_DIAMETER - 1.0e-6 for d in self.diameters.values())

    @property
    def feasible(self) -> bool:
        """Whether the fixed point settled on a buildable set of sections."""
        return self.converged and not self.saturated

    @property
    def peak_bearing_load(self) -> float:
        """Largest crankshaft bearing reaction over the revolution [N]."""
        return float(np.max(np.linalg.norm(self.loads.reaction["R1"], axis=1)))

    def worst_utilisation(self) -> tuple[str, float]:
        """The member closest to its allowable, and by how much.

        Sizing drives every member to a utilisation of 1, so this is a check
        that the solve did its job rather than a design margin.
        """
        worst_name, worst = "", 0.0
        for name, item in self.sizing.items():
            value = max(
                item.static_utilisation,
                item.fatigue_utilisation,
                item.buckling_utilisation,
            )
            if value > worst:
                worst_name, worst = name, value
        return worst_name, worst


def solve_coupled(
    kinematics: Kinematics,
    thermodynamics: Thermodynamics,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    spec: EngineSpec = DEFAULT_SPEC,
    initial_diameters: dict[str, float] | None = None,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    relaxation: float = DEFAULT_RELAXATION,
) -> CoupledResult:
    """Solve the sizing / dynamics coupling by relaxed fixed-point iteration.

    This is the reference implementation: a plain Gauss-Seidel sweep, easy to
    read and to check against.  :mod:`exlink.disciplines` exposes the same
    coupling to GEMSEO as two disciplines joined by an MDA, and
    ``tests/test_coupled.py`` checks the two agree.

    Args:
        kinematics: A solved mechanism.
        thermodynamics: Its solved cycle, supplying the gas load.
        speed_rpm: Crankshaft speed [rev/min].  Zero recovers purely
            quasi-static sizing, with no coupling at all.
        material: The material.
        safety: The design factors.
        spec: Fixed engine data.
        initial_diameters: Starting sections [mm]; a warm start from a nearby
            design cuts the sweep count sharply.
        tolerance: Convergence tolerance on the diameter change [mm].
        max_iterations: Sweep limit.
        relaxation: Under-relaxation in ``(0, 1]``.

    Returns:
        The converged sections and the load case that goes with them.  Check
        :attr:`CoupledResult.feasible` -- a design whose inertia outruns its
        structure will exit either unconverged or saturated at the ceiling.

    Raises:
        ValueError: If ``relaxation`` is outside ``(0, 1]``.
    """
    if not 0.0 < relaxation <= 1.0:
        msg = "relaxation must lie in (0, 1]"
        raise ValueError(msg)

    speed = rpm_to_rad_per_s(speed_rpm)
    crown, piston = piston_mass(thermodynamics, material, safety, spec)

    diameters = dict(initial_diameters or dict.fromkeys(MEMBER_NAMES, INITIAL_DIAMETER))
    history: list[float] = []
    residual = float("inf")
    properties = mass_properties(kinematics, diameters, material.density, piston, spec)
    loads = solve_dynamics(kinematics, thermodynamics.piston_force, properties, speed, spec)
    sizing: dict[str, MemberSizing] = {}

    # ``iteration`` is read after the loop, to report how many sweeps it took.
    iteration = 0
    for iteration in range(1, max_iterations + 1):  # noqa: B007
        properties = mass_properties(kinematics, diameters, material.density, piston, spec)
        loads = solve_dynamics(kinematics, thermodynamics.piston_force, properties, speed, spec)
        sizing = size_members(loads, material, safety)

        updated = {
            name: (1.0 - relaxation) * diameters[name] + relaxation * sizing[name].diameter
            for name in MEMBER_NAMES
        }
        residual = max(abs(updated[n] - diameters[n]) for n in MEMBER_NAMES)
        diameters = updated
        history.append(residual)
        if residual <= tolerance:
            break

    properties = mass_properties(kinematics, diameters, material.density, piston, spec)
    return CoupledResult(
        diameters=diameters,
        sizing=sizing,
        loads=loads,
        mass_properties=properties,
        piston_crown_thickness=crown,
        piston_mass=piston,
        speed=speed,
        iterations=iteration,
        residual=residual,
        converged=residual <= tolerance,
        history=history,
    )


def solve_for_design(
    design: Design,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    samples: int = 360,
    **kwargs: Any,
) -> CoupledResult:
    """Convenience wrapper: analyse a design, then solve its sizing coupling.

    Args:
        design: The mechanism dimensions.
        speed_rpm: Crankshaft speed [rev/min].
        samples: Crank angles per revolution.
        **kwargs: Forwarded to :func:`solve_coupled`.

    Returns:
        The coupled result.

    Raises:
        ValueError: If the design cannot be analysed kinematically at all.
    """
    from .model import analyse

    analysis = analyse(design, samples=samples)
    if not analysis.valid:
        msg = f"cannot size an unanalysable design: {analysis.metrics.reason}"
        raise ValueError(msg)
    solved = analysis.require_solved()
    return solve_coupled(
        solved.kinematics, solved.thermodynamics, speed_rpm=speed_rpm, **kwargs
    )
