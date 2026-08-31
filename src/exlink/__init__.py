"""Optimization of the EX-link Atkinson-cycle engine mechanism.

A multidisciplinary optimization of an extended-expansion (Atkinson) linkage
for a Shell Eco-marathon single-cylinder engine: eleven design variables, three competing
objectives (mechanical efficiency, and the two envelope dimensions), and a set
of constraints most of which exist to make the problem *well posed* rather than
to express a specification.

The original was MATLAB.  Here the physics is NumPy and the optimization is
driven by `GEMSEO <https://gemseo.readthedocs.io>`_, so the same problem can be
handed to a gradient-based solver, a differential-evolution search, an augmented
Lagrangian, or NSGA-II without rewriting anything.

The geometric problem above stops before sizing the parts, because their masses
are not known until they have a shape.  :mod:`exlink.dynamics`, :mod:`exlink.sizing` and
:mod:`exlink.coupled` carry out that next iteration, which closes a loop:
sections set the masses, the masses set the inertia loads, and the loads set the
sections.  That coupling has to be solved, and GEMSEO solves it with an MDA.

Quick start::

    from exlink import analyse, PUBLISHED_DESIGN
    print(analyse(PUBLISHED_DESIGN).metrics.efficiency)

    from exlink.scenarios import maximise_efficiency
    result = maximise_efficiency()

    # size the parts, with inertia in the load path
    from exlink import solve_for_design
    sized = solve_for_design(PUBLISHED_DESIGN, speed_rpm=1000.0)
    print(sized.total_mass_kg, sized.diameters)
"""

from __future__ import annotations

from .constants import (
    DEFAULT_PENALTY,
    DEFAULT_SPEC,
    DEFAULT_TARGETS,
    DesignTargets,
    EngineSpec,
    PenaltyValues,
)
from .coupled import CoupledResult, solve_coupled, solve_for_design
from .cycle import Phase, PhaseError, Phases, Thermodynamics
from .design import (
    GLOBAL_BOUNDS,
    VARIABLE_DESCRIPTIONS,
    VARIABLE_NAMES,
    Bounds,
    Design,
)
from .dynamics import DEFAULT_SPEED_RPM, MEMBER_NAMES, DynamicLoads, MassProperties
from .friction import FrictionLosses
from .gears import GearPair, lattice_neighbours, size_pair
from .kinematics import Kinematics
from .loads import Loads
from .mass_budget import MassBudget
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors
from .metrics import Metrics
from .model import Analysis, SolvedAnalysis, analyse
from .performance import Performance, evaluate, speed_sweep
from .reference import (
    COUPLED_DESIGN,
    GRADIENT_DESIGN,
    PUBLISHED_DESIGN,
    PUBLISHED_METRICS,
    REFINED_DESIGN,
)
from .sizing import MemberSizing
from .vehicle import RangeResult, Vehicle

__all__ = [
    "COUPLED_DESIGN",
    "DEFAULT_MATERIAL",
    "DEFAULT_PENALTY",
    "DEFAULT_SAFETY",
    "DEFAULT_SPEC",
    "DEFAULT_SPEED_RPM",
    "DEFAULT_TARGETS",
    "GLOBAL_BOUNDS",
    "GRADIENT_DESIGN",
    "MEMBER_NAMES",
    "PUBLISHED_DESIGN",
    "PUBLISHED_METRICS",
    "REFINED_DESIGN",
    "VARIABLE_DESCRIPTIONS",
    "VARIABLE_NAMES",
    "Analysis",
    "Bounds",
    "CoupledResult",
    "Design",
    "DesignTargets",
    "DynamicLoads",
    "EngineSpec",
    "FrictionLosses",
    "GearPair",
    "Kinematics",
    "Loads",
    "MassBudget",
    "MassProperties",
    "Material",
    "MemberSizing",
    "Metrics",
    "PenaltyValues",
    "Performance",
    "Phase",
    "PhaseError",
    "Phases",
    "RangeResult",
    "SafetyFactors",
    "SolvedAnalysis",
    "Thermodynamics",
    "Vehicle",
    "analyse",
    "evaluate",
    "lattice_neighbours",
    "size_pair",
    "solve_coupled",
    "solve_for_design",
    "speed_sweep",
]

__version__ = "1.0.0"
