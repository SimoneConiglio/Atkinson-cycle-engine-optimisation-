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
from .formulations import (
    CouplingStrength,
    compare_formulations,
    coupling_curve,
    coupling_strength,
)
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
    RANGE_DESIGN,
    REFINED_DESIGN,
)
from .robustness import ToleranceReport, tolerance_report
from .sizing import MemberSizing
from .slidercrank import (
    SliderCrank,
    evaluate_slidercrank,
    optimise_slidercrank,
    slidercrank_reliability,
)
from .synthesis import TargetMotion, fit_to_target, target_from_design, target_motion
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
    "RANGE_DESIGN",
    "REFINED_DESIGN",
    "VARIABLE_DESCRIPTIONS",
    "VARIABLE_NAMES",
    "Analysis",
    "Bounds",
    "CoupledResult",
    "CouplingStrength",
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
    "SliderCrank",
    "SolvedAnalysis",
    "TargetMotion",
    "Thermodynamics",
    "ToleranceReport",
    "Vehicle",
    "analyse",
    "compare_formulations",
    "coupling_curve",
    "coupling_strength",
    "evaluate",
    "evaluate_slidercrank",
    "fit_to_target",
    "lattice_neighbours",
    "optimise_slidercrank",
    "size_pair",
    "slidercrank_reliability",
    "solve_coupled",
    "solve_for_design",
    "speed_sweep",
    "target_from_design",
    "target_motion",
    "tolerance_report",
]

#: :mod:`exlink.minlp` is deliberately *not* re-exported here.  It needs the
#: optional ``gemseo-bilevel-outer-approximation`` plugin (``pip install
#: exlink-opt[minlp]``), and importing it from the package root would make that
#: plugin a hard dependency of everything.  Import it directly instead::
#:
#:     from exlink.minlp import solve, candidates_from_design
__version__ = "1.0.0"
