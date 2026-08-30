"""Optimization of the EX-link Atkinson-cycle engine mechanism.

A Python reconstruction of a 2015 student project (Universite de Technologie de
Compiegne, TN12 / mechanical optimization), which sized an extended-expansion
linkage for the Shell Eco-marathon: eleven design variables, three competing
objectives (mechanical efficiency, and the two envelope dimensions), and a set
of constraints most of which exist to make the problem *well posed* rather than
to express a specification.

The original was MATLAB.  Here the physics is NumPy and the optimization is
driven by `GEMSEO <https://gemseo.readthedocs.io>`_, so the same problem can be
handed to a gradient-based solver, a differential-evolution search, an augmented
Lagrangian, or NSGA-II without rewriting anything.

Quick start::

    from exlink import analyse, PUBLISHED_DESIGN
    print(analyse(PUBLISHED_DESIGN).metrics.efficiency)

    from exlink.scenarios import maximise_efficiency
    result = maximise_efficiency()
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
from .cycle import Phase, PhaseError, Phases, Thermodynamics
from .design import (
    GLOBAL_BOUNDS,
    VARIABLE_DESCRIPTIONS,
    VARIABLE_NAMES,
    Bounds,
    Design,
)
from .kinematics import Kinematics
from .loads import Loads
from .metrics import Metrics
from .model import Analysis, analyse
from .reference import PUBLISHED_DESIGN, PUBLISHED_METRICS

__all__ = [
    "DEFAULT_PENALTY",
    "DEFAULT_SPEC",
    "DEFAULT_TARGETS",
    "GLOBAL_BOUNDS",
    "PUBLISHED_DESIGN",
    "PUBLISHED_METRICS",
    "VARIABLE_DESCRIPTIONS",
    "VARIABLE_NAMES",
    "Analysis",
    "Bounds",
    "Design",
    "DesignTargets",
    "EngineSpec",
    "Kinematics",
    "Loads",
    "Metrics",
    "PenaltyValues",
    "Phase",
    "PhaseError",
    "Phases",
    "Thermodynamics",
    "analyse",
]

__version__ = "1.0.0"
