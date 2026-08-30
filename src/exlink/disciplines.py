"""GEMSEO wrappers around the mechanism analysis.

The 2015 study was written in MATLAB, with the penalty function, the design
space and the algorithm loop all hand-rolled.  Here GEMSEO owns the problem
formulation: :class:`ExlinkDiscipline` exposes the analysis as a discipline with
a named grammar, and the scenarios in :mod:`exlink.scenarios` attach objectives
and constraints to it declaratively.  Every algorithm in GEMSEO's catalogue --
gradient-based, derivative-free, evolutionary, multi-objective -- then applies
to the same problem without touching the physics.

Two disciplines are provided:

:class:`ExlinkDiscipline`
    Takes the eleven design variables as separate scalar inputs and returns
    every objective and constraint measure as a separate scalar output.  This is
    the natural form for GEMSEO and the one the scenarios use.

:class:`PenalisedExlinkDiscipline`
    Adds the report's external penalty function ``F(X)`` as one extra output, so
    that the historical "penalise, then run an unconstrained solver" workflow
    can be reproduced as-is.
"""

from __future__ import annotations

from typing import Any, ClassVar

import numpy as np
from gemseo.core.discipline import Discipline
from gemseo.typing import StrKeyMapping

from .constants import (
    DEFAULT_PENALTY,
    DEFAULT_SPEC,
    DEFAULT_TARGETS,
    DesignTargets,
    EngineSpec,
    PenaltyValues,
)
from .design import VARIABLE_NAMES, Design
from .kinematics import DEFAULT_SAMPLES
from .model import (
    EQUALITY_NAMES,
    INEQUALITY_NAMES,
    Analysis,
    analyse,
    equality_constraints,
    inequality_constraints,
)
from .reference import PUBLISHED_DESIGN

#: Outputs produced by :class:`ExlinkDiscipline`, in a stable order.
OUTPUT_NAMES: tuple[str, ...] = (
    "neg_efficiency",
    "efficiency",
    "height",
    "width",
    *INEQUALITY_NAMES,
    *EQUALITY_NAMES,
    "expansion_stroke",
    "compression_ratio",
    "rod_angle",
    "compatibility",
    "tdc_gap",
    "clearance",
    "side_load_ratio",
    "mean_torque",
    "valid",
)


class ExlinkDiscipline(Discipline):
    """The EX-link mechanism as a GEMSEO discipline.

    Inputs are the eleven design variables of :data:`exlink.design.VARIABLE_NAMES`,
    each a scalar array.  Outputs are the three objectives, the seven constraint
    residuals and a handful of diagnostics; see :data:`OUTPUT_NAMES`.

    Constraint outputs follow GEMSEO's convention of "feasible when ``<= 0``"
    for inequalities and "feasible when ``== 0``" for equalities, so they can be
    handed straight to :meth:`~gemseo.scenarios.base_scenario.BaseScenario.add_constraint`.

    Args:
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        penalty: Values substituted for an unanalysable design.
        name: Discipline name.
    """

    auto_detect_grammar_files: ClassVar[bool] = False

    def __init__(
        self,
        samples: int = DEFAULT_SAMPLES,
        spec: EngineSpec = DEFAULT_SPEC,
        targets: DesignTargets = DEFAULT_TARGETS,
        penalty: PenaltyValues = DEFAULT_PENALTY,
        name: str = "",
    ) -> None:
        super().__init__(name=name)
        self.samples = samples
        self.spec = spec
        self.targets = targets
        self.penalty = penalty
        self.input_grammar.update_from_names(VARIABLE_NAMES)
        self.output_grammar.update_from_names(OUTPUT_NAMES)
        self.default_input_data = PUBLISHED_DESIGN.to_mapping()

    def analyse_design(self, design: Design) -> Analysis:
        """Run the underlying analysis with this discipline's settings."""
        return analyse(
            design,
            samples=self.samples,
            spec=self.spec,
            targets=self.targets,
            penalty=self.penalty,
        )

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        design = Design.from_mapping(dict(input_data))
        analysis = self.analyse_design(design)
        return to_output_data(analysis, self.targets)


def to_output_data(
    analysis: Analysis, targets: DesignTargets = DEFAULT_TARGETS
) -> dict[str, np.ndarray]:
    """Flatten an :class:`~exlink.model.Analysis` into GEMSEO output data."""
    metrics = analysis.metrics
    inequality = inequality_constraints(analysis, targets)
    equality = equality_constraints(analysis, targets)
    values: dict[str, float] = {
        "neg_efficiency": -metrics.efficiency,
        "efficiency": metrics.efficiency,
        "height": metrics.height,
        "width": metrics.width,
        "expansion_stroke": metrics.expansion_stroke,
        "compression_ratio": metrics.compression_ratio,
        "rod_angle": metrics.rod_angle,
        "compatibility": metrics.compatibility,
        "tdc_gap": metrics.tdc_gap,
        "clearance": metrics.clearance,
        "side_load_ratio": min(metrics.side_load_ratio, 1.0e3),
        "mean_torque": metrics.mean_torque,
        "valid": float(metrics.valid),
    }
    values.update(dict(zip(INEQUALITY_NAMES, inequality, strict=True)))
    values.update(dict(zip(EQUALITY_NAMES, equality, strict=True)))
    return {name: np.array([values[name]], dtype=float) for name in OUTPUT_NAMES}


class PenalisedExlinkDiscipline(ExlinkDiscipline):
    """Adds the report's external penalty function as an output.

    The report converts the constrained problem into an unconstrained one via

    .. math::
        F(X) = -\\eta(X) + \\frac{1}{r^2}\\left(
            c_{eq}^T c_{eq} + \\langle c \\rangle^T \\langle c \\rangle \\right),

    where ``<c>`` keeps only the violated inequalities and ``0 < r < 1`` is the
    penalty parameter.  Smaller ``r`` means a more accurate but worse
    conditioned problem -- the trade-off the report describes, and the reason it
    finishes with an augmented Lagrangian instead.

    The size objectives are handled the way the report handles them when it
    sweeps a Pareto front by hand: as moving limits ``H <= h_max``,
    ``B <= b_max`` folded into the penalty.

    Args:
        penalty_parameter: ``r``, in ``(0, 1]``.
        max_height: Moving limit on ``H`` [mm]; ``inf`` disables it.
        max_width: Moving limit on ``B`` [mm]; ``inf`` disables it.
        **kwargs: Forwarded to :class:`ExlinkDiscipline`.
    """

    def __init__(
        self,
        penalty_parameter: float = 0.1,
        max_height: float = float("inf"),
        max_width: float = float("inf"),
        **kwargs: Any,
    ) -> None:
        if not 0.0 < penalty_parameter <= 1.0:
            msg = "the penalty parameter r must lie in (0, 1]"
            raise ValueError(msg)
        super().__init__(**kwargs)
        self.penalty_parameter = penalty_parameter
        self.max_height = max_height
        self.max_width = max_width
        self.output_grammar.update_from_names(["penalised_objective"])

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        design = Design.from_mapping(dict(input_data))
        analysis = self.analyse_design(design)
        output = to_output_data(analysis, self.targets)
        output["penalised_objective"] = np.array(
            [self.penalised_objective(analysis)], dtype=float
        )
        return output

    def penalised_objective(self, analysis: Analysis) -> float:
        """Evaluate ``F(X)`` for an analysed design."""
        metrics = analysis.metrics
        inequality = list(inequality_constraints(analysis, self.targets))
        equality = list(equality_constraints(analysis, self.targets))
        if np.isfinite(self.max_height):
            inequality.append(metrics.height - self.max_height)
        if np.isfinite(self.max_width):
            inequality.append(metrics.width - self.max_width)

        violated = np.maximum(np.asarray(inequality, dtype=float), 0.0)
        residual = np.asarray(equality, dtype=float)
        penalty = float(residual @ residual + violated @ violated)
        return float(-metrics.efficiency + penalty / self.penalty_parameter**2)
