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
from .coupled import INITIAL_DIAMETER
from .design import VARIABLE_NAMES, Design
from .dynamics import (
    DEFAULT_SPEED_RPM,
    MEMBER_NAMES,
    mass_properties,
    rpm_to_rad_per_s,
)
from .dynamics import solve as solve_dynamics
from .kinematics import DEFAULT_SAMPLES
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors
from .model import (
    EQUALITY_NAMES,
    INEQUALITY_NAMES,
    Analysis,
    analyse,
    equality_constraints,
    inequality_constraints,
)
from .reference import PUBLISHED_DESIGN
from .sizing import (
    MAX_DIAMETER,
    MEMBER_IS_SLENDER,
    STATIONS,
    member_lengths,
    member_loads,
    piston_mass,
    size_from_arrays,
)

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


# =============================================================================
# The coupled sizing problem
#
# Everything above is the report's own, one-way problem: geometry in, efficiency
# and envelope out.  What follows adds the iteration the report deferred, and it
# is genuinely two-way.  :class:`DynamicsDiscipline` needs the member sections
# to know the inertia forces; :class:`StructureDiscipline` needs the resulting
# internal loads to choose those sections.  Neither can run first, so an MDA has
# to resolve them -- see :func:`exlink.scenarios.build_coupled_scenario`.
# =============================================================================

COUPLING_AXIAL = "member_axial"
"""Internal axial force of every member, flattened [N]."""

COUPLING_BENDING = "member_bending"
"""Internal bending moment of every member, flattened [N.mm]."""

COUPLING_DIAMETERS = "diameters"
"""Section diameter of every member, in the order of ``MEMBER_NAMES`` [mm]."""

COUPLED_SAMPLES = 360
"""Crank angles per revolution used inside the coupled disciplines.

This sets the size of the coupling vector -- one value per member, per angle,
per station -- so there is a real incentive to keep it small.  It is set by
``g`` rather than by the loads, which are smooth and would be happy with far
fewer.

``g`` is the difference between two nearly equal maxima of ``lambda``, so its
*absolute* error is what matters against a 0.01 mm bound, and on a coarse grid
that error is not small: for the refined reference, ``g`` reads 0.0086 mm at
180 angles against a converged 0.0060 mm -- 44 % high.  Optimizing against a
constraint measured that badly lands the design outside the real one.  At 360
angles the error is 3e-4 mm, or 3 % of the bound; 720 brings it to 1e-5 mm at
twice the cost.  See ``tests/test_coupled.py``.
"""


def _initial_diameters() -> np.ndarray:
    """Starting sections for the MDA [mm]."""
    return np.full(len(MEMBER_NAMES), INITIAL_DIAMETER)


class DynamicsDiscipline(Discipline):
    """Loads on the mechanism at speed, given the member sections.

    Takes the eleven design variables, the section diameters and the crankshaft
    speed; returns the internal load history of every member together with the
    output torque and the bearing loads.

    Args:
        speed_rpm: Crankshaft speed [rev/min].
        samples: Crank angles per revolution.
        stations: Sections evaluated along each member.
        material: The material, for its density.
        safety: The design factors, used only to size the piston crown.
        spec: Fixed engine data.
        name: Discipline name.
    """

    auto_detect_grammar_files: ClassVar[bool] = False

    def __init__(
        self,
        speed_rpm: float = DEFAULT_SPEED_RPM,
        samples: int = COUPLED_SAMPLES,
        stations: int = STATIONS,
        material: Material = DEFAULT_MATERIAL,
        safety: SafetyFactors = DEFAULT_SAFETY,
        spec: EngineSpec = DEFAULT_SPEC,
        name: str = "",
    ) -> None:
        super().__init__(name=name)
        self.speed_rpm = speed_rpm
        self.samples = samples
        self.stations = stations
        self.material = material
        self.safety = safety
        self.spec = spec

        self.input_grammar.update_from_names([*VARIABLE_NAMES, COUPLING_DIAMETERS])
        self.output_grammar.update_from_names(
            [
                COUPLING_AXIAL,
                COUPLING_BENDING,
                # Named apart from ExlinkDiscipline's quasi-static ``mean_torque``
                # so both can appear in one problem. The two are provably equal:
                # at constant speed the inertia forces do no net work over a
                # closed cycle, so they reshape the torque curve without moving
                # its mean. ``tests/test_dynamics.py`` pins that.
                "dynamic_mean_torque",
                "peak_bearing_load",
                "conditioning",
                "analysable",
                "piston_mass",
            ]
        )
        self.default_input_data = {
            **PUBLISHED_DESIGN.to_mapping(),
            COUPLING_DIAMETERS: _initial_diameters(),
        }

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        data = dict(input_data)
        design = Design.from_mapping(data)
        diameters = np.asarray(data[COUPLING_DIAMETERS], dtype=float).ravel()
        shape = (len(MEMBER_NAMES), self.samples, self.stations)

        analysis = analyse(design, samples=self.samples, spec=self.spec)
        if not analysis.valid:
            # A design the kinematics cannot even close has no load case.  Return
            # zeros and flag it: the structural discipline then sizes to its floor
            # and the optimizer is steered by the constraints it can still see.
            return {
                COUPLING_AXIAL: np.zeros(shape).ravel(),
                COUPLING_BENDING: np.zeros(shape).ravel(),
                "dynamic_mean_torque": np.zeros(1),
                "peak_bearing_load": np.zeros(1),
                "conditioning": np.array([1.0e12]),
                "analysable": np.zeros(1),
                "piston_mass": np.zeros(1),
            }

        solved = analysis.require_solved()
        _, piston = piston_mass(solved.thermodynamics, self.material, self.safety, self.spec)
        properties = mass_properties(
            solved.kinematics,
            dict(zip(MEMBER_NAMES, diameters, strict=True)),
            self.material.density,
            piston,
            self.spec,
        )
        loads = solve_dynamics(
            solved.kinematics,
            solved.thermodynamics.piston_force,
            properties,
            rpm_to_rad_per_s(self.speed_rpm),
            self.spec,
        )
        per_member = member_loads(loads, stations=self.stations)
        axial = np.stack([per_member[n][0] for n in MEMBER_NAMES])
        bending = np.stack([per_member[n][1] for n in MEMBER_NAMES])
        bearing = float(np.max(np.linalg.norm(loads.reaction["R1"], axis=1)))
        return {
            COUPLING_AXIAL: axial.ravel(),
            COUPLING_BENDING: bending.ravel(),
            "dynamic_mean_torque": np.array([loads.mean_torque]),
            "peak_bearing_load": np.array([bearing]),
            "conditioning": np.array([loads.conditioning]),
            "analysable": np.ones(1),
            "piston_mass": np.array([1000.0 * piston]),
        }


class StructureDiscipline(Discipline):
    """Section sizes that survive a given load history.

    Takes the design variables and the internal load histories; returns the
    diameter each member needs to satisfy yield, fatigue and buckling, together
    with the mass that follows.

    Args:
        samples: Crank angles per revolution, matching the dynamics discipline.
        stations: Sections evaluated along each member.
        material: The material.
        safety: The design factors.
        spec: Fixed engine data.
        name: Discipline name.
    """

    auto_detect_grammar_files: ClassVar[bool] = False

    def __init__(
        self,
        samples: int = COUPLED_SAMPLES,
        stations: int = STATIONS,
        material: Material = DEFAULT_MATERIAL,
        safety: SafetyFactors = DEFAULT_SAFETY,
        spec: EngineSpec = DEFAULT_SPEC,
        name: str = "",
    ) -> None:
        super().__init__(name=name)
        self.samples = samples
        self.stations = stations
        self.material = material
        self.safety = safety
        self.spec = spec

        self.input_grammar.update_from_names(
            [*VARIABLE_NAMES, COUPLING_AXIAL, COUPLING_BENDING, "piston_mass"]
        )
        self.output_grammar.update_from_names(
            [
                COUPLING_DIAMETERS,
                "member_mass",
                "structural_mass",
                "total_mass",
                "max_utilisation",
                "saturation_margin",
                "slenderness_margin",
            ]
        )
        shape = (len(MEMBER_NAMES), samples, stations)
        self.default_input_data = {
            **PUBLISHED_DESIGN.to_mapping(),
            COUPLING_AXIAL: np.zeros(shape).ravel(),
            COUPLING_BENDING: np.zeros(shape).ravel(),
            "piston_mass": np.zeros(1),
        }

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        data = dict(input_data)
        design = Design.from_mapping(data)
        shape = (len(MEMBER_NAMES), self.samples, self.stations)
        axial = np.asarray(data[COUPLING_AXIAL], dtype=float).reshape(shape)
        bending = np.asarray(data[COUPLING_BENDING], dtype=float).reshape(shape)
        lengths = member_lengths(design)

        sizing = size_from_arrays(axial, bending, lengths, self.material, self.safety)
        diameters = np.array([sizing[n].diameter for n in MEMBER_NAMES])
        masses = np.array([sizing[n].mass for n in MEMBER_NAMES])
        utilisation = max(
            max(s.static_utilisation, s.fatigue_utilisation, s.buckling_utilisation)
            for s in sizing.values()
        )
        structural = 1000.0 * float(masses.sum())
        piston = float(np.ravel(data["piston_mass"])[0])
        return {
            COUPLING_DIAMETERS: diameters,
            "member_mass": masses,
            "structural_mass": np.array([structural]),
            "total_mass": np.array([structural + piston]),
            "max_utilisation": np.array([utilisation]),
            # Negative while every member stays clear of the diameter ceiling;
            # positive means the loop has run away and no section is thick enough.
            "saturation_margin": np.array([float(diameters.max()) - 0.98 * MAX_DIAMETER]),
            # A connecting link thicker than a third of its own length is no
            # longer a rod and its beam idealisation has stopped being credible.
            # Crank throws are exempt -- see MEMBER_IS_SLENDER.
            "slenderness_margin": np.array(
                [
                    float(
                        np.max(
                            (diameters / lengths)[MEMBER_IS_SLENDER],
                            initial=0.0,
                        )
                    )
                    - 0.34
                ]
            ),
        }
