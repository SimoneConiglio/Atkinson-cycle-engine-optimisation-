"""GEMSEO wrappers around the mechanism analysis.

GEMSEO owns the problem formulation, rather than the penalty function, the
design space and the algorithm loop being hand-rolled: :class:`ExlinkDiscipline`
exposes the analysis as a discipline with a named grammar, and the scenarios in
:mod:`exlink.scenarios` attach objectives and constraints to it declaratively.
Every algorithm in GEMSEO's catalogue -- gradient-based, derivative-free,
evolutionary, multi-objective -- then applies to the same problem without
touching the physics.

Two disciplines are provided:

:class:`ExlinkDiscipline`
    Takes the eleven design variables as separate scalar inputs and returns
    every objective and constraint measure as a separate scalar output.  This is
    the natural form for GEMSEO and the one the scenarios use.

:class:`PenalisedExlinkDiscipline`
    Adds an external penalty function ``F(X)`` as one extra output, so that a
    "penalise, then run an unconstrained solver" workflow can be run as-is.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
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
from .dynamics_jacobian import (
    DIAMETER_SLICE,
    N_PARAMETERS,
    coupled_jacobian,
    member_length_jacobian,
    sizing_jacobian,
)
from .friction import losses as friction_losses
from .gears import MAX_WIDTH_FACTOR
from .jacobian import kinematic_jacobian, metric_jacobian
from .kinematics import DEFAULT_SAMPLES
from .mass_budget import SPEED_FLUCTUATION, assemble
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
from .vehicle import Vehicle, best_strategy, brake_efficiency, heat_release

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

    def _compute_jacobian(
        self,
        input_names: Iterable[str] = (),
        output_names: Iterable[str] = (),
    ) -> None:
        """Differentiate the analysis with respect to the design vector.

        The tight constraints -- the two equalities and ``W``, ``mra``, ``g``,
        ``gamma`` -- are differentiated exactly, by chaining through the closed
        forms and applying the envelope theorem at each extremum over the crank
        angle.  That matters more than it sounds: those metrics are *maxima*
        over the revolution, so the sample attaining them switches as the design
        moves and a difference quotient taken across the switch is simply wrong.
        Measured on ``gamma`` at the reference design, a step of 1e-4 mm gives a
        gradient 25 % off, because the argmax of the side load moves one sample.

        The remaining outputs -- efficiency, the two envelope dimensions and the
        clearance -- are filled by central differences.  None is tight, all are
        smooth in the design, and efficiency in particular would need the crank
        angle of top dead centre differentiated too, because the combustion
        pressure jump puts moving-boundary terms in its integral.
        """
        self._init_jacobian(input_names, output_names)
        design = Design.from_mapping(dict(self.io.data))
        analysis = self.analyse_design(design)

        exact: dict[str, np.ndarray] = {}
        if analysis.valid:
            derivatives = kinematic_jacobian(design, analysis.kinematics, self.spec)
            exact = metric_jacobian(design, analysis, derivatives, self.spec)

        wanted = set(output_names) or set(self.jac)
        approximate = [name for name in wanted if name not in exact]
        numerical = self._difference_outputs(design, approximate) if approximate else {}

        for output in self.jac:
            row = exact.get(output)
            if row is None:
                row = numerical.get(output)
            if row is None:
                row = np.zeros(len(VARIABLE_NAMES))
            for index, variable in enumerate(VARIABLE_NAMES):
                if variable in self.jac[output]:
                    self.jac[output][variable] = np.array([[row[index]]])

    def _difference_outputs(
        self, design: Design, outputs: Sequence[str], step: float = 1.0e-6
    ) -> dict[str, np.ndarray]:
        """Central differences for the outputs left to numerical treatment.

        Args:
            design: The point to differentiate at.
            outputs: Output names needing a gradient.
            step: Relative step; scaled by each variable's magnitude.

        Returns:
            ``{output name: gradient}``, each shaped ``(11,)``.
        """
        base = design.to_array()
        gradients = {name: np.zeros(len(VARIABLE_NAMES)) for name in outputs}
        for index in range(len(VARIABLE_NAMES)):
            offset = step * max(abs(base[index]), 1.0)
            forward, backward = base.copy(), base.copy()
            forward[index] += offset
            backward[index] -= offset
            plus = to_output_data(self.analyse_design(Design.from_array(forward)), self.targets)
            minus = to_output_data(
                self.analyse_design(Design.from_array(backward)), self.targets
            )
            for name in outputs:
                gradients[name][index] = (float(plus[name][0]) - float(minus[name][0])) / (
                    2.0 * offset
                )
        return gradients


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
    """Adds an external penalty function as an output.

    This converts the constrained problem into an unconstrained one via

    .. math::
        F(X) = -\\eta(X) + \\frac{1}{r^2}\\left(
            c_{eq}^T c_{eq} + \\langle c \\rangle^T \\langle c \\rangle \\right),

    where ``<c>`` keeps only the violated inequalities and ``0 < r < 1`` is the
    penalty parameter.  Smaller ``r`` means a more accurate but worse
    conditioned problem -- the trade-off that motivates finishing with an
    augmented Lagrangian instead.

    The size objectives are handled the way a hand-swept Pareto front handles
    them: as moving limits ``H <= h_max``,
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
# Everything above is the one-way, geometric problem: geometry in, efficiency
# and envelope out.  What follows adds the structural iteration, and it is
# genuinely two-way.  :class:`DynamicsDiscipline` needs the member sections
# to know the inertia forces; :class:`StructureDiscipline` needs the resulting
# internal loads to choose those sections.  Neither can run first, so an MDA has
# to resolve them -- see :func:`exlink.scenarios.build_coupled_scenario`.
# =============================================================================

RANGE_OUTPUTS: tuple[str, ...] = (
    "km_per_litre",
    "neg_range",
    "engine_mass",
    "flywheel_mass",
    "crankcase_mass",
    "brake_efficiency",
    "mechanical_efficiency",
    "brake_power",
    "runs_margin",
    "gear_margin",
    "runs_violation",
    "gear_violation",
)
"""Everything :class:`RangeDiscipline` reports, in a stable order.

The two margins appear twice, once in each sign convention.  ``*_margin`` is
positive when satisfied, which reads naturally and is what the single-level
scenarios attach with ``positive=True``.  ``*_violation`` is the negation,
positive when violated, matching the coupled margins.  The bi-level formulation
needs the second form because its scenario adapter identifies a constraint by
its *output* name, and ``positive=True`` renames it to ``-runs_margin``, which
is not an output of anything.
"""

#: What the range chain returns for a design the kinematics cannot close.
RANGE_PENALTY: dict[str, float] = {
    "km_per_litre": 0.0,
    "neg_range": 0.0,
    "engine_mass": 1000.0,
    "flywheel_mass": 1000.0,
    "crankcase_mass": 1000.0,
    "brake_efficiency": 0.0,
    "mechanical_efficiency": 0.0,
    "brake_power": 0.0,
    "runs_margin": -1.0,
    "gear_margin": -1.0,
    "runs_violation": 1.0,
    "gear_violation": 1.0,
}


COUPLING_AXIAL = "member_axial"
"""Internal axial force of every member, flattened [N]."""

COUPLING_BENDING = "member_bending"
"""Internal bending moment of every member, flattened [N.mm]."""

GEAR_MODULE = "gear_module"
"""Chosen gear module [mm], an input so a mixed-integer master can set it."""

GEAR_TEETH = "gear_teeth"
"""Teeth on the small gear, likewise."""

_DEFAULT_MODULE = 1.5
_DEFAULT_TEETH = 25

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
        self._state: tuple[Any, ...] | None = None

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
        # Cached so ``_compute_jacobian`` can differentiate the very point that
        # was evaluated, without repeating the analysis.
        self._state = (
            design,
            analysis,
            dict(zip(MEMBER_NAMES, diameters, strict=True)),
            properties,
            loads,
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

    def _compute_jacobian(
        self,
        input_names: Iterable[str] = (),
        output_names: Iterable[str] = (),
    ) -> None:
        """Report the local Jacobian of the loads, exactly.

        Differentiating this discipline is what makes the coupled problem
        tractable.  Without it the whole MDA has to be differenced -- eleven
        converged fixed points per gradient, some fifty sweeps each.  With it
        GEMSEO assembles the coupled derivative from the two local Jacobians and
        one small linear solve.

        See :mod:`exlink.dynamics_jacobian` for the three ideas that carry it:
        the spectral operator is linear, the 18x18 solve is differentiated
        through its own factorisation, and the internal loads are closed form.
        """
        self._init_jacobian(input_names, output_names)
        if self._state is None:
            self.execute(dict(self.io.data))
        assert self._state is not None
        design, analysis, diameters, properties, loads = self._state

        size = len(MEMBER_NAMES) * self.samples * self.stations
        if not analysis.valid:
            return

        derivatives = coupled_jacobian(
            design,
            analysis,
            diameters,
            properties,
            loads,
            self.stations,
            self.material,
            self.spec,
        )
        flat_axial = derivatives.axial.reshape(size, N_PARAMETERS)
        flat_bending = derivatives.bending.reshape(size, N_PARAMETERS)

        rows = {
            COUPLING_AXIAL: flat_axial,
            COUPLING_BENDING: flat_bending,
            "dynamic_mean_torque": derivatives.mean_torque[None, :],
            "peak_bearing_load": derivatives.peak_bearing_load[None, :],
        }
        for output, block in rows.items():
            if output not in self.jac:
                continue
            for index, variable in enumerate(VARIABLE_NAMES):
                if variable in self.jac[output]:
                    self.jac[output][variable] = block[:, index : index + 1]
            if COUPLING_DIAMETERS in self.jac[output]:
                self.jac[output][COUPLING_DIAMETERS] = block[:, DIAMETER_SLICE]


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
        self._state: tuple[Any, ...] | None = None

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        data = dict(input_data)
        design = Design.from_mapping(data)
        shape = (len(MEMBER_NAMES), self.samples, self.stations)
        axial = np.asarray(data[COUPLING_AXIAL], dtype=float).reshape(shape)
        bending = np.asarray(data[COUPLING_BENDING], dtype=float).reshape(shape)
        lengths = member_lengths(design)

        sizing = size_from_arrays(axial, bending, lengths, self.material, self.safety)
        diameters = np.array([sizing[n].diameter for n in MEMBER_NAMES])
        self._state = (design, axial, bending, lengths, diameters)
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

    def _compute_jacobian(  # type: ignore[misc]
        self,
        input_names: Iterable[str] = (),
        output_names: Iterable[str] = (),
    ) -> None:
        """Report the local Jacobian of the sizing, exactly.

        The diameters are defined implicitly -- they are the sections that drive
        the worst utilisation to one -- so the implicit function theorem supplies
        their derivative and the bisection is never differentiated.  Everything
        downstream (member masses, the totals, the two margins) follows by chain
        rule from that.
        """
        self._init_jacobian(input_names, output_names)
        if self._state is None:
            self.execute(dict(self.io.data))
        assert self._state is not None
        design, axial, bending, lengths, diameters = self._state

        sizing = sizing_jacobian(axial, bending, diameters, lengths, self.material, self.safety)
        n_members = len(MEMBER_NAMES)
        size = n_members * self.samples * self.stations
        length_rows = member_length_jacobian(design)

        # d(diameter)/d(loads): block diagonal, one member per row.
        d_axial = np.zeros((n_members, size))
        d_bending = np.zeros((n_members, size))
        block = self.samples * self.stations
        for member in range(n_members):
            span = slice(member * block, (member + 1) * block)
            d_axial[member, span] = sizing.d_axial[member].ravel()
            d_bending[member, span] = sizing.d_bending[member].ravel()
        # d(diameter)/dX, through the member lengths only.
        d_design = sizing.d_length[:, None] * length_rows

        area = np.pi * diameters**2 / 4.0
        density = self.material.density
        # m_k = rho (pi d_k^2 / 4) L_k, and d_k itself depends on the loads.
        mass_from_diameter = density * np.pi * diameters / 2.0 * lengths
        d_mass_axial = mass_from_diameter[:, None] * d_axial
        d_mass_bending = mass_from_diameter[:, None] * d_bending
        d_mass_design = (
            mass_from_diameter[:, None] * d_design + density * area[:, None] * length_rows
        )

        slender = MEMBER_IS_SLENDER
        ratio = np.where(slender, diameters / lengths, -np.inf)
        critical = int(np.argmax(ratio))
        d_slender_axial = d_axial[critical] / lengths[critical]
        d_slender_bending = d_bending[critical] / lengths[critical]
        d_slender_design = (
            d_design[critical] / lengths[critical]
            - diameters[critical] * length_rows[critical] / lengths[critical] ** 2
        )
        widest = int(np.argmax(diameters))

        rows: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {
            COUPLING_DIAMETERS: (d_axial, d_bending, d_design),
            "member_mass": (d_mass_axial, d_mass_bending, d_mass_design),
            "structural_mass": (
                1000.0 * d_mass_axial.sum(axis=0)[None, :],
                1000.0 * d_mass_bending.sum(axis=0)[None, :],
                1000.0 * d_mass_design.sum(axis=0)[None, :],
            ),
            "total_mass": (
                1000.0 * d_mass_axial.sum(axis=0)[None, :],
                1000.0 * d_mass_bending.sum(axis=0)[None, :],
                1000.0 * d_mass_design.sum(axis=0)[None, :],
            ),
            "saturation_margin": (
                d_axial[widest][None, :],
                d_bending[widest][None, :],
                d_design[widest][None, :],
            ),
            "slenderness_margin": (
                d_slender_axial[None, :],
                d_slender_bending[None, :],
                d_slender_design[None, :],
            ),
        }
        for output, (by_axial, by_bending, by_design) in rows.items():
            if output not in self.jac:
                continue
            if COUPLING_AXIAL in self.jac[output]:
                self.jac[output][COUPLING_AXIAL] = by_axial
            if COUPLING_BENDING in self.jac[output]:
                self.jac[output][COUPLING_BENDING] = by_bending
            for index, variable in enumerate(VARIABLE_NAMES):
                if variable in self.jac[output]:
                    self.jac[output][variable] = by_design[:, index : index + 1]
        # The piston mass enters the total directly and nothing else.
        if "total_mass" in self.jac and "piston_mass" in self.jac["total_mass"]:
            self.jac["total_mass"]["piston_mass"] = np.ones((1, 1))


# ---------------------------------------------------------------------------------
# The vehicle-level objective
#
# Everything above answers "what does this mechanism do?".  This discipline
# answers "what is that worth?", which is the question the competition asks.
#
# It sits *downstream* of the MDA rather than inside it.  That is not a
# convenience: given converged section diameters, the load history is a direct
# function of the design, so the range chain is genuinely feed-forward and
# putting it inside the fixed point would add iterations that change nothing.
# It consumes ``diameters`` -- the MDA's converged coupling variable -- so every
# point the optimizer sees is a range computed on a self-consistent structure.
# ---------------------------------------------------------------------------------


class RangeDiscipline(Discipline):
    """Kilometres per litre, from the design and its converged sections.

    Chains friction, the engine mass budget and the burn-and-coast vehicle
    model.  The objective is exposed as ``neg_range`` so that a minimiser can
    take it directly, alongside the quantities a reader needs to interpret it.

    Derivatives are finite differences, and deliberately so.  The chain here
    contains a bisection on the active average-speed constraint, a maximum over
    the crank revolution, and a choice among discrete gear modules; the first
    two are differentiable almost everywhere but tedious to differentiate
    exactly, and the third is not differentiable at all.  What makes finite
    differences *safe* here -- and they are not safe on the geometric metrics,
    where :mod:`exlink.jacobian` exists precisely because they are not -- is
    that this discipline is cheap: one load solve, no fixed point, about 0.3 s.
    A full gradient costs eighteen of those, against a single MDA at 2.1 s.

    Args:
        speed_rpm: Crankshaft speed [rev/min].
        samples: Crank angles per revolution.
        vehicle: The car; a default Prototype-class entry if omitted.
        module: Gear module [mm]; the lightest workable standard one if omitted.
        teeth: Teeth on the small gear; derived from ``I`` if omitted.
        fluctuation: Allowed cyclic speed fluctuation, for the flywheel.
        material: The material.
        safety: The design factors.
        spec: Fixed engine data.
        step: Relative finite-difference step.
        name: Discipline name.
    """

    auto_detect_grammar_files: ClassVar[bool] = False

    def __init__(
        self,
        speed_rpm: float = DEFAULT_SPEED_RPM,
        samples: int = COUPLED_SAMPLES,
        vehicle: Vehicle | None = None,
        module: float | None = None,
        teeth: int | None = None,
        fluctuation: float = SPEED_FLUCTUATION,
        material: Material = DEFAULT_MATERIAL,
        safety: SafetyFactors = DEFAULT_SAFETY,
        spec: EngineSpec = DEFAULT_SPEC,
        step: float = 1.0e-6,
        name: str = "",
    ) -> None:
        super().__init__(name=name)
        self.speed_rpm = speed_rpm
        self.samples = samples
        self.vehicle = vehicle if vehicle is not None else Vehicle()
        self.module = module
        self.teeth = teeth
        self.fluctuation = fluctuation
        self.material = material
        self.safety = safety
        self.spec = spec
        self.step = step

        # The gear pair enters as *inputs* rather than as fixed construction
        # data, so that a mixed-integer master can choose it: under the Benders
        # formulation of :mod:`exlink.minlp` a catalogue interpolation maps the
        # master's one-hot selection onto these, and the sub-problem sees them
        # as constants.  Passing ``module``/``teeth`` to the constructor simply
        # sets their defaults, which is what every single-level scenario does.
        self.input_grammar.update_from_names(
            [*VARIABLE_NAMES, COUPLING_DIAMETERS, GEAR_MODULE, GEAR_TEETH]
        )
        self.output_grammar.update_from_names(list(RANGE_OUTPUTS))
        self.default_input_data = {
            **PUBLISHED_DESIGN.to_mapping(),
            COUPLING_DIAMETERS: _initial_diameters(),
            GEAR_MODULE: np.array([_DEFAULT_MODULE if module is None else module]),
            GEAR_TEETH: np.array([float(_DEFAULT_TEETH if teeth is None else teeth)]),
        }

    def _gear_pair(self, data: StrKeyMapping) -> tuple[float, int]:
        """The gear pair in force, from the inputs if given and the defaults otherwise."""
        module = data.get(GEAR_MODULE)
        teeth = data.get(GEAR_TEETH)
        chosen = self.module if module is None else float(np.ravel(module)[0])
        count = self.teeth if teeth is None else round(float(np.ravel(teeth)[0]))
        return (
            _DEFAULT_MODULE if chosen is None else chosen,
            _DEFAULT_TEETH if count is None else count,
        )

    def _evaluate(
        self,
        design: Design,
        diameters: np.ndarray,
        module: float | None = None,
        teeth: int | None = None,
    ) -> dict[str, float]:
        """The whole range chain, as plain floats.

        A design the kinematics cannot close returns the penalty row rather
        than raising, so that the optimizer can walk through infeasible ground
        the same way it does everywhere else in this package.
        """
        analysis = analyse(design, samples=self.samples, spec=self.spec)
        if not analysis.valid:
            return dict(RANGE_PENALTY)

        solved = analysis.require_solved()
        sections = {
            name: float(value)
            for name, value in zip(MEMBER_NAMES, np.asarray(diameters).ravel(), strict=True)
        }
        _crown, piston = piston_mass(
            solved.thermodynamics, self.material, self.safety, self.spec
        )
        properties = mass_properties(
            solved.kinematics, sections, self.material.density, piston, self.spec
        )
        loads = solve_dynamics(
            solved.kinematics,
            solved.thermodynamics.piston_force,
            properties,
            rpm_to_rad_per_s(self.speed_rpm),
            self.spec,
        )

        loss = friction_losses(loads, sections)
        metrics = analysis.metrics
        thermo = solved.thermodynamics
        quantity = heat_release(
            self.spec.dead_volume,
            thermo.p_compression_end,
            thermo.p_combustion,
            self.spec.heat_capacity_ratio,
        )
        budget = assemble(
            loads,
            sections,
            properties.member_mass,
            piston,
            metrics.height,
            metrics.width,
            metrics.expansion_stroke + metrics.compression_stroke,
            float(np.max(thermo.gauge_pressure)),
            module=self.module if module is None else module,
            teeth=self.teeth if teeth is None else teeth,
            fluctuation=self.fluctuation,
            spec=self.spec,
            material=self.material,
            safety=self.safety,
        )
        efficiency = brake_efficiency(loss.brake_work, quantity)
        power = loss.brake_work / 1000.0 * (self.speed_rpm / 60.0)
        outcome = best_strategy(self.vehicle, budget.total_kg, power, efficiency)

        pair = budget.gears
        return {
            "km_per_litre": outcome.km_per_litre,
            "neg_range": -outcome.km_per_litre,
            "engine_mass": budget.total_kg,
            "flywheel_mass": 1000.0 * budget.items.get("flywheel", 0.0),
            "crankcase_mass": 1000.0 * budget.items.get("crankcase", 0.0),
            "brake_efficiency": efficiency,
            "mechanical_efficiency": loss.mechanical_efficiency,
            "brake_power": power,
            # Constraint form: positive means the engine produces net work.
            "runs_margin": loss.brake_work / max(loss.indicated_work, 1.0),
            # Positive means the gear pair is within its face-width limit.
            "gear_margin": (MAX_WIDTH_FACTOR - pair.width_factor if pair is not None else -1.0),
            # The same two negated: positive means violated.  The bi-level
            # formulation needs this form; see RANGE_OUTPUTS.
            "runs_violation": -loss.brake_work / max(loss.indicated_work, 1.0),
            "gear_violation": (
                pair.width_factor - MAX_WIDTH_FACTOR if pair is not None else 1.0
            ),
        }

    def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
        data = dict(input_data)
        design = Design.from_mapping(data)
        diameters = np.asarray(data[COUPLING_DIAMETERS], dtype=float).ravel()
        module, teeth = self._gear_pair(data)
        values = self._evaluate(design, diameters, module, teeth)
        return {name: np.array([values[name]]) for name in RANGE_OUTPUTS}

    def _compute_jacobian(
        self,
        input_names: Sequence[str] = (),
        output_names: Sequence[str] = (),
    ) -> None:
        """Finite-difference Jacobian of the range chain.

        Central differences on a relative step.  Eighteen inputs, two chain
        evaluations each; the chain has no fixed point in it, so the whole
        gradient costs less than one MDA.
        """
        self._init_jacobian(input_names, output_names)
        data = dict(self.io.data)
        design = Design.from_mapping(data)
        base_vector = design.to_array()
        diameters = np.asarray(data[COUPLING_DIAMETERS], dtype=float).ravel()
        module, teeth = self._gear_pair(data)

        def perturbed(index: int, delta: float) -> dict[str, float]:
            if index < len(VARIABLE_NAMES):
                vector = base_vector.copy()
                vector[index] += delta
                return self._evaluate(Design.from_array(vector), diameters, module, teeth)
            if index < len(VARIABLE_NAMES) + diameters.size:
                sections = diameters.copy()
                sections[index - len(VARIABLE_NAMES)] += delta
                return self._evaluate(design, sections, module, teeth)
            # The module column: the gear pair is a master-level choice, so the
            # sub-problem holds it fixed and its derivative is not required.
            return self._evaluate(design, diameters, module, teeth)

        columns = [*base_vector, *diameters]
        gradients: dict[str, list[float]] = {name: [] for name in RANGE_OUTPUTS}
        for index, value in enumerate(columns):
            delta = self.step * max(abs(float(value)), 1.0)
            forward = perturbed(index, delta)
            backward = perturbed(index, -delta)
            for name in RANGE_OUTPUTS:
                gradients[name].append((forward[name] - backward[name]) / (2.0 * delta))

        for output in RANGE_OUTPUTS:
            if output not in self.jac:
                continue
            row = gradients[output]
            for index, variable in enumerate(VARIABLE_NAMES):
                if variable in self.jac[output]:
                    self.jac[output][variable] = np.array([[row[index]]])
            if COUPLING_DIAMETERS in self.jac[output]:
                self.jac[output][COUPLING_DIAMETERS] = np.array([row[len(VARIABLE_NAMES) :]])
