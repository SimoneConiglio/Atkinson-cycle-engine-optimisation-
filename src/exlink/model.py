"""One-call analysis of an EX-link design.

:func:`analyse` chains kinematics -> Atkinson cycle -> quasi-static loads ->
metrics, and -- crucially for the optimizer -- never raises on a bad design.
Following the report, a design that cannot be analysed is *penalised* rather
than rejected:

* the compatibility conditions (4a)/(6a) fail, so the crankshaft only rocks; or
* ``lambda(theta_1)`` does not have four monotone phases, so no Atkinson cycle
  can be fitted onto it.

In either case the report substitutes ``eta = 0``, ``H = 1000``, ``B = 1000``.
The constraint measures that *can* still be computed (notably ``W``) are
returned unchanged, so a gradient-free optimizer still has a signal telling it
which way to go -- that signal is what makes the global search converge.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import cycle
from . import kinematics as kinematics_module
from . import loads as loads_module
from . import metrics as metrics_module
from .constants import (
    DEFAULT_PENALTY,
    DEFAULT_SPEC,
    DEFAULT_TARGETS,
    DesignTargets,
    EngineSpec,
    PenaltyValues,
)
from .design import Design
from .kinematics import DEFAULT_SAMPLES, Kinematics
from .metrics import Metrics


@dataclass(frozen=True)
class SolvedAnalysis:
    """An :class:`Analysis` that succeeded, with none of its members optional.

    The plotting and animation code needs the kinematics, the cycle and the
    loads, all of which are ``None`` on a penalised design.  Narrowing once, at
    the entry point, keeps every downstream helper free of ``None`` checks and
    lets the type checker verify that.
    """

    design: Design
    spec: EngineSpec
    metrics: Metrics
    kinematics: Kinematics
    thermodynamics: cycle.Thermodynamics
    loads: loads_module.Loads


@dataclass(frozen=True)
class Analysis:
    """Everything computed about one design."""

    design: Design
    spec: EngineSpec
    metrics: Metrics
    kinematics: Kinematics | None
    thermodynamics: cycle.Thermodynamics | None
    loads: loads_module.Loads | None

    @property
    def valid(self) -> bool:
        """Whether the design could be analysed at all."""
        return self.metrics.valid

    def require_solved(self) -> SolvedAnalysis:
        """Return a non-optional view of this analysis.

        Returns:
            The same analysis, with every member guaranteed present.

        Raises:
            ValueError: If the design was penalised rather than analysed.
        """
        if (
            not self.valid
            or self.kinematics is None
            or self.thermodynamics is None
            or self.loads is None
        ):
            msg = f"design was not analysable: {self.metrics.reason}"
            raise ValueError(msg)
        return SolvedAnalysis(
            design=self.design,
            spec=self.spec,
            metrics=self.metrics,
            kinematics=self.kinematics,
            thermodynamics=self.thermodynamics,
            loads=self.loads,
        )


def _penalised(
    design: Design,
    spec: EngineSpec,
    penalty: PenaltyValues,
    reason: str,
    kinematics: Kinematics | None,
) -> Analysis:
    """Build the penalised outcome for a design that cannot be analysed."""
    compatibility = kinematics.compatibility if kinematics is not None else float("inf")
    rod_angle = (
        metrics_module.rod_angle_deviation(kinematics)
        if kinematics is not None and kinematics.feasible
        else 180.0
    )
    return Analysis(
        design=design,
        spec=spec,
        metrics=Metrics(
            efficiency=penalty.efficiency,
            torque_pressure_ratio=0.0,
            lever_arm=0.0,
            height=penalty.height,
            width=penalty.width,
            expansion_stroke=0.0,
            compression_stroke=0.0,
            compression_ratio=1.0,
            rod_angle=rod_angle,
            compatibility=compatibility,
            tdc_gap=penalty.height,
            clearance=0.0,
            side_load_ratio=float("inf"),
            mean_torque=0.0,
            mean_piston_force=0.0,
            valid=False,
            reason=reason,
        ),
        kinematics=kinematics,
        thermodynamics=None,
        loads=None,
    )


def analyse(
    design: Design,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    penalty: PenaltyValues = DEFAULT_PENALTY,
) -> Analysis:
    """Analyse a design end to end.

    Args:
        design: The mechanism dimensions.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides (used only to gate the analysis).
        penalty: Values substituted for an unanalysable design.

    Returns:
        The full analysis; check :attr:`Analysis.valid` before reading the
        physical quantities.
    """
    kin = kinematics_module.solve(design, samples=samples, spec=spec)

    # Gate 1: can the crankshaft turn through a full revolution at all?
    if not kin.feasible:
        return _penalised(design, spec, penalty, "kinematically incompatible (W >= 1)", kin)

    # Gate 2: does the piston motion have four monotone phases with two equal
    # top dead centres and two distinct bottom dead centres?
    try:
        thermo = cycle.solve(kin.lam, spec=spec)
    except cycle.PhaseError as error:
        return _penalised(design, spec, penalty, str(error), kin)

    loads = loads_module.solve(kin, thermo.piston_force, spec=spec)
    phases = thermo.phases
    eta, phi, lever = metrics_module.efficiency(
        loads, phases.expansion_stroke, phases.compression_stroke, spec=spec
    )
    height, width = metrics_module.envelope(kin, spec=spec)

    metrics = Metrics(
        efficiency=eta,
        torque_pressure_ratio=phi,
        lever_arm=lever,
        height=height,
        width=width,
        expansion_stroke=phases.expansion_stroke,
        compression_stroke=phases.compression_stroke,
        compression_ratio=thermo.compression_ratio,
        rod_angle=metrics_module.rod_angle_deviation(kin),
        compatibility=kin.compatibility,
        tdc_gap=phases.tdc_gap,
        clearance=metrics_module.cylinder_clearance(kin, spec=spec),
        side_load_ratio=loads.side_load_ratio,
        mean_torque=loads.mean_torque,
        mean_piston_force=loads.mean_piston_force,
        valid=True,
    )
    return Analysis(
        design=design,
        spec=spec,
        metrics=metrics,
        kinematics=kin,
        thermodynamics=thermo,
        loads=loads,
    )


def objectives(analysis: Analysis) -> np.ndarray:
    """``f(X) = (-eta, H, B)^T`` of the report's final formulation."""
    m = analysis.metrics
    return np.array([-m.efficiency, m.height, m.width], dtype=float)


def inequality_constraints(
    analysis: Analysis, targets: DesignTargets = DEFAULT_TARGETS
) -> np.ndarray:
    """``c(X) <= 0`` of the report's final formulation.

    Order: ``mra - 10``, ``W - 0.985``, ``g - 0.01``, ``10 - d``,
    ``gamma - 0.02``.
    """
    m = analysis.metrics
    side_load = m.side_load_ratio
    if not np.isfinite(side_load):
        side_load = 1.0e3
    return np.array(
        [
            m.rod_angle - targets.max_rod_angle,
            m.compatibility - targets.max_transmission,
            m.tdc_gap - targets.max_tdc_gap,
            targets.min_clearance - m.clearance,
            side_load - targets.max_side_load,
        ],
        dtype=float,
    )


def equality_constraints(
    analysis: Analysis, targets: DesignTargets = DEFAULT_TARGETS
) -> np.ndarray:
    """``c_eq(X) = 0`` of the report's final formulation.

    Order: ``STE - 74``, ``epsilon - 16``.
    """
    m = analysis.metrics
    return np.array(
        [
            m.expansion_stroke - targets.expansion_stroke,
            m.compression_ratio - targets.compression_ratio,
        ],
        dtype=float,
    )


#: Names of the inequality constraints, in the order returned above.
INEQUALITY_NAMES: tuple[str, ...] = (
    "rod_angle_margin",
    "compatibility_margin",
    "tdc_gap_margin",
    "clearance_margin",
    "side_load_margin",
)

#: Names of the equality constraints, in the order returned above.
EQUALITY_NAMES: tuple[str, ...] = ("stroke_error", "compression_ratio_error")
