"""End-to-end analysis, penalisation, and the constraint vectors."""

from __future__ import annotations

import numpy as np
import pytest

from exlink import DEFAULT_TARGETS, PUBLISHED_DESIGN, analyse
from exlink.metrics import cylinder_clearance, envelope, rod_angle_deviation
from exlink.model import (
    EQUALITY_NAMES,
    INEQUALITY_NAMES,
    equality_constraints,
    inequality_constraints,
    objectives,
)
from exlink.reference import REFINED_DESIGN


def test_the_published_design_analyses(published_analysis) -> None:
    assert published_analysis.valid
    assert published_analysis.metrics.efficiency > 0.0
    assert published_analysis.metrics.mean_torque > 0.0


def test_an_incompatible_design_is_penalised_not_raised() -> None:
    """The report's rule: ``eta = 0``, ``H = B = 1000``, never an exception."""
    analysis = analyse(PUBLISHED_DESIGN.replace(a=25.0, c=25.0), samples=180)
    assert not analysis.valid
    assert analysis.metrics.efficiency == 0.0
    assert analysis.metrics.height == 1000.0
    assert analysis.metrics.width == 1000.0
    assert "incompatible" in analysis.metrics.reason


def test_an_otto_motion_is_penalised() -> None:
    """A linkage that closes fine but lifts the piston only once per revolution.

    Shrinking the eccentric crank collapses the second top dead centre, leaving
    a plain Otto motion -- the exact failure the report describes. The linkage
    is still assemblable, so only the phase count catches it.
    """
    analysis = analyse(REFINED_DESIGN.replace(q_2=2.0), samples=360)
    assert not analysis.valid
    assert "2 monotone phases" in analysis.metrics.reason
    assert analysis.metrics.efficiency == 0.0
    assert analysis.metrics.height == 1000.0


def test_a_six_phase_motion_is_penalised() -> None:
    """The other failure mode: an over-long crank adds a spurious pair of phases."""
    analysis = analyse(REFINED_DESIGN.replace(q_1=40.0), samples=360)
    assert not analysis.valid
    assert "6 monotone phases" in analysis.metrics.reason


def test_penalised_designs_still_report_compatibility() -> None:
    """The optimizer needs a gradient-free signal out of a failed design."""
    analysis = analyse(PUBLISHED_DESIGN.replace(a=25.0, c=25.0), samples=180)
    assert np.isfinite(analysis.metrics.compatibility)
    assert analysis.metrics.compatibility >= 1.0


def test_objectives_are_negated_efficiency_and_the_two_sizes(published_analysis) -> None:
    f = objectives(published_analysis)
    assert f.shape == (3,)
    assert f[0] == pytest.approx(-published_analysis.metrics.efficiency)
    assert f[1] == pytest.approx(published_analysis.metrics.height)
    assert f[2] == pytest.approx(published_analysis.metrics.width)


def test_constraint_vectors_have_the_declared_lengths(published_analysis) -> None:
    assert inequality_constraints(published_analysis).shape == (len(INEQUALITY_NAMES),)
    assert equality_constraints(published_analysis).shape == (len(EQUALITY_NAMES),)


def test_constraints_are_residuals_against_the_targets(refined_analysis) -> None:
    metrics = refined_analysis.metrics
    ineq = inequality_constraints(refined_analysis)
    eq = equality_constraints(refined_analysis)
    assert ineq[0] == pytest.approx(metrics.rod_angle - DEFAULT_TARGETS.max_rod_angle)
    assert ineq[1] == pytest.approx(metrics.compatibility - DEFAULT_TARGETS.max_transmission)
    assert ineq[3] == pytest.approx(DEFAULT_TARGETS.min_clearance - metrics.clearance)
    assert eq[0] == pytest.approx(metrics.expansion_stroke - DEFAULT_TARGETS.expansion_stroke)
    assert eq[1] == pytest.approx(metrics.compression_ratio - DEFAULT_TARGETS.compression_ratio)


def test_the_refined_reference_is_feasible(refined_analysis) -> None:
    """The design shipped as the reference must satisfy every constraint."""
    from exlink.scenarios import is_feasible

    assert is_feasible(refined_analysis)


def test_envelope_contains_every_joint(published_analysis) -> None:
    height, width = envelope(published_analysis.kinematics)
    kinematics = published_analysis.kinematics
    cloud = np.concatenate([kinematics.Q, kinematics.A, kinematics.D, kinematics.E])
    assert width >= np.ptp(cloud[:, 0])
    assert height >= np.ptp(cloud[:, 1])


def test_envelope_contains_the_gear_primitives(published_analysis) -> None:
    height, width = envelope(published_analysis.kinematics)
    assert width >= 2.0 * published_analysis.design.r_1
    assert height >= 2.0 * published_analysis.design.r_1


def test_rod_angle_is_the_deviation_from_vertical(published_analysis) -> None:
    theta_e = published_analysis.kinematics.theta_e
    expected = np.degrees(np.max(np.abs(theta_e - np.pi / 2.0)))
    assert rod_angle_deviation(published_analysis.kinematics) == pytest.approx(expected)


def test_clearance_is_non_negative(published_analysis) -> None:
    assert cylinder_clearance(published_analysis.kinematics) >= 0.0


def test_clearance_falls_when_the_cylinder_moves_over_the_linkage() -> None:
    """Sliding the cylinder axis toward the trigonal link must reduce ``d``."""
    far = analyse(PUBLISHED_DESIGN, samples=360)
    near = analyse(PUBLISHED_DESIGN.replace(x_1=40.0), samples=360)
    if not near.valid:
        pytest.skip("the shifted design is not analysable")
    assert near.metrics.clearance < far.metrics.clearance


def test_analysis_is_reproducible(published_analysis) -> None:
    again = analyse(PUBLISHED_DESIGN, samples=1440)
    assert again.metrics.efficiency == pytest.approx(published_analysis.metrics.efficiency)


def test_metrics_converge_with_resolution() -> None:
    """Doubling the crank-angle grid must not move the answers much."""
    coarse = analyse(PUBLISHED_DESIGN, samples=720).metrics
    fine = analyse(PUBLISHED_DESIGN, samples=2880).metrics
    assert fine.efficiency == pytest.approx(coarse.efficiency, rel=5e-3)
    assert fine.expansion_stroke == pytest.approx(coarse.expansion_stroke, abs=0.05)
