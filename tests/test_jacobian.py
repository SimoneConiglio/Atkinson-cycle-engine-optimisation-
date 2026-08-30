"""Exact derivatives of the analysis chain.

Every check here compares against a *converged* central difference. The step
matters: several of these metrics are maxima over the crank angle, and at a step
large enough to move the sample attaining the maximum, the difference quotient
straddles the switch and is simply wrong. That failure mode is itself pinned
below, because it is the reason these derivatives are computed analytically.
"""

from __future__ import annotations

import numpy as np
import pytest

from exlink import Design, analyse
from exlink.design import VARIABLE_NAMES
from exlink.disciplines import ExlinkDiscipline
from exlink.jacobian import (
    ANALYTIC_OUTPUTS,
    kinematic_jacobian,
    metric_jacobian,
)
from exlink.kinematics import solve
from exlink.reference import REFINED_DESIGN

SAMPLES = 720
STEP = 1.0e-6


@pytest.fixture(scope="module")
def analysis():
    return analyse(REFINED_DESIGN, samples=SAMPLES)


@pytest.fixture(scope="module")
def kinematic(analysis):
    return kinematic_jacobian(REFINED_DESIGN, analysis.kinematics)


@pytest.fixture(scope="module")
def metrics(analysis, kinematic):
    return metric_jacobian(REFINED_DESIGN, analysis, kinematic)


def _difference(index: int, extract, step: float = STEP) -> float:
    """Central difference of a scalar metric with respect to variable ``index``."""
    base = REFINED_DESIGN.to_array()
    forward, backward = base.copy(), base.copy()
    forward[index] += step
    backward[index] -= step
    return (
        extract(analyse(Design.from_array(forward), samples=SAMPLES).metrics)
        - extract(analyse(Design.from_array(backward), samples=SAMPLES).metrics)
    ) / (2.0 * step)


@pytest.mark.parametrize("field", ["lam", "theta_e", "theta_T", "theta_a"])
def test_kinematic_derivatives_match_finite_differences(kinematic, field) -> None:
    """The forward-mode chain reproduces the whole history, angle by angle."""
    base = REFINED_DESIGN.to_array()
    for index in range(len(VARIABLE_NAMES)):
        step = STEP * max(abs(base[index]), 1.0)
        forward, backward = base.copy(), base.copy()
        forward[index] += step
        backward[index] -= step
        numerical = (
            getattr(solve(Design.from_array(forward), samples=SAMPLES), field)
            - getattr(solve(Design.from_array(backward), samples=SAMPLES), field)
        ) / (2.0 * step)
        analytic = getattr(kinematic, field)[:, index]
        scale = max(float(np.max(np.abs(numerical))), 1e-9)
        assert np.max(np.abs(analytic - numerical)) < 1e-6 * scale, VARIABLE_NAMES[index]


METRIC_EXTRACTORS = {
    "expansion_stroke": (lambda m: m.expansion_stroke, 74.0),
    "compression_ratio": (lambda m: m.compression_ratio, 16.0),
    "compatibility": (lambda m: m.compatibility, 1.0),
    "rod_angle": (lambda m: m.rod_angle, 10.0),
    "side_load_ratio": (lambda m: m.side_load_ratio, 0.02),
    "tdc_gap_margin": (lambda m: m.tdc_gap, 0.01),
}


@pytest.mark.parametrize("name", list(METRIC_EXTRACTORS))
def test_metric_derivatives_match_finite_differences(metrics, name) -> None:
    """Each tight metric, against a converged difference.

    Compared on the scale of the metric's own constraint bound rather than
    relatively: ``g`` is a difference of two nearly equal maxima, so its
    relative error is dominated by cancellation while its absolute error stays
    far below anything the optimizer can see.
    """
    extract, scale = METRIC_EXTRACTORS[name]
    for index in range(len(VARIABLE_NAMES)):
        numerical = _difference(index, extract)
        assert metrics[name][index] == pytest.approx(numerical, abs=1e-4 * scale), (
            VARIABLE_NAMES[index]
        )


def test_every_declared_analytic_output_is_produced(metrics) -> None:
    assert set(ANALYTIC_OUTPUTS) == set(metrics)
    for value in metrics.values():
        assert value.shape == (len(VARIABLE_NAMES),)
        assert np.all(np.isfinite(value))


def test_a_coarse_step_makes_the_difference_quotient_wrong(metrics) -> None:
    """The reason these are analytic rather than differenced.

    ``gamma`` is the peak side load over the peak gas load. At a step of 1e-4 mm
    on the swing rod the crank angle attaining the peak side load moves by one
    sample, and the central difference -- taken across that switch -- lands about
    20 % away from the true derivative. The envelope theorem has no such
    failure: it evaluates the partial derivative *at* the maximiser.
    """
    index = VARIABLE_NAMES.index("a")
    extract, _ = METRIC_EXTRACTORS["side_load_ratio"]
    exact = metrics["side_load_ratio"][index]

    converged = _difference(index, extract, step=1e-6)
    coarse = _difference(index, extract, step=1e-4)

    assert exact == pytest.approx(converged, rel=1e-3)
    assert abs(coarse - exact) > 0.1 * abs(exact)


def test_gemseo_accepts_the_discipline_jacobian() -> None:
    """GEMSEO's own checker, on the outputs claimed analytic."""
    discipline = ExlinkDiscipline(samples=SAMPLES)
    assert discipline.check_jacobian(
        REFINED_DESIGN.to_mapping(),
        threshold=1e-4,
        step=1e-6,
        output_names=[
            "compatibility_margin",
            "stroke_error",
            "compression_ratio_error",
            "rod_angle_margin",
            "side_load_margin",
            "tdc_gap_margin",
        ],
    )


def test_the_discipline_fills_every_jacobian_entry() -> None:
    """Outputs left to differences still get a gradient, not a hole."""
    discipline = ExlinkDiscipline(samples=360)
    discipline.linearize(REFINED_DESIGN.to_mapping(), compute_all_jacobians=True)
    for output in ("neg_efficiency", "height", "width", "clearance"):
        row = np.array([float(discipline.jac[output][name][0, 0]) for name in VARIABLE_NAMES])
        assert np.all(np.isfinite(row))
        assert np.any(row != 0.0)
