"""Phase detection and the approximate Atkinson cycle."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.cycle import Phase, PhaseError, find_phases


def test_four_phases_are_found(published_analysis) -> None:
    phases = published_analysis.thermodynamics.phases
    assert set(np.unique(phases.labels)) == {int(p) for p in Phase}


def test_strokes_are_measured_from_top_dead_centre(published_analysis) -> None:
    kinematics = published_analysis.kinematics
    phases = published_analysis.thermodynamics.phases
    assert phases.lam_tdc == pytest.approx(kinematics.lam.max(), abs=1e-3)
    assert phases.expansion_stroke == pytest.approx(
        phases.lam_tdc - kinematics.lam.min(), abs=1e-3
    )
    assert phases.expansion_stroke > phases.compression_stroke


def test_a_single_up_and_down_is_rejected() -> None:
    """A plain Otto motion has two phases, not four, and must be rejected."""
    angles = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False)
    otto = 100.0 + 37.0 * np.cos(angles)
    with pytest.raises(PhaseError, match="2 monotone phases"):
        find_phases(otto)


def test_a_six_phase_motion_is_rejected() -> None:
    angles = np.linspace(0.0, 2.0 * np.pi, 720, endpoint=False)
    with pytest.raises(PhaseError, match="6 monotone phases"):
        find_phases(100.0 + 37.0 * np.cos(3.0 * angles))


def test_a_symmetric_motion_gives_equal_strokes() -> None:
    """Two identical minima are four-phase, but Otto rather than Atkinson.

    ``find_phases`` accepts it -- rejecting it would need an arbitrary threshold
    on "different enough" minima. What rules it out is the pair of equality
    constraints downstream, which cannot both hold when STE == STC: a 74 mm
    stroke on both sides would give a compression ratio far from 16.
    """
    angles = np.linspace(0.0, 2.0 * np.pi, 720, endpoint=False)
    phases = find_phases(100.0 + 37.0 * np.cos(2.0 * angles))
    assert phases.expansion_stroke == pytest.approx(phases.compression_stroke)
    assert phases.tdc_gap == pytest.approx(0.0, abs=1e-9)


def test_the_tdc_gap_is_measured() -> None:
    """A motion with deliberately unequal maxima must report the gap as ``g``."""
    angles = np.linspace(0.0, 2.0 * np.pi, 2880, endpoint=False)
    lam = 100.0 + 30.0 * np.cos(2.0 * angles) + 5.0 * np.cos(angles)
    phases = find_phases(lam)
    assert phases.tdc_gap == pytest.approx(10.0, abs=0.05)


def test_parabolic_refinement_beats_the_sample_grid() -> None:
    """The gap tolerance is 0.01 mm on a 0.5 deg grid, so refinement matters."""
    angles = np.linspace(0.0, 2.0 * np.pi, 720, endpoint=False)
    lam = 100.0 + 30.0 * np.cos(2.0 * angles) + 0.5 * np.cos(angles)
    phases = find_phases(lam)
    assert phases.tdc_gap == pytest.approx(1.0, abs=1e-3)


def test_pressure_follows_the_specified_cycle(published_analysis) -> None:
    thermo = published_analysis.thermodynamics
    spec = published_analysis.spec
    labels = thermo.phases.labels

    # Intake and exhaust sit at the plenum pressure, so the piston is unloaded.
    breathing = np.isin(labels, [int(Phase.INTAKE), int(Phase.EXHAUST)])
    assert thermo.pressure[breathing] == pytest.approx(spec.p_intake)
    assert thermo.gauge_pressure[breathing] == pytest.approx(0.0)

    # Compression and expansion are above it.
    working = ~breathing
    assert np.all(thermo.pressure[working] >= spec.p_intake - 1e-12)


def test_the_combustion_jump_uses_the_explosion_ratio(published_analysis) -> None:
    thermo = published_analysis.thermodynamics
    spec = published_analysis.spec
    assert thermo.p_combustion == pytest.approx(spec.explosion_ratio * thermo.p_compression_end)
    assert thermo.p_compression_end == pytest.approx(
        spec.p_intake * thermo.compression_ratio**spec.heat_capacity_ratio
    )


def test_compression_ratio_follows_from_the_compression_stroke(published_analysis) -> None:
    thermo = published_analysis.thermodynamics
    spec = published_analysis.spec
    expected = (
        spec.dead_volume + spec.piston_area * thermo.phases.compression_stroke
    ) / spec.dead_volume
    assert thermo.compression_ratio == pytest.approx(expected)


def test_volume_never_falls_below_the_dead_volume(published_analysis) -> None:
    """``V >= V0`` everywhere, and it touches ``V0`` at top dead centre.

    ``lam_tdc`` is the parabolically refined maximum, which sits a hair above
    the highest sampled point, so the minimum volume lands just above ``V0``
    rather than exactly on it.
    """
    thermo = published_analysis.thermodynamics
    dead_volume = published_analysis.spec.dead_volume
    assert thermo.volume.min() >= dead_volume
    assert thermo.volume.min() == pytest.approx(dead_volume, abs=1.0)


def test_expansion_is_longer_than_compression(published_analysis) -> None:
    """The defining property of the Atkinson cycle."""
    thermo = published_analysis.thermodynamics
    expansion = thermo.volume[thermo.phases.labels == int(Phase.EXPANSION)]
    compression = thermo.volume[thermo.phases.labels == int(Phase.COMPRESSION)]
    assert expansion.max() > compression.max()


def test_a_design_meeting_the_target_reaches_the_target_ratio(refined_analysis) -> None:
    assert refined_analysis.metrics.compression_ratio == pytest.approx(16.0, abs=0.05)
    assert refined_analysis.metrics.expansion_stroke == pytest.approx(74.0, abs=0.05)
