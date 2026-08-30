"""Fitting the approximate Atkinson cycle onto the piston motion.

The brief specifies an idealised cycle: intake and exhaust losses neglected,
compression and expansion adiabatic with a common exponent ``gamma``,
combustion modelled as an instantaneous pressure jump ``k = P3 / P2`` at top
dead centre, and blow-down modelled as an instantaneous fall to ``P0`` at the
end of expansion.

.. code-block:: text

    0 -> 1  intake       P = P0                       (TDC -> shallow BDC)
    1 -> 2  compression  P = P0 (V1 / V)^gamma        (shallow BDC -> TDC)
    2 -> 3  combustion   P3 = k P2                    (instantaneous, at TDC)
    3 -> 5  expansion    P = P3 (V0 / V)^gamma        (TDC -> deep BDC)
    5 -> 6  blow-down    P -> P0                      (instantaneous, at BDC)
    6 -> 0  exhaust      P = P0                       (deep BDC -> TDC)

What makes this an *Atkinson* cycle is that the piston reaches top dead centre
twice per crankshaft revolution while the two bottom dead centres differ: the
short one sets the compression stroke, the long one the expansion stroke.  So
before any pressure can be assigned, ``lambda(theta_1)`` must be shown to have
exactly four monotone phases -- two maxima and two distinct minima.  That test
is :func:`find_phases`, and it is the numerical heart of the report's
well-posedness argument: designs that fail it are penalised rather than
analysed.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec

FloatArray = NDArray[np.float64]


class Phase(IntEnum):
    """The four strokes, in the order they follow one another."""

    EXPANSION = 0
    EXHAUST = 1
    INTAKE = 2
    COMPRESSION = 3


class PhaseError(ValueError):
    """Raised when ``lambda(theta_1)`` is not a valid four-phase Atkinson motion."""


@dataclass(frozen=True)
class Phases:
    """Decomposition of one revolution into the four strokes."""

    labels: NDArray[np.int_]
    """Phase index of every crank angle, shaped ``(n_angles,)``."""

    lam_tdc: float
    """Piston height at top dead centre (the higher of the two maxima) [mm]."""

    tdc_gap: float
    """``g``, the gap between the two top dead centres [mm]."""

    expansion_stroke: float
    """``STE``, from top dead centre to the deeper bottom dead centre [mm]."""

    compression_stroke: float
    """``STC``, from top dead centre to the shallower bottom dead centre [mm]."""

    maxima_indices: tuple[int, int]
    minima_indices: tuple[int, int]


def _refine_extremum(lam: FloatArray, index: int) -> float:
    """Parabolic refinement of an extremum sampled on a uniform grid.

    The two top dead centres are compared to a tolerance of 0.01 mm while the
    default grid is 0.5 deg, so the sampled maximum is not accurate enough on
    its own; fitting a parabola through the three points around it recovers
    roughly three more digits.
    """
    n = lam.size
    y0, y1, y2 = lam[(index - 1) % n], lam[index], lam[(index + 1) % n]
    denominator = y0 - 2.0 * y1 + y2
    if denominator == 0.0:
        return float(y1)
    offset = 0.5 * (y0 - y2) / denominator
    if abs(offset) > 1.0:
        return float(y1)
    return float(y1 - 0.25 * (y0 - y2) * offset)


def find_phases(lam: FloatArray) -> Phases:
    """Split one revolution of ``lambda(theta_1)`` into the four strokes.

    Args:
        lam: Piston height at uniformly spaced crank angles over ``[0, 2 pi)``.

    Returns:
        The phase decomposition and the two strokes.

    Raises:
        PhaseError: If the motion does not consist of exactly four monotone
            phases with two maxima and two minima -- i.e. the design does not
            realise an Atkinson cycle at all (a plain Otto motion with a single
            up-and-down per revolution fails here, as do the "stair-stepped"
            motions shown in the report).
    """
    n = lam.size
    slope = np.sign(np.diff(np.concatenate([lam, lam[:1]])))
    if np.any(slope == 0.0):
        # A plateau makes the phase count ill-defined; nudge it to the previous
        # sign so that a flat sample does not read as an extra phase.
        for index in np.nonzero(slope == 0.0)[0]:
            slope[index] = slope[index - 1]
    turning = np.nonzero(slope != np.roll(slope, 1))[0]
    if turning.size != 4:
        msg = f"lambda(theta_1) has {turning.size} monotone phases, expected 4"
        raise PhaseError(msg)

    maxima = [int(i) for i in turning if slope[i] < 0]
    minima = [int(i) for i in turning if slope[i] > 0]
    if len(maxima) != 2 or len(minima) != 2:
        msg = "lambda(theta_1) must have exactly two maxima and two minima"
        raise PhaseError(msg)

    max_values = [_refine_extremum(lam, i) for i in maxima]
    min_values = [_refine_extremum(lam, i) for i in minima]
    lam_tdc = max(max_values)
    tdc_gap = abs(max_values[0] - max_values[1])

    deep = minima[int(np.argmin(min_values))]
    shallow = minima[int(np.argmax(min_values))]
    expansion_stroke = lam_tdc - min(min_values)
    compression_stroke = lam_tdc - max(min_values)
    # Two *equal* minima are not rejected here.  Such a motion is four-phase but
    # symmetric, so STE == STC and the cycle is Otto, not Atkinson -- and that is
    # caught downstream by the equality constraints, which demand STE = 74 mm and
    # epsilon = 16 simultaneously.  Rejecting it here would need an arbitrary
    # threshold and would blind the optimizer to a region it may legitimately
    # pass through.

    # Walk the revolution and label each sample by the stroke it belongs to.
    # Expansion runs from a top dead centre down to the deeper bottom dead
    # centre; the rest follows in order.
    labels = np.empty(n, dtype=int)
    order = sorted(turning)
    for start, end in zip(order, order[1:] + order[:1], strict=True):
        indices = np.arange(start + 1, end + 1 if end > start else end + 1 + n) % n
        if end == deep:
            phase = Phase.EXPANSION
        elif start == deep:
            phase = Phase.EXHAUST
        elif end == shallow:
            phase = Phase.INTAKE
        else:
            phase = Phase.COMPRESSION
        labels[indices] = int(phase)

    return Phases(
        labels=labels,
        lam_tdc=lam_tdc,
        tdc_gap=tdc_gap,
        expansion_stroke=expansion_stroke,
        compression_stroke=compression_stroke,
        maxima_indices=(maxima[0], maxima[1]),
        minima_indices=(deep, shallow),
    )


@dataclass(frozen=True)
class Thermodynamics:
    """In-cylinder state over one crankshaft revolution."""

    phases: Phases
    volume: FloatArray
    """Instantaneous cylinder volume ``V(theta_1)`` [mm^3]."""

    pressure: FloatArray
    """Absolute in-cylinder pressure ``P(theta_1)`` [MPa]."""

    gauge_pressure: FloatArray
    """``P - P0``, the pressure that actually loads the piston [MPa]."""

    piston_force: FloatArray
    """Gas force on the piston crown, ``(P - P0) pi phi^2 / 4`` [N]."""

    compression_ratio: float
    """``epsilon = V1 / V0`` realised by this design."""

    p_compression_end: float
    """``P2`` at the end of compression [MPa]."""

    p_combustion: float
    """``P3 = k P2`` just after combustion [MPa]."""


def solve(
    lam: FloatArray,
    spec: EngineSpec = DEFAULT_SPEC,
    phases: Phases | None = None,
) -> Thermodynamics:
    """Apply the approximate Atkinson cycle to a piston motion.

    Args:
        lam: Piston height at uniformly spaced crank angles over ``[0, 2 pi)``.
        spec: Fixed engine data.
        phases: A precomputed decomposition; recomputed from ``lam`` if omitted.

    Returns:
        Volume, pressure and piston force at every crank angle.

    Raises:
        PhaseError: Propagated from :func:`find_phases`.
    """
    phases = find_phases(lam) if phases is None else phases

    area = spec.piston_area
    volume = spec.dead_volume + area * (phases.lam_tdc - lam)

    v1 = spec.dead_volume + area * phases.compression_stroke
    compression_ratio = v1 / spec.dead_volume
    p2 = spec.p_intake * compression_ratio**spec.heat_capacity_ratio
    p3 = spec.explosion_ratio * p2

    pressure = np.full_like(lam, spec.p_intake)
    is_compression = phases.labels == int(Phase.COMPRESSION)
    is_expansion = phases.labels == int(Phase.EXPANSION)
    gamma = spec.heat_capacity_ratio
    pressure[is_compression] = spec.p_intake * (v1 / volume[is_compression]) ** gamma
    pressure[is_expansion] = p3 * (spec.dead_volume / volume[is_expansion]) ** gamma

    gauge_pressure = pressure - spec.p_intake
    return Thermodynamics(
        phases=phases,
        volume=volume,
        pressure=pressure,
        gauge_pressure=gauge_pressure,
        piston_force=gauge_pressure * area,
        compression_ratio=compression_ratio,
        p_compression_end=p2,
        p_combustion=p3,
    )
