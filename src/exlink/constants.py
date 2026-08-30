"""Physical constants and specification data for the Honda EX-link derived engine.

Unit system used throughout the package
---------------------------------------
============  =====================================
quantity      unit
============  =====================================
length        mm
angle         rad (degrees only at the API surface)
area          mm^2
volume        mm^3
pressure      MPa (= N/mm^2); 1 bar = 0.1 MPa
force         N
torque        N.mm
============  =====================================

The numbers below are the ones handed to the students in the TN12 design
brief (Universite de Technologie de Compiegne, 2014-2015) and reproduced in
the 2015 report *Exlink Motor Mechanism Optimization*.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

BAR = 0.1
"""One bar expressed in MPa."""


@dataclass(frozen=True)
class EngineSpec:
    """Fixed data of the thermodynamic cycle and of the cylinder.

    Everything that is *not* a design variable lives here, so that a study can
    be repeated with a different brief by passing another :class:`EngineSpec`
    around instead of editing the code.
    """

    bore: float = 32.0
    """Piston diameter ``phi`` [mm]."""

    piston_length: float = 16.0
    """Distance ``p`` from the wrist pin ``P`` to the piston crown ``H`` [mm]."""

    piston_skirt: float = 16.0
    """Guided length ``h`` of the piston inside the liner [mm].

    Only enters the (unused) reaction moment ``M_D``; the report never states a
    value, so it defaults to the piston length.
    """

    dead_volume: float = 3000.0
    """Clearance volume ``V0`` at top dead centre [mm^3] (3 cc)."""

    p_intake: float = 1.2 * BAR
    """Intake / exhaust plenum pressure ``P0`` [MPa] (1.2 bar)."""

    heat_capacity_ratio: float = 1.22
    """Polytropic exponent ``gamma`` of compression and expansion."""

    explosion_ratio: float = 1.7
    """``k = P3 / P2``, the instantaneous pressure jump modelling combustion."""

    gear_ratio: float = 2.0
    """``r1 / r2`` between crankshaft and eccentric-shaft gears."""

    pressure_angle: float = math.radians(20.0)
    """Gear pressure angle ``alpha``.

    The report writes ``alpha`` symbolically but never gives a value; 20 deg is
    the standard involute profile and is used as the default.
    """

    @property
    def piston_area(self) -> float:
        """Cross-section ``pi phi^2 / 4`` on which the gas pressure acts [mm^2]."""
        return math.pi * self.bore**2 / 4.0


@dataclass(frozen=True)
class DesignTargets:
    """Right-hand sides of the equality and inequality constraints.

    These are the "table of design constraints" of the report plus the
    well-posedness constraints introduced in its optimization chapter.
    """

    expansion_stroke: float = 74.0
    """Required expansion stroke ``STE`` [mm] (equality constraint)."""

    compression_ratio: float = 16.0
    """Required compression ratio ``epsilon`` (equality constraint)."""

    max_rod_angle: float = 10.0
    """Upper bound on ``mra``, the piston-rod tilt away from vertical [deg]."""

    max_transmission: float = 0.985
    """Upper bound ``C`` on ``W = max(delta_c1, delta_c2)``.

    ``C = 0.9848`` in the report, i.e. the transmission angle ``T`` is kept
    inside ``[10 deg, 170 deg]``; it is rounded to 0.985 in the final
    formulation, which is the value used here.
    """

    max_tdc_gap: float = 0.01
    """Upper bound on ``g``, the gap between the two top dead centres [mm]."""

    min_clearance: float = 10.0
    """Lower bound on ``d``, trigonal-link to cylinder clearance [mm]."""

    max_side_load: float = 0.02
    """Upper bound on ``gamma = max(D) / max(P)``, the piston side-load ratio."""

    max_bearing_load: float = 25_000.0
    """Upper bound on the peak crankshaft bearing reaction [N].

    Only meaningful once inertia is in the load path: at rest this mechanism
    peaks near 7.7 kN, but the reaction grows as the square of engine speed and
    is what a plain bearing actually has to survive.
    """


@dataclass(frozen=True)
class PenaltyValues:
    """Values substituted when a design cannot be analysed at all.

    The report fixes ``eta(X) = 0``, ``H(X) = 1000``, ``B(X) = 1000`` whenever
    the kinematic compatibility conditions (4a)/(6a) fail, or whenever
    ``lambda(theta_1)`` does not have four monotone phases.
    """

    efficiency: float = 0.0
    height: float = 1000.0
    width: float = 1000.0


DEFAULT_SPEC = EngineSpec()
DEFAULT_TARGETS = DesignTargets()
DEFAULT_PENALTY = PenaltyValues()
