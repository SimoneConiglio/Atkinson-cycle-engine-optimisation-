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

The numbers below define a single-cylinder Shell Eco-marathon engine and are
held fixed throughout; only the linkage geometry and the member sections are
designed.
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

    Only enters the (unused) reaction moment ``M_D``, so it simply defaults to
    the piston length.
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
    """``r1 / r2`` between the analysis shaft's gear and the output shaft's.

    The kinematics of :mod:`exlink.kinematics` are parametrised on ``theta_1``,
    the shaft carrying ``q_1`` and the *large* gear, and the four strokes
    complete in one turn of it.  The second shaft therefore turns
    ``theta_2 = -2 theta_1 + theta_f``: twice per cycle, or **720 deg**, which is
    what a conventional four-stroke crankshaft does.  Power is taken from that
    shaft, so ``theta_1`` is the half-speed member of the pair -- Honda's own
    arrangement -- and every engine speed the study quotes is the output
    shaft's, at :attr:`output_revolutions_per_cycle` times the ``speed_rpm``
    the analysis functions take.
    """

    pressure_angle: float = math.radians(20.0)
    """Gear pressure angle ``alpha``.

    20 deg is the standard involute profile and is used throughout; ``alpha``
    enters only the gear-mesh separating force.
    """

    @property
    def piston_area(self) -> float:
        """Cross-section ``pi phi^2 / 4`` on which the gas pressure acts [mm^2]."""
        return math.pi * self.bore**2 / 4.0

    @property
    def output_revolutions_per_cycle(self) -> float:
        """Turns of the output shaft per four-stroke cycle [-].

        Two, from :attr:`gear_ratio`, and the same as a conventional engine --
        which is what makes the comparison of section 6.3 a comparison at equal
        speed and equal firing rate rather than one needing a correction.
        """
        return self.gear_ratio

    def output_speed_rpm(self, speed_rpm: float) -> float:
        """Output-shaft speed for an analysis speed ``theta_1`` [rev/min]."""
        return float(speed_rpm) * self.gear_ratio

    def output_torque(self, torque: float) -> float:
        """Torque referred from ``theta_1`` to the output shaft [N.mm].

        The virtual-work identity of :mod:`exlink.loads` makes ``M_r`` the
        *whole* engine torque referred to ``theta_1``, so referring it to a
        shaft turning ``gear_ratio`` times faster divides it by that ratio and
        leaves the power untouched.  Nothing physical depends on which shaft
        the reading is taken from; only the numbers on the datasheet do.
        """
        return float(torque) / self.gear_ratio


@dataclass(frozen=True)
class DesignTargets:
    """Right-hand sides of the equality and inequality constraints.

    The performance requirements the engine must meet, plus the well-posedness
    constraints that keep the linkage analysable (see :mod:`exlink.kinematics`).
    """

    expansion_stroke: float = 74.0
    """Required expansion stroke ``STE`` [mm] (equality constraint)."""

    compression_ratio: float = 16.0
    """Required compression ratio ``epsilon`` (equality constraint)."""

    max_rod_angle: float = 10.0
    """Upper bound on ``mra``, the piston-rod tilt away from vertical [deg]."""

    max_transmission: float = 0.985
    """Upper bound ``C`` on ``W = max(delta_c1, delta_c2)``.

    ``C = 0.9848`` keeps the transmission angle ``T`` inside
    ``[10 deg, 170 deg]``; it is rounded to 0.985 here.
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

    ``eta(X) = 0``, ``H(X) = 1000``, ``B(X) = 1000`` are substituted whenever
    the kinematic compatibility conditions (4a)/(6a) fail, or whenever
    ``lambda(theta_1)`` does not have four monotone phases.
    """

    efficiency: float = 0.0
    height: float = 1000.0
    width: float = 1000.0


DEFAULT_SPEC = EngineSpec()
DEFAULT_TARGETS = DesignTargets()
DEFAULT_PENALTY = PenaltyValues()
