"""Mechanical losses, and the efficiency that actually reaches the crankshaft.

The geometric problem carries a quantity called ``eta``, defined as a ratio of
the crank torque's work to the gas force's work.  With no friction in the model
those two works are equal at every crank angle -- that is the virtual-work
identity :mod:`exlink.loads` is verified against -- so ``eta`` is not an
efficiency in the thermodynamic sense at all.  It is a *kinematic quality
measure*: how much mean torque a given piston motion converts into, normalised
by the stroke.  Useful, but it cannot be the objective of an engine study,
because nothing is lost in it.

This module supplies the losses that are really there, from quantities the
dynamic solve already produces:

**Journal bearings.**  Every revolute joint carries a known reaction ``R`` and
rotates through a known relative angle.  A pin of radius ``r`` under boundary
or mixed lubrication dissipates ``dW = mu |R| r |d(theta_rel)|``.  Seven joints
contribute, and the two crank throws contribute the most because they turn
through a full revolution while the others only oscillate.

**Piston and rings.**  The liner reaction ``D`` is already solved for -- it is
the same quantity the side-load constraint ``gamma`` bounds.  The rings add a
roughly constant radial tension on top of it, which at the low loads of an
economy engine is the *dominant* term.  Both slide through ``|d(lambda)|``.

**Gear mesh.**  A well-cut spur pair loses one to two per cent of what it
transmits; this is carried as a flat mesh efficiency rather than modelled.

Why this closes the argument
-----------------------------
Friction is what makes the constraint set mean something.  Without it the
side-load bound ``gamma <= 0.02`` is an assertion about wear that never appears
in any objective, and the bearing-load bound is likewise inert.  With it, both
feed straight into the quantity being maximised: a design that leans on the
liner burns its fuel on the liner.

It is also what makes size cost something in the *right* place.  Large ``H``
and ``B`` buy a long lever arm and a flattering ``eta``, but they do it with
long members carrying large bearing reactions -- and now that shows up as
friction work as well as mass.

Model fidelity
--------------
This is a Coulomb model with constant coefficients, not a hydrodynamic one.  It
gets the *scaling* right (loss proportional to load and to sliding distance)
and it responds correctly to design changes, which is what an optimizer needs.
It will not predict absolute friction mean effective pressure to better than
about 30 %, and no conclusion here rests on the absolute value -- only on
comparisons between designs evaluated with the same coefficients.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .derivatives import ramp_derivative, spectral_derivative
from .dynamics import DynamicLoads
from .materials import FloatArray

JOURNAL_FRICTION = 0.008
"""Equivalent Coulomb coefficient at a hydrodynamic steel journal [-].

A plain bearing running on a full oil film has an effective coefficient of
0.005 to 0.015; 0.008 is a warm, mid-speed value.  It is *not* the boundary
figure of 0.04, which applies only at start-up -- using the boundary value
would predict an engine that cannot turn itself over, and this one does.

The absolute loss scales linearly with this number, so every conclusion drawn
from it is stated as a comparison between designs at the same coefficient.
:func:`exlink.friction.sensitivity` sweeps it.
"""

PISTON_FRICTION = 0.10
"""Coulomb coefficient between the piston and ring pack and the liner [-]."""

RING_TENSION = 30.0
"""Radial force the ring pack presses on the liner with, independent of gas load [N].

For a 32 mm bore with a two-compression-plus-oil pack.  At the low cylinder
pressures of an economy engine this term is larger than the gas-driven side
load, which is why it cannot be dropped.
"""

MESH_EFFICIENCY = 0.985
"""Fraction of transmitted power surviving one spur-gear mesh [-]."""

#: Relative rotation at each joint: the two bodies whose angles differ there.
JOINT_BODIES: dict[str, tuple[str, str]] = {
    "R1": ("crank_1", "frame"),
    "Q": ("crank_1", "swing_rod"),
    "A": ("swing_rod", "trigonal"),
    "D": ("trigonal", "crank_2"),
    "R2": ("crank_2", "frame"),
    "E": ("trigonal", "piston_rod"),
    "P": ("piston_rod", "piston"),
}

#: Which sized members meet at each joint; the pin is as stout as the stoutest.
JOINT_MEMBERS: dict[str, tuple[str, ...]] = {
    "R1": ("crank_1",),
    "Q": ("crank_1", "swing_rod"),
    "A": ("swing_rod", "trigonal_ad", "trigonal_ae"),
    "D": ("trigonal_ad", "trigonal_de", "crank_2"),
    "R2": ("crank_2",),
    "E": ("trigonal_ae", "trigonal_de", "piston_rod"),
    "P": ("piston_rod",),
}


@dataclass(frozen=True)
class FrictionLosses:
    """Mechanical loss breakdown over one crankshaft revolution."""

    joint_work: dict[str, float]
    """Friction work at each journal [N.mm per revolution]."""

    piston_work: float
    """Friction work at the piston and rings [N.mm per revolution]."""

    mesh_work: float
    """Work lost in the gear mesh [N.mm per revolution]."""

    indicated_work: float
    """Work delivered to the crankshaft before friction [N.mm per revolution]."""

    @property
    def bearing_work(self) -> float:
        """Total journal friction work [N.mm per revolution]."""
        return float(sum(self.joint_work.values()))

    @property
    def total_work(self) -> float:
        """All mechanical losses [N.mm per revolution]."""
        return self.bearing_work + self.piston_work + self.mesh_work

    @property
    def brake_work(self) -> float:
        """Work leaving the crankshaft [N.mm per revolution]."""
        return self.indicated_work - self.total_work

    @property
    def mechanical_efficiency(self) -> float:
        """``W_brake / W_indicated``, the real mechanical efficiency [-].

        Negative or zero when friction exceeds the indicated work, which is a
        genuine outcome for a badly proportioned linkage at low load, not a
        numerical failure.  Callers should treat a non-positive value as an
        engine that will not run.
        """
        if self.indicated_work <= 0.0:
            return 0.0
        return self.brake_work / self.indicated_work

    @property
    def runs(self) -> bool:
        """Whether the engine produces net positive work."""
        return self.brake_work > 0.0

    def breakdown_kj(self) -> dict[str, float]:
        """Loss shares as fractions of indicated work, for reporting."""
        if self.indicated_work <= 0.0:
            return {}
        return {
            "bearings": self.bearing_work / self.indicated_work,
            "piston": self.piston_work / self.indicated_work,
            "mesh": self.mesh_work / self.indicated_work,
            "brake": self.brake_work / self.indicated_work,
        }


def _relative_rotation_rate(loads: DynamicLoads) -> dict[str, FloatArray]:
    """``|d(theta_rel) / d(theta_1)|`` at each joint.

    The body angles are taken from the mass properties, where they are already
    assembled; ``frame`` is the fixed reference and contributes nothing.
    """
    angles = dict(loads.mass_properties.body_angle)
    angles["frame"] = np.zeros_like(loads.kinematics.theta_1)
    rates: dict[str, FloatArray] = {}
    for joint, (first, second) in JOINT_BODIES.items():
        relative = angles[first] - angles[second]
        rates[joint] = np.abs(ramp_derivative(relative, 1))
    return rates


def pin_diameters(diameters: dict[str, float]) -> dict[str, float]:
    """Pin diameter at each joint, from the members that meet there.

    A pin is taken as stout as the stoutest member it joins.  This is a
    modelling choice, not a sizing calculation: the pins are not separately
    sized, and making them follow the members keeps the friction radius
    responding to the structural design in the right direction.

    Args:
        diameters: Sized member diameters [mm].

    Returns:
        ``{joint: diameter}`` [mm].
    """
    return {
        joint: max(float(diameters[name]) for name in members if name in diameters)
        for joint, members in JOINT_MEMBERS.items()
    }


def losses(
    loads: DynamicLoads,
    diameters: dict[str, float],
    journal_friction: float = JOURNAL_FRICTION,
    piston_friction: float = PISTON_FRICTION,
    ring_tension: float = RING_TENSION,
    mesh_efficiency: float = MESH_EFFICIENCY,
) -> FrictionLosses:
    """Integrate the mechanical losses over one crankshaft revolution.

    Args:
        loads: A solved dynamic load case, supplying joint reactions, the
            liner force and the crank torque.
        diameters: Sized member diameters [mm], for the pin radii.
        journal_friction: Coulomb coefficient at the journals.
        piston_friction: Coulomb coefficient at the liner.
        ring_tension: Ring-pack radial load [N].
        mesh_efficiency: Fraction of transmitted power surviving the mesh.

    Returns:
        The loss breakdown.  All works are per crankshaft revolution [N.mm].
    """
    n = loads.kinematics.theta_1.size
    step = 2.0 * np.pi / n
    rates = _relative_rotation_rate(loads)
    pins = pin_diameters(diameters)

    joint_work: dict[str, float] = {}
    for joint, rate in rates.items():
        reaction = np.linalg.norm(loads.reaction[joint], axis=1)
        radius = 0.5 * pins[joint]
        joint_work[joint] = float(journal_friction * radius * np.sum(reaction * rate) * step)

    # The piston slides through |d(lambda)| against the liner reaction plus the
    # ring tension, which is there whether or not the gas is pushing.
    slide = np.abs(spectral_derivative(loads.kinematics.lam, 1))
    normal = np.abs(loads.liner_force) + ring_tension
    piston_work = float(piston_friction * np.sum(normal * slide) * step)

    indicated = float(2.0 * np.pi * loads.mean_torque)
    mesh_work = max(indicated, 0.0) * (1.0 - mesh_efficiency)

    return FrictionLosses(
        joint_work=joint_work,
        piston_work=piston_work,
        mesh_work=mesh_work,
        indicated_work=indicated,
    )


def sensitivity(
    loads: DynamicLoads,
    diameters: dict[str, float],
    coefficients: tuple[float, ...] = (0.004, 0.006, 0.008, 0.012, 0.020),
) -> dict[float, float]:
    """Mechanical efficiency against the assumed journal coefficient.

    The single most uncertain number in the loss model is the journal
    coefficient, so its influence is reported rather than assumed away.  A
    conclusion that survives this sweep does not depend on the friction model
    being right in absolute terms.

    Args:
        loads: A solved dynamic load case.
        diameters: Sized member diameters [mm].
        coefficients: Journal coefficients to evaluate.

    Returns:
        ``{coefficient: mechanical_efficiency}``.
    """
    return {
        mu: losses(loads, diameters, journal_friction=mu).mechanical_efficiency
        for mu in coefficients
    }
