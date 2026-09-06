"""Hollow sections, and the four shape factors that is all they are.

Every member in this study was a solid round bar, and section 7.2 named that
the largest single modelling conservatism: a tube of the same outer diameter
carries nearly the same bending for appreciably less mass, because the metal
near the neutral axis contributes to weight and almost nothing else.  That
conservatism is not neutral here.  Section 6.1's finding is a balance between
the load a member must carry and the inertia its own mass creates, so removing
mass from the members moves the balance -- and it moves it in the direction the
finding is most sensitive to.

What a tube changes, and what it does not
-----------------------------------------
Write ``k = d_i / d_o`` for the bore ratio and ``R^2`` for ``r_o^2 + r_i^2``.
Against a solid bar of the same outer diameter:

============  =========================  ==================
quantity      formula                    factor
============  =========================  ==================
area          ``pi d^2 (1 - k^2) / 4``   ``1 - k^2``
section mod.  ``pi d^3 (1 - k^4) / 32``  ``1 - k^4``
second mom.   ``pi d^4 (1 - k^4) / 64``  ``1 - k^4``
transverse I  ``m (3 R^2 + L^2) / 12``   ``1 + k^2`` on R^2
============  =========================  ==================

Four constants, and the whole of the change.  Nothing about the *form* of any
stress, utilisation or derivative moves: every quantity that was proportional to
a power of ``d`` is still proportional to the same power, scaled.  That is why
the analytic derivatives of :mod:`exlink.dynamics_jacobian` survive untouched
in their structure -- the logarithmic derivatives ``d(1/A)/dd = -2/(A d)`` and
``d(1/Z)/dd = -3/(Z d)`` are independent of the shape factor, which cancels.

Which members may be hollow
---------------------------
Not all of them.  The two crank throws are cantilevered off their shafts and
are, on any real engine, solid forgings continuous with the web; boring them
would be a manufacturing fiction.  The rods and the trigonal link are exactly
what is made from tube in practice.  So the bore ratio is applied by member
kind, and :data:`HOLLOW_KINDS` names the two kinds it applies to.

Why a wall thickness floor is needed
------------------------------------
With ``k`` fixed, the wall thickness ``d (1 - k) / 2`` shrinks with the
diameter, and a lightly loaded member sized down to a few millimetres would be
asked to have a wall of tenths.  :attr:`Section.min_wall` puts a floor under it
by raising the outer diameter of any member whose wall would fall below it,
which is the manufacturing constraint a real drawing would carry and is also
what stops the optimizer discovering free mass at small diameters.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

HOLLOW_KINDS: frozenset[str] = frozenset({"link", "truss"})
"""Member kinds a bore may be put through; see the module docstring."""

MAX_BORE_RATIO = 0.95
"""Largest ``d_i / d_o`` accepted [-].

Beyond this the section is foil and the beam theory used here -- which ignores
local buckling of the wall entirely -- stops being an approximation of anything.
Practical engine tubing sits near 0.7 to 0.85.
"""

DEFAULT_MIN_WALL = 1.5
"""Thinnest wall the members may be drawn with [mm]."""


@dataclass(frozen=True)
class Section:
    """The member cross-section shape, as a bore ratio and a wall floor."""

    bore_ratio: float = 0.0
    """``d_i / d_o`` for the members allowed a bore; zero is a solid bar."""

    min_wall: float = DEFAULT_MIN_WALL
    """Thinnest wall accepted [mm]."""

    def __post_init__(self) -> None:
        if not 0.0 <= self.bore_ratio <= MAX_BORE_RATIO:
            msg = f"bore ratio must lie in [0, {MAX_BORE_RATIO}], got {self.bore_ratio}"
            raise ValueError(msg)
        if self.min_wall <= 0.0:
            msg = f"the wall floor must be positive, got {self.min_wall}"
            raise ValueError(msg)

    @property
    def solid(self) -> bool:
        """Whether this is the solid round bar the study started from."""
        return self.bore_ratio == 0.0

    def ratios(self, kinds: Sequence[str]) -> FloatArray:
        """The bore ratio of each member, zero where a bore is not allowed.

        Args:
            kinds: One member kind per member, in the member order.

        Returns:
            ``(n_members,)`` bore ratios.
        """
        return np.asarray(
            [self.bore_ratio if kind in HOLLOW_KINDS else 0.0 for kind in kinds],
            dtype=np.float64,
        )

    def minimum_diameter(self, kinds: Sequence[str]) -> FloatArray:
        """Smallest outer diameter each member may be drawn at [mm].

        A wall of :attr:`min_wall` at the member's own bore ratio; zero
        constraint on a solid member, which has no wall.
        """
        ratios = self.ratios(kinds)
        with np.errstate(divide="ignore"):
            return np.where(ratios > 0.0, 2.0 * self.min_wall / (1.0 - ratios), 0.0)


SOLID = Section()
"""The solid round bar of sections 6.1 to 6.5."""

TUBULAR = Section(bore_ratio=0.6)
"""The tube section 6.6 measures, and it is a *measured* choice, not a habit.

Two effects fight over the bore ratio on this engine and the winner changes
sides part-way:

* raising ``k`` takes metal out of the members, which is the point;
* raising ``k`` also raises the wall floor's outer diameter, ``2 t / (1 - k)``,
  and past some ``k`` the floor -- not the load -- is what sizes the light
  members, so the section grows back.  At ``k = 0.75`` and a 1.5 mm wall no
  member may be drawn below 12 mm, and most of this linkage's are smaller than
  that when solid.

Optimising range over ``k`` and speed together lands at ``k = 0.579`` and
2276 rpm.  0.6 is that answer to one figure, and it is inside the 0.5 to 0.7
band tube stock is actually drawn in.
"""


def area_factor(ratios: FloatArray) -> FloatArray:
    """``1 - k^2``: cross-sectional area against a solid bar."""
    return 1.0 - np.asarray(ratios, dtype=np.float64) ** 2


def modulus_factor(ratios: FloatArray) -> FloatArray:
    """``1 - k^4``: section modulus, and second moment, against a solid bar."""
    return 1.0 - np.asarray(ratios, dtype=np.float64) ** 4


def gyration_factor(ratios: FloatArray) -> FloatArray:
    """``1 + k^2``: the radial term of the transverse inertia, per unit mass.

    A tube's mass sits further from its own axis than a bar's does, so per
    kilogram it has *more* transverse inertia even though it has less in
    absolute terms.
    """
    return 1.0 + np.asarray(ratios, dtype=np.float64) ** 2
