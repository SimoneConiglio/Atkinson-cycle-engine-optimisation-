"""Manufacturability: what a shop can actually cut, and what that costs in mass.

A sizing loop that reports ``d = 7.3184 mm`` has answered a question nobody
asked.  Round bar arrives in preferred sizes, gear cutters exist for standard
modules, and a wall thinner than a couple of millimetres cannot be cast or
reliably machined.  This module holds those rules in one place so that the rest
of the package can state a *buildable* design rather than a mathematical one.

Three kinds of rule appear here:

**Preferred sizes.**  :data:`STOCK_DIAMETERS` is the R20 renard series of
metric round bar, :data:`STANDARD_MODULES` the ISO 54 series 1 gear modules.
Both are ordered ascending, and :func:`round_up_to_stock` moves a continuous
requirement up to the next available size.  Rounding *up* is the only safe
direction: it can never turn a safe section unsafe.

**Minimum thicknesses.**  A dimension can be safe against every stress and
still be unmakeable.  :data:`MIN_WALL_THICKNESS` and friends are floors applied
after strength sizing, never before.

**The manufacturing premium.**  Rounding up adds mass that no stress analysis
asked for.  :func:`stock_premium` measures it, because it is the honest answer
to "how much did buildability cost?" and because it is not negligible: on a
member whose strength requirement lands just above a stock size the premium can
exceed 30 %.

Why this is applied *outside* the sizing fixed point
----------------------------------------------------
Rounding is a step function, and a step function inside a fixed-point iteration
turns a smooth contraction into something that can limit-cycle between two
stock sizes forever.  The coupled solve therefore converges on continuous
diameters, and discretisation is applied once at the end -- then the loads are
re-evaluated on the discretised sections to confirm the design is still safe.
That "size continuous, discretise, verify" order is standard practice and it is
what :func:`exlink.coupled.solve_coupled` does when asked for a buildable
result.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

#: Preferred metric round-bar diameters, R20 renard series [mm].
STOCK_DIAMETERS: FloatArray = np.array(
    [
        3.0,
        3.5,
        4.0,
        4.5,
        5.0,
        5.5,
        6.0,
        7.0,
        8.0,
        9.0,
        10.0,
        11.0,
        12.0,
        14.0,
        16.0,
        18.0,
        20.0,
        22.0,
        25.0,
        28.0,
        30.0,
        32.0,
        35.0,
        40.0,
        45.0,
        50.0,
        56.0,
        63.0,
        70.0,
        80.0,
        90.0,
        100.0,
    ]
)

#: ISO 54 series-1 gear modules [mm].  Cutters exist for these; a module
#: between two of them means a bespoke hob.
STANDARD_MODULES: FloatArray = np.array(
    [0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
)

MIN_WALL_THICKNESS = 2.5
"""Thinnest reliably castable aluminium wall [mm].

Sand and gravity die casting of a small crankcase; a thinner wall risks
misruns and cold shuts regardless of what the stress says.
"""

MIN_MACHINED_THICKNESS = 1.5
"""Thinnest wall to leave on a machined-from-solid part [mm].

Below this the part distorts as the clamping releases and the tolerance is
not holdable.
"""

MIN_CROWN_THICKNESS = 3.0
"""Thinnest piston crown [mm].

Set by thermal, not mechanical, considerations: the crown has to conduct
combustion heat to the rings, and a thin crown runs hot enough to pre-ignite.
"""

MIN_FACE_WIDTH = 4.0
"""Narrowest gear face a small hob will cut cleanly [mm]."""

MIN_TEETH = 17
"""Fewest teeth on a 20 deg involute pinion before undercutting [-].

The standard result: ``z_min = 2 / sin^2(alpha)`` rounded up, which is 17 at a
20 deg pressure angle.  Fewer teeth need profile shift, which this model does
not carry.
"""


def round_up_to_stock(
    diameter: FloatArray | float, catalogue: FloatArray = STOCK_DIAMETERS
) -> FloatArray:
    """Move each diameter up to the next available stock size.

    Args:
        diameter: Continuous requirement [mm], any shape.
        catalogue: Ascending preferred sizes [mm].

    Returns:
        The next catalogue entry at or above each input, elementwise.  A
        requirement above the largest entry is returned unchanged, so that a
        runaway design is still visible as a runaway rather than being silently
        clipped to the top of the catalogue.
    """
    d = np.asarray(diameter, dtype=float)
    index = np.searchsorted(catalogue, d, side="left")
    clipped = np.minimum(index, catalogue.size - 1)
    rounded = catalogue[clipped]
    return np.where(d > catalogue[-1], d, rounded)


def round_to_module(module: float, catalogue: FloatArray = STANDARD_MODULES) -> float:
    """Move a gear module to the nearest standard one.

    Unlike bar stock this rounds to *nearest*, not up: a module is not a safety
    margin, and the tooth strength is recovered by the face width instead.

    Args:
        module: Continuous module [mm].
        catalogue: Ascending standard modules [mm].

    Returns:
        The closest catalogue entry.
    """
    return float(catalogue[int(np.argmin(np.abs(catalogue - float(module))))])


def stock_premium(continuous: FloatArray | float, discrete: FloatArray | float) -> float:
    """Fractional mass added by rounding sections up to stock.

    Mass goes as ``d^2`` at fixed length, so the premium on one member is
    ``(d_stock / d_required)^2 - 1``.  Reported over a whole set of members as
    the ratio of summed areas, which weights each member by how much material
    it actually carries.

    Args:
        continuous: Strength-driven diameters [mm].
        discrete: The same members rounded to stock [mm].

    Returns:
        ``mass_discrete / mass_continuous - 1``, a fraction.  Zero means every
        member landed exactly on a stock size.
    """
    required = np.asarray(continuous, dtype=float)
    supplied = np.asarray(discrete, dtype=float)
    if required.size == 0:
        return 0.0
    return float(np.sum(supplied**2) / np.sum(required**2) - 1.0)
