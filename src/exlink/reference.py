"""The design published in the 2015 report, and notes on reproducing it.

The report closes with a table of the eleven design variables of its chosen
solution and a table of that solution's properties.  :data:`PUBLISHED_DESIGN`
is the first table; :data:`PUBLISHED_METRICS` is the second.

Reproducing the second table from the first
-------------------------------------------
It does not reproduce, and the reason is instructive rather than alarming.

The published design sits at ``W = 0.982``, i.e. a hair away from the
transmission-angle singularity ``W = 1`` that condition (4a) exists to keep it
clear of.  Near that boundary the piston motion is violently sensitive to the
link lengths: re-analysing the table as printed gives ``STE = 79.6 mm`` against
the reported 73.98, and -- the telling one -- a top-dead-centre gap of
``g = 8.5 mm`` against the reported 0.0069.  Since ``g`` is precisely the
quantity the optimizer drove to zero, a design that misses it by three orders of
magnitude cannot be the design that produced the second table.

The table is printed to four significant figures.  Perturbing each variable by
its rounding half-width moves ``STE`` by less than 0.2 mm, so rounding alone
does not explain the gap either; the printed values are a rounded record of a
design that the report's own MATLAB code would also fail to reproduce.

What *does* hold up is the model.  Sweeping the two crank angles about the
published values recovers the reported figures closely -- ``theta_f = 86 deg``,
``theta_r = -30 deg`` gives ``W = 0.9817`` (reported 0.9817) and
``mra = 6.7 deg`` (reported 7.55) -- and the force chain reproduces the
virtual-work torque to machine precision at every crank angle.

So this package treats the published table as a *starting point*, not as an
answer to be matched.  :data:`REFINED_DESIGN` is what the framework produces
when the augmented Lagrangian is run from it; it is genuinely feasible, and it
is what the tests and examples use.

One correction is folded into the code: the report's printed inversion of the
trigonal-link moment equation has a sign slip in the swing-rod load ``A``.  See
:mod:`exlink.loads`.
"""

from __future__ import annotations

from .design import Design

PUBLISHED_DESIGN = Design(
    a=112.7,
    c=91.1,
    I=57.0,
    x_b=103.0,
    y_b=-74.1,
    x_1=-18.8,
    e=156.2,
    q_1=8.0,
    q_2=23.1,
    theta_f=111.1,
    theta_r=-43.1,
)
"""The design vector tabulated at the end of the 2015 report."""

PUBLISHED_METRICS: dict[str, float] = {
    "efficiency": 0.2776,
    "compatibility": 0.9817,
    "torque_pressure_ratio": 0.9446,
    "rod_angle": 7.5499,
    "tdc_gap": 0.0069,
    "width": 156.2,
    "height": 256.7,
    "expansion_stroke": 73.98,
    "compression_ratio": 15.98,
    "side_load_ratio": 0.02,
    "clearance": 68.9,
}
"""The properties tabulated for that design.

``torque_pressure_ratio`` is reported as "94.46%".  The report defines
``phi = M_r,ave / p_ave``, which has the dimension of a volume; no normalisation
turning it into a percentage is given, so the figure is recorded here as
printed but is not compared against.
"""


REFINED_DESIGN = Design(
    a=109.892413,
    c=88.629613,
    I=56.547758,
    x_b=99.798999,
    y_b=-72.803270,
    x_1=-18.691005,
    e=137.158796,
    q_1=8.081368,
    q_2=23.008134,
    theta_f=106.385886,
    theta_r=-23.199421,
)
"""A genuinely feasible design, obtained from :data:`PUBLISHED_DESIGN`.

Produced by running :func:`exlink.scenarios.refine` (augmented Lagrangian,
Nelder-Mead sub-problems) from the published table, in successively tighter
boxes, with a small margin held on the active constraints so that the design
stays feasible at any crank-angle resolution.

It lands remarkably close to the report's own published *properties*, which is
the real evidence that this reconstruction is faithful:

==========  ==========  ==========
quantity    this        report 2015
==========  ==========  ==========
``eta``     27.87 %     27.76 %
``W``       0.9811      0.9817
``mra``     8.68 deg    7.55 deg
``g``       0.0060 mm   0.0069 mm
``STE``     74.000 mm   73.98 mm
``eps``     16.000      15.98
``gamma``   0.0181      0.02
``B``       151.9 mm    156.2 mm
``H``       238.6 mm    256.7 mm
==========  ==========  ==========

``d`` is not compared: the report gives the constraint but not the geometric
construction behind it, so :func:`exlink.metrics.cylinder_clearance` is a
reconstruction (see its docstring).

Reproduce it with::

    exlink refine --design published --save refined.json
"""

REFINED_METRICS: dict[str, float] = {
    "efficiency": 0.27874,
    "compatibility": 0.98115,
    "rod_angle": 8.679,
    "tdc_gap": 0.0060,
    "width": 151.86,
    "height": 238.55,
    "expansion_stroke": 74.000,
    "compression_ratio": 16.000,
    "side_load_ratio": 0.0181,
    "clearance": 57.65,
}
"""Properties of :data:`REFINED_DESIGN`, at ``samples=1440``."""


GRADIENT_DESIGN = Design(
    a=104.115553,
    c=95.061896,
    I=50.326843,
    x_b=84.708674,
    y_b=-87.482327,
    x_1=-19.922767,
    e=186.341296,
    q_1=7.777963,
    q_2=24.567851,
    theta_f=91.100000,
    theta_r=-38.844734,
)
"""What a gradient-based search finds, starting from the published table.

SLSQP with the exact Jacobians of :mod:`exlink.jacobian`, run on the report's
own single-objective problem (maximise ``eta`` subject to every constraint),
from :data:`PUBLISHED_DESIGN`.  It reaches **eta = 30.77 %** against the
report's 27.76 % and the 27.87 % of :data:`REFINED_DESIGN`, and it is feasible
at every crank-angle resolution from 720 samples up.

It is not a strictly better engine, and the comparison should not be read that
way: with no limit on the envelope in this formulation, it buys the efficiency
partly with size, growing to ``H = 320 mm`` against 239 mm.  What it does show
is that the report's augmented-Lagrangian step stopped well short of the
efficiency available to it -- three points of it -- and that the reason was the
method rather than the mechanism.

Reproduce it with::

    exlink optimize --design published --algorithm SLSQP --gradients
"""

GRADIENT_METRICS: dict[str, float] = {
    "efficiency": 0.30768,
    "compatibility": 0.98500,
    "rod_angle": 9.151,
    "tdc_gap": 0.00952,
    "width": 159.34,
    "height": 319.79,
    "expansion_stroke": 74.001,
    "compression_ratio": 16.001,
    "side_load_ratio": 0.0199,
    "clearance": 103.92,
}
"""Properties of :data:`GRADIENT_DESIGN`, at ``samples=2880``."""
