"""Reference designs: one historical baseline and three results.

:data:`PUBLISHED_DESIGN` is a **historical baseline** -- a design vector carried
over from an earlier, unpublished MATLAB study of this mechanism (Universite de
Technologie de Compiegne, 2015) by the same author.  It is kept as a *starting
point*, not as a result, and everything needed to use it is stated here.

Why it is a starting point and not a result
-------------------------------------------
Re-analysed as printed it violates five constraints.  Most tellingly it gives a
top-dead-centre gap of ``g = 8.5 mm`` against a 0.01 mm bound, and
``STE = 79.6 mm`` against a required 74 mm.  Since ``g`` is precisely the
quantity such an optimization drives to zero, this vector cannot be the design
that produced the properties recorded alongside it in that study.

Rounding does not explain it.  The vector is given to four significant figures,
and perturbing each variable by its rounding half-width moves ``STE`` by less
than 0.2 mm.

The real reason is conditioning, and it generalises.  The vector sits at
``W = 0.982``, a hair from the transmission-angle singularity at ``W = 1``, where
the piston motion is violently sensitive to the link lengths.  A design in that
region is not reproducible to four digits across two independent
implementations -- a property of where the optimum sits, not of either
implementation.  That same proximity to the singularity is what
:mod:`exlink.dynamics` later shows to be the expensive place to be once inertia
enters the load path.

Sweeping the two crank angles about the baseline values does recover its
recorded figures closely -- ``theta_f = 86 deg``, ``theta_r = -30 deg`` gives
``W = 0.9817`` and ``mra = 6.7 deg`` -- so the *model* agrees; it is the vector
that is not reproducible.

The results
-----------
============================  ==========================================
:data:`REFINED_DESIGN`        augmented Lagrangian, ``eta = 27.87 %``
:data:`GRADIENT_DESIGN`       SLSQP with exact gradients, ``eta = 30.77 %``
:data:`COUPLED_DESIGN`        coupled sizing, ``0.234 kg`` of moving mass
:data:`RANGE_DESIGN`          maximum range, constraints checked at the end
:data:`RELIABLE_DESIGN`       **the study's result**, ``3395 km/L`` at
                              ``P_f = 1.3e-3``
============================  ==========================================

Speeds are quoted at the output shaft throughout, which turns twice per cycle
(:attr:`exlink.constants.EngineSpec.output_revolutions_per_cycle`); the
``speed_rpm`` the analysis functions take is half of it.
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
"""Historical baseline: the design vector from the earlier MATLAB study.

Infeasible when re-analysed -- see the module docstring.  Used as a starting
point for the optimizers, and as the "before" in every comparison here.
"""

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

``torque_pressure_ratio`` was recorded there as "94.46%".  With
``phi = M_r,ave / p_ave`` that quantity has the dimension of a volume and no
normalisation turning it into a percentage was given, so it is kept here for
completeness but never compared against.
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

It lands close to the *properties* recorded in the earlier study, which is
worth noting even though those numbers are not independently verifiable:

==========  ==========  ==========
quantity    this        earlier study
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

``d`` is not compared: no geometric construction for the clearance was recorded,
so :func:`exlink.metrics.cylinder_clearance` defines its own (see its
docstring).

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

SLSQP with the exact Jacobians of :mod:`exlink.jacobian`, on the single-objective
geometric problem (maximise ``eta`` subject to every constraint), started from
:data:`PUBLISHED_DESIGN`.  It reaches **eta = 30.77 %** against the 27.87 % of
:data:`REFINED_DESIGN`, and is feasible at every crank-angle resolution from 720
samples up.

It is not a strictly better engine, and the comparison should not be read that
way: with no limit on the envelope in this formulation, it buys the efficiency
partly with size, growing to ``H = 320 mm`` against 239 mm.  What it does show
is that a derivative-free polish stops well short of the efficiency available at
comparable feasibility -- three points of it -- and that the limit is the method
rather than the mechanism.

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


COUPLED_DESIGN = Design(
    a=85.841929,
    c=66.472210,
    I=57.784469,
    x_b=80.562685,
    y_b=-73.481770,
    x_1=-23.363749,
    e=102.869097,
    q_1=7.057353,
    q_2=19.807655,
    theta_f=86.385886,
    theta_r=-33.363656,
)
"""The lightest design the *coupled* problem yields, at 2000 rev/min.

SLSQP on the MDF formulation, with exact Jacobians through the sizing/dynamics
MDA (:mod:`exlink.dynamics_jacobian`), minimising total moving mass subject to
every constraint and a floor of 25 % on efficiency.

Against :data:`REFINED_DESIGN` at the same speed:

==================  =========  ==========
quantity            refined    coupled
==================  =========  ==========
total moving mass   1.039 kg   **0.234 kg**
peak bearing load   12 629 N   7 504 N
``H``               238.5 mm   205.7 mm
``W``               0.9811     0.9372
``eta``             28.20 %    25.00 %
==================  =========  ==========

Four times lighter, a third off the bearing load, and a smaller envelope, for
three points of efficiency -- and the mechanism has moved well away from the
transmission-angle singularity, which is where the mass was going.

Feasible at every crank-angle resolution from 360 to 2880 samples.  The search
held a small margin on the constraints SLSQP drives hardest (``g``, ``mra``,
``gamma``), because a gradient method converges *onto* its active bounds and
those metrics shift slightly with resolution.

Reproduce it with::

    exlink optimize-mass --rpm 1000 --gradients
"""

COUPLED_METRICS: dict[str, float] = {
    "efficiency": 0.24933,
    "compatibility": 0.93717,
    "rod_angle": 9.701,
    "tdc_gap": 0.00697,
    "expansion_stroke": 73.970,
    "compression_ratio": 15.970,
    "total_mass_kg": 0.2343,
    "peak_bearing_load": 7504.0,
}
"""Properties of :data:`COUPLED_DESIGN`, at ``samples=1440``, 2000 rev/min."""


RANGE_DESIGN = Design(
    a=76.0262,
    c=71.0370,
    I=57.6000,
    x_b=72.9272,
    y_b=-74.1048,
    x_1=-31.5411,
    e=127.6415,
    q_1=8.2337,
    q_2=22.6621,
    theta_f=87.4537,
    theta_r=-45.2847,
)
"""The best design the vehicle-level problem produced, and what it does not achieve.

SLSQP on ``neg_range`` with the gear pair pinned to ``m = 0.8``, ``z = 48``
(hence ``I = 57.6``) at 2000 rev/min, started from :data:`COUPLED_DESIGN`.  It
reaches **3388 km/L** against 3338, on an engine of 12.47 kg, and it satisfies
every inequality constraint including the top-dead-centre gap, at
``g = 0.0009 mm``.

It is nonetheless **not** reported as feasible, and the reason is worth stating
rather than smoothing over.

What it misses, and by how much
--------------------------------
The two equalities are relaxed to ``|STE - 74| <= 0.05`` and
``|epsilon - 16| <= 0.05`` for the coupled problem.  This design finishes at
``0.0501`` and ``0.0501`` -- outside the band by 1.5e-4 mm and 6.1e-5, which is
SLSQP stopping within its own convergence tolerance of the constraint it was
given.  For scale, :mod:`exlink.robustness` puts the machining standard
deviation of ``STE`` at 0.020 mm, some 130 times larger.  No real part would
distinguish this design from one exactly on the band.

Why it is not simply projected back
-------------------------------------
Because doing so breaks something else.
:func:`exlink.scenarios.project_onto_equalities` restores both equalities to
machine precision with a minimum-norm step of a few hundredths of a
millimetre -- and that step moves ``g`` from 0.0009 mm to 0.0201 mm, twice its
bound.

That is the same hypersensitivity three other measurements find:

============================================  ===========================
perturbation                                  effect on ``g``
============================================  ===========================
IT8 machining tolerance on the members        ``sigma = 0.013 mm``
snapping ``I`` 0.18 mm onto the gear lattice  ``0.003 -> 0.058 mm``
minimum-norm equality projection              ``0.0009 -> 0.0201 mm``
============================================  ===========================

against a bound of 0.01 mm.  All three are larger than the bound, so which
number is written down decides whether a design on the equality manifold exists
at all: at ``g <= 0.01`` the manifold and the feasible region intersect in a
sliver thinner than the Newton step that reaches the manifold, and at
``g <= 0.1`` -- 0.1 mm of dead-centre mismatch being 2.7 % of the clearance
volume, worth 0.47 % of the range -- the projected design is strictly feasible.

What to do about it
--------------------
Take the wider bound, which :mod:`exlink.robustness` shows is also where ``g``
stops governing the reliability of the mechanism.  The projected design is then
the result at 3388 km/L.  Under the 0.01 mm bound as written the best *strictly*
feasible design remains :data:`COUPLED_DESIGN` at 3338 km/L, and the 1.5 %
difference is the price of a bound set finer than the model can resolve.
"""


RELIABLE_DESIGN = Design(
    a=72.2352,
    c=71.4395,
    I=57.6000,
    x_b=70.2676,
    y_b=-71.2069,
    x_1=-31.4620,
    e=129.7043,
    q_1=9.0337,
    q_2=23.5170,
    theta_f=76.6957,
    theta_r=-52.5420,
)
"""The study's result: the best design that also meets a reliability target.

From :func:`exlink.synthesis.maximise_range_from_target` with ``beta_target =
3``, started at :data:`COUPLED_DESIGN`, with the gear pair pinned to ``m = 0.8``,
``z = 48`` and every constraint -- geometric, coupled, vehicle and reliability --
imposed at every step rather than checked at the end.  1352 evaluations, 61
minutes, stopped at the iteration cap, so it is a lower bound on what the
formulation reaches.

==============================  ==========
range                           3395 km/L
engine mass                     12.94 kg
output speed                    2000 rev/min (1000 cycles/min)
system ``P_f`` at IT8           1.3e-3
system ``beta``                 3.00, on its target
worst constraint                -2.2e-7
==============================  ==========

The bounds it was solved against are the relaxed ones :mod:`exlink.robustness`
identifies -- the top-dead-centre gap at 0.054 mm and both equality bands at
``+/- 0.15`` -- and the result is insensitive to the first of those: re-scored
at a 0.1 mm gap the failure probability moves from 1.344e-3 to 1.339e-3, because
the gap sits at ``beta = 8.1`` either way.  What binds is the band on the
expansion stroke.

It is the design the figures in the documentation are drawn from.
"""

RELIABLE_METRICS: dict[str, float] = {
    "efficiency": 0.251373,
    "expansion_stroke": 74.0957,
    "compression_ratio": 16.1179,
    "rod_angle": 9.67497,
    "compatibility": 0.936426,
    "tdc_gap": 1.06805e-04,
    "clearance": 55.5741,
    "side_load_ratio": 0.0180587,
    "height": 237.876,
    "width": 128.432,
}
"""Properties of :data:`RELIABLE_DESIGN` as shipped, at ``samples=1440``.

Of the design *rounded to four decimals*, which is how it is stored.  That
rounding moves ``g`` by 12 % -- from 9.6e-05 mm to 1.07e-04 -- and nothing else
by more than a part in ten thousand, which is :mod:`exlink.robustness`'s point
about that quantity restated on the result itself.  Both figures are four
orders of magnitude inside any bound the study considers, so it changes nothing
here; it would matter if ``g`` were still bounded at 0.01 mm.
"""
