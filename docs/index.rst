Range-optimal design of an extended-expansion engine linkage
============================================================

.. rubric:: Abstract

Linkage design problems are conventionally posed as the maximisation of a
kinematic quality measure subject to bounds on envelope dimensions, with
performance requirements imposed as equality constraints. Applied to an
extended-expansion (Atkinson) engine mechanism, that formulation is shown to be
inadequate in three respects: it establishes no exchange rate between its
competing objectives and so admits no solution, only a front; its objective is
identically equal to the work done by the gas under the virtual-work identity
and therefore quantifies no loss; and it determines no cross-section, so the
inertia loads that govern the design at operating speed cannot be represented
within it.

An alternative formulation is developed in which the objective is the
application-level quantity of merit -- distance travelled per unit of fuel
consumed -- and in which the competing geometric quantities are priced by the
physical chain connecting them to it. The resulting problem is simultaneously
multidisciplinary, mixed-integer, and, as specified, possessed of a feasible set
of Lebesgue measure zero. Each property is shown to determine an element of the
admissible solution method rather than merely to complicate it. The measure-zero
feasible set, in particular, obliges the relaxation of the equality requirements
into tolerance bands; those bands are then shown to be comparable in width with
the manufacturing scatter of the components, so that reliability-based
optimization arises as a consequence of the relaxation rather than as an
addition to the problem.

Five quantitative results follow. The quasi-statically optimal geometry
coincides with the transmission-angle singularity at which the inertia loads are
largest and admits no feasible structure above 2000 rev/min, whereas a geometry
displaced from it attains lower mass and greater range simultaneously. A
tolerance study conducted against the specification rather than against a design
identifies which of the stated bounds govern reliability and prices their
relaxation at 0.47 % of the objective in exchange for a fall in the probability
of violating a dimensional requirement from 0.664 to 1.9e-5. Power is taken from
the shaft that turns twice per cycle, as on any four-stroke, so the mechanism's
advantage over a conventional slider-crank -- sized by identical structural and
tribological models, and optimised over its own degrees of freedom rather than
proportioned by hand -- is measured at equal speed and equal firing rate and
amounts to 17.6 %; it decomposes into an indicated efficiency of 0.480 against
0.457, a mechanical efficiency of 0.865 against 0.787 and an engine mass of
12.9 kg against 16.9, so that extended expansion, the feature the topology
exists for, is the smallest of the three. Imposing the coupled and structural
constraints throughout the search rather than verifying them on its result is
worth 4.9 % of the objective under an identical algorithm, the deterministic
optima of both mechanisms being shown to be dominated by designs standing
slightly off their active constraints. Finally, the first-order reliability
model is shown to fail in two independent ways that no output of it signals: it
prices only what its uncertain vector contains, so that the constraint governing
the design -- a gear face width at a reliability index of zero -- lies outside
the model entirely; and it linearises the expansion stroke at a tie between the
two top dead centres that the optimizer itself has closed to 0.107 um, giving a
failure probability seven times too small against 150 000 sampled builds. Both
corrections are free.

The topology is finally removed from the assumptions altogether: a
spring-connected synthesis, in which every candidate member is a spring whose
stiffness is a design variable and is penalised towards rigid or absent, is run
on the same requirements with no mechanism given. It returns discrete linkages
and no alternative. The one start reaching a four-stroke motion returns an Otto
engine whose two strokes are equal to six figures, for a reason that is an
identity rather than an accident -- a piston reached from one shaft alone is
periodic in that shaft's angle, so the asymmetry, which lives entirely in the
first harmonic, is exactly zero. Extended expansion therefore requires the piston
train to be reached from both shafts, which is what the three-cornered link of
the studied mechanism does.

The methodological contribution is a demonstration, on a problem small enough to
be verified throughout, that the choice of objective and the geometry of the
feasible set determine both the design obtained and the class of algorithms that
may legitimately be applied to obtain it.

.. rubric:: Keywords

Multidisciplinary design optimization; mixed-integer nonlinear programming;
reliability-based design optimization; mechanism synthesis; analytic
sensitivities; extended-expansion engine.

.. rubric:: Contents

.. toctree::
   :maxdepth: 2
   :caption: The study

   introduction
   state_of_the_art
   formulation
   method
   results
   conclusions
   references

.. toctree::
   :maxdepth: 2
   :caption: Appendices

   theory
   implementation
   supporting
   api

.. rubric:: At a glance

.. code-block:: python

   from exlink import COUPLED_DESIGN, evaluate
   from exlink.robustness import failure_probability, format_reliability

   # speed_rpm is the half-speed shaft's; the crankshaft turns twice per cycle
   outcome = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
   print(outcome.output_speed_rpm)      # 2000.0, the speed the study quotes
   print(outcome.km_per_litre)          # 3338.3
   print(outcome.budget.kilograms())    # where the mass actually is

   print(format_reliability(failure_probability(COUPLED_DESIGN)))

.. rubric:: Indices

* :ref:`genindex`
* :ref:`modindex`
