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
the manufacturing scatter of the components, so that a deterministic solution
lying within a band is not a design that satisfies the requirement.
Reliability-based optimization thereby arises as a consequence of the relaxation
rather than as an addition to the problem.

Three quantitative results follow. The quasi-statically optimal geometry
coincides with the transmission-angle singularity at which the inertia loads are
largest, and admits no feasible structure above 1000 rev/min, whereas a geometry
displaced from the singularity attains lower mass and greater range
simultaneously. The specified top-dead-centre tolerance of 0.01 mm is
unattainable at any ISO 286 machining grade, and the reference design exhibits a
0.645 probability of violating at least one of its dimensional requirements. The 15.6 % range
advantage of the mechanism over a conventional slider-crank -- sized by
identical structural and tribological models, and optimised over its own
degrees of freedom rather than proportioned by hand -- becomes a 4.3 % deficit
when its firing-frequency difference is removed, locating the advantage in the
cycle rate rather than in extended expansion and reversing its sign.

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
   methodology
   framework
   use_case
   results
   conclusions
   references

.. toctree::
   :maxdepth: 2
   :caption: Reference

   implementation
   api
   theory

.. rubric:: At a glance

.. code-block:: python

   from exlink import COUPLED_DESIGN, evaluate
   from exlink.robustness import failure_probability, format_reliability

   outcome = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
   print(outcome.km_per_litre)          # 3338.3
   print(outcome.budget.kilograms())    # where the mass actually is

   print(format_reliability(failure_probability(COUPLED_DESIGN)))

.. rubric:: Indices

* :ref:`genindex`
* :ref:`modindex`
