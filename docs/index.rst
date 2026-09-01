Range-optimal design of an extended-expansion engine linkage
============================================================

.. rubric:: Abstract

An extended-expansion (Atkinson) engine linkage for a Shell Eco-marathon car is
designed by multidisciplinary optimization with **range** -- distance on a given
quantity of fuel -- as the objective. The formulation conventionally used for
this mechanism maximises a lever-arm quality measure subject to envelope bounds;
it prices nothing, its central quantity is not an efficiency, and it cannot see
the parts. Range makes efficiency, envelope, torque ripple and structural mass
commensurable at rates the physics fixes.

Two structural properties then determine the methodology. The section sizes,
the masses and the inertia loads form a fixed point, making the problem
multidisciplinary in the strict sense, with a measured Gauss-Seidel contraction
factor rising from 0 at rest to 0.68 at 1500 rpm. And the feasible set, as
specified, has measure zero -- which rules out every sampling method, forces the
two equality requirements into tolerance bands, and, because those bands are
only 1.7 standard deviations wide against the scatter of the parts, obliges a
reliability-based rather than a deterministic treatment.

Three results follow. The design a quasi-static formulation prefers sits at the
transmission-angle singularity, exactly where inertia is worst: it has no
feasible structure above 1000 rpm, while a design backed off weighs half as much
and goes further. The specified top-dead-centre bound of 0.01 mm cannot be held
at any ISO grade, and the reference design has a 64.5 % probability of missing at
least one requirement. And the linkage's 24 % range advantage over a
slider-crank sized by identical code falls to 2.8 % once its one-revolution
firing frequency is removed, so the advantage is firing frequency rather than
extended expansion.

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
