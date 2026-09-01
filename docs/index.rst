exlink
======

Multidisciplinary design of an EX-link Atkinson-cycle engine mechanism for a
Shell Eco-marathon car, optimised for the only thing that competition scores:
distance on a given quantity of fuel.

Eleven geometric design variables, a discrete gear choice, a strongly coupled
structure/dynamics analysis solved as an MDA, and a chain carrying all of it
through to kilometres per litre -- so that efficiency, envelope size, torque
ripple and structural mass are priced against each other by physics rather than
by weights.

Three findings
--------------

**The quasi-static optimum is the worst place to be.**
   Maximising efficiency without inertia drives the linkage to its
   transmission-angle singularity, which is exactly where the accelerations, the
   bearing loads and hence the structure are worst.  Backing off costs nothing
   and halves the engine.  See :doc:`coupling`.

**A specified constraint cannot be manufactured.**
   The top-dead-centre gap is bounded at 0.01 mm, and the dimensions producing it
   scatter by more than that.  Four independent routes agree, and a
   reliability formulation puts the number on it: a 64.5 % chance of missing at
   least one requirement.  See :doc:`reliability`.

**The Atkinson linkage's advantage is firing frequency, not extended expansion.**
   Against a slider-crank sized by identical code it wins 24 %; remove the
   one-revolution cycle and the advantage falls to 2.8 %.  See :doc:`validation`.

Start here
----------

.. code-block:: python

   from exlink import analyse, COUPLED_DESIGN, evaluate

   # geometry, cycle, quasi-static loads, metrics
   print(analyse(COUPLED_DESIGN).metrics.efficiency)

   # the whole chain, to kilometres per litre
   outcome = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
   print(outcome.km_per_litre, outcome.engine_mass_kg)
   print(outcome.budget.kilograms())

Guide
-----

.. toctree::
   :maxdepth: 2
   :caption: The study

   problem
   formulation
   coupling
   gradients
   discrete
   reliability
   search
   validation
   results

.. toctree::
   :maxdepth: 2
   :caption: Reference

   theory
   api

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
