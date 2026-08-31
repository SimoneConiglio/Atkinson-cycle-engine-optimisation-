exlink
======

Multidisciplinary design of an EX-link Atkinson-cycle engine mechanism for a
Shell Eco-marathon car, optimised for the only thing that competition scores:
distance on a given quantity of fuel.

.. toctree::
   :maxdepth: 2

   theory

Start here
----------

.. code-block:: python

   from exlink import analyse, PUBLISHED_DESIGN
   from exlink.scenarios import format_analysis, refine

   print(format_analysis(analyse(PUBLISHED_DESIGN)))
   outcome = refine(PUBLISHED_DESIGN)

   # size the parts, with inertia in the load path
   from exlink import solve_for_design
   sized = solve_for_design(PUBLISHED_DESIGN, speed_rpm=1000.0)
   print(sized.total_mass_kg, sized.feasible)

   # carry it all the way to kilometres per litre
   from exlink import COUPLED_DESIGN, evaluate
   outcome = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
   print(outcome.km_per_litre, outcome.engine_mass_kg)
   print(outcome.budget.kilograms())

API
---

Physics
~~~~~~~

.. automodule:: exlink.design
.. automodule:: exlink.constants
.. automodule:: exlink.kinematics
.. automodule:: exlink.cycle
.. automodule:: exlink.loads
.. automodule:: exlink.metrics
.. automodule:: exlink.model

Sizing and dynamics
~~~~~~~~~~~~~~~~~~~

The second iteration: inertia in the load path, parts sized against
static, fatigue and buckling failure, and the fixed point that couples them.

.. automodule:: exlink.derivatives
.. automodule:: exlink.materials
.. automodule:: exlink.dynamics
.. automodule:: exlink.sizing
.. automodule:: exlink.coupled

Derivatives
~~~~~~~~~~~

Exact derivatives of the whole chain, which is what makes a gradient-based
optimizer applicable to a feasible set this thin.

.. automodule:: exlink.jacobian
.. automodule:: exlink.dynamics_jacobian

The vehicle-level problem
~~~~~~~~~~~~~~~~~~~~~~~~~

What turns three competing objectives into one: friction makes the constraint
set cost something, the mass budget converts the envelope and the torque ripple
into kilograms, and the vehicle converts kilograms and efficiency into range.

.. automodule:: exlink.friction
.. automodule:: exlink.gears
.. automodule:: exlink.manufacturing
.. automodule:: exlink.mass_budget
.. automodule:: exlink.vehicle
.. automodule:: exlink.performance

Study of the result
~~~~~~~~~~~~~~~~~~~

Whether the optimum survives manufacturing tolerance, how strongly coupled the
problem actually is, and whether the finding generalises past one linkage.

.. automodule:: exlink.robustness
.. automodule:: exlink.formulations
.. automodule:: exlink.slidercrank

Mixed-integer design
~~~~~~~~~~~~~~~~~~~~

The gear choice is discrete, and it pins the inter-axle distance.  This states
that as a mixed-integer program and solves it by bi-level outer approximation.
Needs the optional ``gemseo-bilevel-outer-approximation`` plugin
(``pip install exlink-opt[minlp]``).

.. automodule:: exlink.minlp

Optimization
~~~~~~~~~~~~

.. automodule:: exlink.disciplines
.. automodule:: exlink.scenarios

Visualisation
~~~~~~~~~~~~~

.. automodule:: exlink.plots
.. automodule:: exlink.animation

Command line
~~~~~~~~~~~~

.. automodule:: exlink.cli

Reference designs
~~~~~~~~~~~~~~~~~

.. automodule:: exlink.reference

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
