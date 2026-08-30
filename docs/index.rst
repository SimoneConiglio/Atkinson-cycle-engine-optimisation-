exlink
======

Multi-objective optimization of an EX-link Atkinson-cycle engine mechanism,
a Python and GEMSEO reconstruction of a 2015 study.

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

The iteration the report defers: inertia in the load path, parts sized against
static, fatigue and buckling failure, and the fixed point that couples them.

.. automodule:: exlink.derivatives
.. automodule:: exlink.materials
.. automodule:: exlink.dynamics
.. automodule:: exlink.sizing
.. automodule:: exlink.coupled

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
