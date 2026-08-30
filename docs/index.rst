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
