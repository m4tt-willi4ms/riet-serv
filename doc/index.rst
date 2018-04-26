.. PROTO Rietveld Refinement documentation master file, created by
   sphinx-quickstart on Wed Sep 13 15:17:44 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

+++++++++++++++++++++++++++++++++++++
PROTO Rietveld Refinement
+++++++++++++++++++++++++++++++++++++

.. Welcome to PROTO Rietveld Refinement's documentation!
.. =====================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

Rietveld Phases Module
======================

This module is used to store input data and build representations
of a particular crystal phases's powder profile, based on some fitting
parameters. The module is meant to be extendable, so that more fine-grained models of powder profiles can be implemented at a later time.

      .. automodule:: src.RietveldPhases
      .. autoclass:: RietveldPhases
         :members:
         :member-order: bysource

Rietveld Refinery Module
========================

The module is designed to take in some list of Rietveld Phases and produce the
ingredients necessary to run the Rietveld Refinement engine. This engine can
then be run a number of times, and can be used to update plot data during the
refinement.

      .. automodule:: src.RietveldRefinery
      .. autoclass:: RietveldRefinery
         :members:
         :member-order: bysource
         :undoc-members:
         :inherited-members:
         :show-inheritance: