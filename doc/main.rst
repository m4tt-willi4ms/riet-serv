Main Modules
++++++++++++

Rietveld Phases Module
======================

This module is used to store input data and build representations
of a particular crystal phases's powder profile, based on some fitting
parameters. The module is meant to be extensible, so that more fine-grained models of powder profiles can be implemented at a later time.

      .. automodule:: src.rietveld_phases
      .. autoclass:: RietveldPhases
         :members:
         :member-order: bysource

Rietveld Refinery Module
========================

The module is designed to take in some list of Rietveld Phases and produce the
ingredients necessary to run the Rietveld Refinement engine. This engine can
then be run a number of times, and can be used to update plot data during the
refinement.

      .. automodule:: src.rietveld_refinery
      .. autoclass:: RietveldRefinery
         :members:
         :member-order: bysource
         :undoc-members:
         :inherited-members:
         :show-inheritance: