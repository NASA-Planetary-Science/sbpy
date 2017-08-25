SBPy Documentation
==================

SBPy - A Python Module for Small-Body Planetary Astronomy

**Please note that this package is under heavy development. The current documentation is only intended to provide an outline for the API to be used in SBPy.**


Module Structure
----------------

We foresee the following module structure for `sbpy`:

.. figure:: static/structure.png
   :alt: sbpy module structure	    

   `sbpy` design schematic. Modules are shown as rectangular boxes,
   important classes as rounded colored boxes. The left-hand side of
   the schematic is mainly populated with support modules that act as
   data containers and query functions. The right-hand side of the
   schematic shows modules that focus on small body-related
   functionality. Colored symbols match the colors and symbols of
   classes and modules they are using.

The current status of the individual modules and code elements can be
inquired from the this :ref:`status page`. 
   





API Outline
-----------

Note that this outline is not up-to-date. All structures are still in
development. Examples provide only a rough guideline to the intended
structure. 

.. toctree::
   :maxdepth: 1


   sbpy/data.rst	     
   sbpy/activity.rst
   sbpy/photometry.rst
   sbpy/shape.rst
   sbpy/spectroscopy.rst
   sbpy/imageanalysis.rst
   sbpy/thermal.rst
   sbpy/obsutil.rst
   sbpy/bib.rst
  
