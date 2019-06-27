.. doctest-skip-all

##################
sbpy Documentation
##################


`sbpy` is a community effort to build a Python package for small-body
planetary astronomy in the form of an `astropy`_ affiliated package.

.. Important:: Sbpy is functional, incomplete, and under development.  Expect API changes between v0.1 and v0.2.  However, the code's functionality will remain the same.  Starting with v0.2, we expect most modules will have a stable API.  For an overview on the expected structure and functionality of `sbpy`, please refer to :doc:`structure` page; the :doc:`status` provides an overview on the implementation status of all modules and functions.


***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   install.rst
   structure.rst
   whatsnew/0.1
   Tutorials <https://github.com/NASA-Planetary-Science/sbpy-tutorial>
   Get Help
   Report Problems <https://github.com/astropy/astropy/issues>
   About the Sbpy Project <https://sbpy.org/about.html>
   status.rst


******************
User Documentation
******************

Data Structures: Orbits, Ephemerides, and Physical Properties
-------------------------------------------------------------

.. toctree::
   :maxdepth: 1

   sbpy/data/index.rst
   sbpy/data/fieldnames.rst

Photometry and Spectroscopy
---------------------------

.. toctree::
   :maxdepth: 1

   sbpy/spectroscopy/index.rst
   sbpy/spectroscopy/sources.rst
   sbpy/photometry.rst

Activity
--------

.. toctree::
   :maxdepth: 1

   sbpy/activity/index.rst

Other Stuff
-----------

.. toctree::
   :maxdepth: 1

   sbpy/bib.rst
   sbpy/utils.rst

=======================
Developer Documentation
=======================

.. toctree::
   :maxdepth: 1

   development/index.rst

*************************
Other Project Information
*************************

External packages that have been modified as part of `sbpy`
-----------------------------------------------------------

* `pyoorb <https://github.com/oorb/oorb/tree/master/python>`__: additional functionality for ephemerides computation, orbit transformation, and orbit propagation
* `astroquery <https://github.com/astropy/astroquery>`__: added submodules ``jplhorizons``, ``jplsbdb``, ``jplspec``, ``imcce``, and modified ``mpc``

Acknowledgments
---------------

`sbpy` is supported by NASA PDART Grant No. 80NSSC18K0987.

If you use `sbpy` in your work, please acknowledge it using the following line:

   "*This work made use of sbpy (http://sbpy.org), a community-driven Python package for small-body planetary astronomy supported by NASA PDART Grant No. 80NSSC18K0987.*"

and also please consider using the `~sbpy.bib` reference tracking
system to properly acknowledge and reference the methods you used in
the preparation of your manuscript.
