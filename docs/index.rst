.. doctest-skip-all

##################
sbpy Documentation
##################

`sbpy` is an `Astropy`_ affiliated package for small-body planetary astronomy.
       
For an overview on the expected structure and functionality of `sbpy`, please
refer to the :doc:`about` page; the :doc:`status` page provides an overview on
the implementation status of all modules and functions.

The initial development of `sbpy` is expected to conclude with version 1.0 in
2025.

	       
**************
Current Status
**************

.. image:: https://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: https://www.astropy.org
    :alt: Powered by Astropy Badge

.. image:: https://readthedocs.org/projects/sbpy/badge/?version=latest
    :target: https://sbpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation status

.. image:: https://joss.theoj.org/papers/10.21105/joss.01426/status.svg
    :target: https://doi.org/10.21105/joss.01426
    :alt: JOSS documentation

.. image:: https://github.com/NASA-Planetary-Science/sbpy/actions/workflows/ci_cron_weekly.yml/badge.svg
    :target: https://github.com/NASA-Planetary-Science/sbpy/actions
    :alt: GitHub testing status

.. image:: https://codecov.io/gh/NASA-Planetary-Science/sbpy/branch/main/graph/badge.svg
    :target: https://app.codecov.io/gh/NASA-Planetary-Science/sbpy
    :alt: codecov status

	     
***************
Getting Started
***************

.. toctree::
   :maxdepth: 2

   about.rst
   install.rst
   Tutorials <https://github.com/NASA-Planetary-Science/sbpy-tutorial>
   status.rst


******************
User Documentation
******************

Data Structures: Orbits, Ephemerides, Observations, and Physical Properties
---------------------------------------------------------------------------

.. toctree::
   :maxdepth: 1
   :glob:

   sbpy/data/index.rst
   sbpy/data/fieldnames.rst
   sbpy/data/*

Photometry and Spectroscopy
---------------------------

.. toctree::
   :maxdepth: 2

   sbpy/calib.rst
   sbpy/spectroscopy/index.rst
   sbpy/spectroscopy/sources.rst
   sbpy/photometry.rst

Dynamics and time
-----------------

.. toctree::
   :maxdepth: 2

   sbpy/dynamics
   sbpy/time

Activity
--------

.. toctree::
   :maxdepth: 3

   sbpy/activity/index.rst

Miscellaneous
-------------

.. toctree::
   :maxdepth: 1

   sbpy/bib.rst
   sbpy/exceptions.rst
   sbpy/units.rst
   sbpy/utils.rst

External packages that have been modified as part of `sbpy`
-----------------------------------------------------------

* `pyoorb <https://github.com/oorb/oorb/tree/master/python>`__: additional functionality for ephemerides computation, orbit transformation, and orbit propagation
* `astroquery <https://github.com/astropy/astroquery>`__: added submodules ``jplhorizons``, ``jplsbdb``, ``jplspec``, ``imcce``, and modified ``mpc``
   
   
*********
Reference
*********

.. toctree::
   :maxdepth: 2

   api_reference.rst
   sbpy/data/fieldnames.rst


***********************
Developer Documentation
***********************

.. toctree::
   :maxdepth: 2

   development/index.rst
   development/design-principles.rst

   
***************  
Acknowledgments
***************

`sbpy` is supported by NASA Planetary Data Archiving, Restoration, and Tools
(PDART) Grant Numbers 80NSSC18K0987 and 80NSSC22K0143.

If you use `sbpy` in your work, please acknowledge it by citing

    `Mommert, Kelley, de-Val Borro, Li et al., (2019). sbpy: A Python module for small-body planetary astronomy. Journal of Open Source Software, 4(38), 1426 <https://joss.theoj.org/papers/10.21105/joss.01426>`_

and also please consider using the :doc:`bib </sbpy/bib>` reference tracking system to
properly acknowledge and reference the methods you used in the preparation of
your manuscript.

The core `sbpy` Team is introduced `here <https://sbpy.org/team.html>`__.
