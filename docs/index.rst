sbpy Documentation
==================

sbpy - A Python Module for Small-Body Planetary Astronomy

`sbpy` is a community effort to build a Python package for small-body
planetary astronomy in the form of an `astropy`_ affiliated package.

The goal is to collect and implement well-tested and well-documented
code for the scientific study of asteroids and comets, including (but
not limited to):

* observation planning tools tailored to moving objects,
* photometry models for resolved and unresolved observations,
* wrappers and tools for astrometry and orbit fitting,
* spectroscopy analysis tools and models for reflected solar light and
  emission from gas,
* cometary gas and dust coma simulation and analysis tools,
* asteroid thermal models for flux estimation and size/albedo estimation,
* image enhancement tools for comet comae and PSF subtraction tools,
* lightcurve and shape analysis tools, and
* access tools for various databases for orbital and physical data, as well as
  ephemerides services.


Please note that this package is under heavy development. The current
documentation is only intended to provide an outline for the API to be
used in sbpy.

For an overview on the expected structure and functionality of `sbpy`,
please refer to :doc:`structure` page; the :doc:`status` provides an overview
on the implementation status of all modules and functions.


Content
-------
 
.. toctree::
   :maxdepth: 1  
	      
   structure.rst
   status.rst
   sbpy/index.rst


Current Status
--------------

.. image:: https://travis-ci.org/mommermi/sbpy.svg?branch=master
    :target: https://travis-ci.org/mommermi/sbpy
    :alt: Travis-CI status

.. image:: https://coveralls.io/repos/github/mommermi/sbpy/badge.svg?branch=master
    :target: https://coveralls.io/github/mommermi/sbpy?branch=master
    :alt: coverage
	 
.. image:: https://readthedocs.org/projects/sbpy/badge/?version=latest
    :target: http://sbpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
