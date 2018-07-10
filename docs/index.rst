sbpy Documentation
==================

`sbpy` is a community effort to build a Python package for small-body
planetary astronomy in the form of an `astropy`_ affiliated package.

Overview
--------

The goal of `sbpy` is to provide a standard library for algorithms
used by asteroid and comet researchers. The functionality of `sbpy`
will include the following:

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

Please note that this package is currently under heavy development.


Installation
------------

The current development version of `sbpy` can be obtained from `github
<https://github.com/mommermi/sbpy>`_ using

    ``git clone https://github.com/mommermi/sbpy.git``

This will create a new directory (``sbpy/``). In this directory, run

    ``python setup.py install``

in order to use `sbpy` in your default Python environment.

A `pip` installer will be provided with the first official release of
`sbpy` (presumably at the end of 2018).


Learning how to use `sbpy`
--------------------------

The `sbpy` team maintains a `tutorial repository
<https://github.com/NASA-Planetary-Science/sbpy-tutorial>`_ providing
tutorials and learning materials used in workshops. Upcoming workshops
are also announced on this website.
   

Current Status
--------------

This package is currently under heavy development. For an overview on
the expected structure and functionality of `sbpy`, please refer to
:doc:`structure` page; the :doc:`status` provides an overview on the
implementation status of all modules and functions.
	  
The current development version status is as follows:

.. image:: https://travis-ci.org/mommermi/sbpy.svg?branch=master
    :target: https://travis-ci.org/mommermi/sbpy
    :alt: Travis-CI status

.. image:: https://coveralls.io/repos/github/mommermi/sbpy/badge.svg?branch=master
    :target: https://coveralls.io/github/mommermi/sbpy?branch=master
    :alt: coverage
	 
.. image:: https://readthedocs.org/projects/sbpy/badge/?version=latest
    :target: http://sbpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


More information on `sbpy`
--------------------------

.. toctree::
   :maxdepth: 1  

   contributing.rst
   structure.rst
   status.rst

   
Reference/API
-------------

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



Acknowledgments
---------------

`sbpy` is supported by NASA PDART Grant No. 80NSSC18K0987.

If you use `sbpy` in your work, please acknowledge it using the following line:
   
   "*This work made use of sbpy (http://sbpy.org), a community-driven Python package for small-body planetary astronomy supported by NASA PDART Grant No. 80NSSC18K0987.*"

and also please consider using the `~sbpy.bib` reference tracking
system to properly acknowledge and reference the methods you used in
the preparation of your manuscript.
