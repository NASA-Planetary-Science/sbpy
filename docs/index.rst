.. doctest-skip-all

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

`sbpy` requires Python 3.x - compatibility with Python 2.x is not
supported. We hence recommend that you install the latest version of
`Anaconda Python 3.x <https://www.anaconda.com/download/>`_ on your
system before installing `sbpy`. Make sure that Anaconda Python is
your default Python (this will be asked during the installation process).

Requirements
^^^^^^^^^^^^

Most requirements should be resolved during the installation process. However, we recommend to install the latest development version of `astroquery` using

    >>> pip install git+https://github.com/astropy/astroquery.git

Also, if you want to use `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_, you will have to
install it using the instructions provided on that page.
    

Using pip
^^^^^^^^^

The latest development version of `sbpy` can be easily installed using

    >>> pip install git+https://github.com/NASA-Planetary-Science/sbpy.git


Using GitHub
^^^^^^^^^^^^

This way of installing `sbpy` is recommended if you plan to contribute
to the module. The current development version of `sbpy` can be
obtained from `GitHub <https://github.com/mommermi/sbpy>`_ using

    >>> git clone https://github.com/mommermi/sbpy.git

This will create a new directory (``sbpy/``). In this directory, run

    >>> python setup.py install --user

in order to use `sbpy` in your default Python environment. If you plan to work on the code and always want to use the latest version of your code, you can install it with

    >>> python setup.py develop --user

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

.. image:: https://coveralls.io/repos/github/NASA-Planetary-Science/sbpy/badge.svg?branch=master
    :target: https://coveralls.io/github/NASA-Planetary-Science/sbpy?branch=master
    :alt: Coveralls status
	
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
