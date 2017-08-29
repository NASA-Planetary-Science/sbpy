A Python Module for Small-Body Planetary Astronomy
--------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

	  
`sbpy` is a community effort to build a Python package for small-body
planetary astronomy in the form of an `astropy`_ affiliated package.

The goal is to collect and implement well-tested and well-documented
code for the scientific study of asteroids and comets, including (but
not limited to):

* observation planning tools tailored to moving objects,
* photometry models for resolved and unresolved bodies,
* wrappers and tools for astrometry and orbit fitting,
* spectroscopy analysis tools and models for reflected light and emission
  from gas,
* cometary gas and dust coma simulation and analysis tools,
* asteroid thermal models for flux estimation and size/albedo estimation,
* image enhancement tools for comet comae and PSF subtraction tools,
* lightcurve and shape analysis tools, and
* access tools for various databases for orbital and physical data, as well as
  ephemerides services.


This package is under heavy development and currently provides only
very limited functionality. A funding proposal to support the
development of `sbpy` is currently pending.



Documentation
-------------

The official documentation is available `here`_.


Contributing
------------

If you are interested in contributing to `sbpy`, you can do that
through testing existing code (somewhat easy) or by contributing
Python code (somewhat hard). We appreciate every contribution and
acknowledge them by adding you to the list of `sbpy` authors.

Testing
~~~~~~~

We are looking for people to test our routines. In order to make this
task worthwhile for you, we only ask you to test routines that are
stable and which you can already use for your research. We will
indicate which functions are ready to be tested on this dedicated
status page (TBD). Keep in mind that there might be issues, so feel free
to compare the function results with your own results. If you find a
problem, please issue an issue report on this github page.

Coding
~~~~~~

If you would like to implement (or modify) some functionality
yourself, you are welcome to do so. Please check the wiki and the
documentation if such functionality already exists or is
planned. Before working on the code, please issue an issue report,
explaining your contribution and how it would fit into existing `sbpy`
functionality - please clearly state if existing functionality would
have to be modified. Once your contribution report has been accepted,
please fork the `sbpy` repository and work on that code. Once you are
done, please issue a pull request, linking to your original issue
report. Please follow the `astropy code of conduct`_.

	 

License
-------

This project is Copyright (c) sbpy team and licensed under the terms of the BSD 3-Clause license. See the licenses folder for more information.


.. _astropy: http://www.astropy.org/
.. _here: http://sbpy.readthedocs.io/en/latest/
.. _astropy conde of conduct: http://docs.astropy.org/en/latest/development/codeguide.html
