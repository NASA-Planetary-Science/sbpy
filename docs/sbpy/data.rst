Data Module (`sbpy.data`)
=========================

Introduction
------------

`sbpy.data` provides classes for dealing with orbital elements
(`sbpy.data.Orbit`), ephemerides (`sbpy.data.Ephem`), and physical
properties (`sbpy.data.Phys`). `Ephem`, `Orbit`, and `Phys` objects
act as containers for such parameters and can (and should) be used to
provide these to functions in sbpy. Each class is based on an
`astropy.Table`, providing the same functionality and features.

`sbpy.data` also provides additional interfaces to a number of
different services. Finally, `sbpy.data.Names` provides functions
related to naming conventions for asteroids and comets.


Reference/API
-------------
.. automodapi:: sbpy.data
    :no-heading:
		


 
