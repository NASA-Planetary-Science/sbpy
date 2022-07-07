=========================
Data Module (`sbpy.data`)
=========================

Introduction
------------

`sbpy.data` provides classes for dealing with all kinds of data that
planetary astronomers encounter: orbital elements
(`~sbpy.data.Orbit`), ephemerides (`~sbpy.data.Ephem`), observational
data (`~sbpy.data.Obs`), and physical properties
(`~sbpy.data.Phys`). These objects act as containers that bundle
similar data and have to be used in `sbpy` to maintain a transparent
and clean data flow between functions. Each of these classes is based
on the `~sbpy.data.DataClass` base class, which internally uses an
`~astropy.table.QTable` object in combination with some other
features. Hence, each of these objects can actually be thought of as a
table with different *fields* (or columns) and *rows*.

However, in contrast to simple tables, `~sbpy.data` objects provide
additional functionality that is detailed below.

Content
-------

.. toctree::
   :maxdepth: 2

   dataclass.rst
   ephem.rst
   obs.rst
   orbit.rst
   phys.rst
   names.rst


Reference/API
-------------
.. automodapi:: sbpy.data
    :no-heading:
