******************
Time (`sbpy.time`)
******************

Ephemeris time
==============

References to times in `sbpy` use the `astropy.time.Time` class.  For consistency with other tools, especially JPL Horizons and NAIF SPICE, `sbpy` defines an "Ephemeris Time" format as the number of seconds from the J2000 epoch in the TDB scale.  To use this enhancement, simply import the `sbpy.time` sub-module, and pass the format "et" when initializing a ``Time`` object:

.. doctest::

    >>> from astropy.time import Time
    >>> import sbpy.time
    >>>
    >>> j2000 = Time(0, format="et")
    >>> j2000.iso
    '2000-01-01 12:00:00.000'
    >>> Time(759307746.954761, format="et")
    <Time object: scale='tdb' format='et' value=759307746.954761>

Or, use the ``et`` property to transform dates to ephemeris time:

.. doctest::

    >>> Time("2010-11-04 13:59:47.31", scale="utc").et
    342151253.4925505

The conversion from UTC is thought to be good to the 0.1 ms level.  See the `sbpy` tests for more information.


Reference/API
=============

.. automodapi:: sbpy.time
    :no-main-docstr: