# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy time
=========

Solar System time keeping.

"""

__all__ = [
    "SpiceEphemerisTime",
]

from astropy.time import TimeFromEpoch


class SpiceEphemerisTime(TimeFromEpoch):
    """Number of seconds since J2000.0 epoch.

    This is equivalent to the NAIF SPICE Ephemeris Time when the Barycentric
    Dynamical Time (TDB) scale is used.

    Examples
    --------
    >>> from astropy.time import Time
    >>> import sbpy.time
    >>>
    >>> t = Time("2023-11-21")
    >>> print(t.et)
    753796869.1828325

    >>> t = Time(0, format="et")
    >>> print(t.iso)
    2000-01-01 12:00:00.000

    """

    name = "et"
    unit = 1.0 / 86400.0  # in days (1 day == 86400 seconds)
    epoch_val = "2000-01-01 12:00:00"
    epoch_val2 = None
    epoch_scale = "tdb"
    epoch_format = "iso"
