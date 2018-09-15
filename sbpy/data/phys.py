# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Phys Module
=====================

Class for storing and querying physical properties

created on June 04, 2017
"""


from numpy import ndarray
from astroquery.jplsbdb import SBDB

from .core import DataClass, conf

__all__ = ['Phys']


class Phys(DataClass):
    """Class for storing and querying physical properties"""

    @classmethod
    def from_sbdb(cls, targetid):
        """Load physical properties from JPL Small-Body Database using
        `~astroquery.jplsbdb` for one target. Builds a `~Phys` object
        from the output of `'phys_par'` from SBDB.

        The current implementation of this function is limited to
        single object queries. Future implementations will enable the
        query of multiple objects in one call.

        Parameters
        ----------
        targetid : str, mandatory 
            Target identifier to be queried; use object numbers, names, or
            designations that as unambiguous as possible.

        Returns
        -------
        `~Phys` object

        Examples
        --------
        >>> from sbpy.data import Phys
        >>> ceres = Phys.from_sbdb('Ceres')
        >>> print(ceres['H'])  # doctest: +SKIP
           H    
        --------
        3.34 mag

        """
        sbdb = SBDB.query(targetid, phys=True)
        return cls.from_dict(sbdb['phys_par'])

    @classmethod
    def from_lowell(cls, targetid):
        """Load physical properties from Lowell Observatory
        (http://asteroid.lowell.edu/).

        The Lowell database will provide a database of physical
        properties which is a compilation of a number of different sources.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Phys # doctest: +SKIP
        >>> phys = Phys.from_astorb('Ceres') # doctest: +SKIP

        not yet implemented

        """
