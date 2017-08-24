# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
================
SBPy Data Module
================

created on June 22, 2017
"""

from astropy.table import Table, Column
from astropy.time import Time
import astropy.units as u
import callhorizons
from .data import DataClass


class Phys():
    """Class for storing and querying physical properties"""

    @classmethod
    def from_horizons(cls, targetid, bib=None):
        """Load physical properties from `JPL Horizons`_

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Phys
        >>> phys = Phys.from_horizons('ceres'(

        not yet implemented

        .. _JPL Horizons: https://ssd.jpl.nasa.gov/horizons.cgi

        """

    @classmethod
    def from_lowell(cls, targetid, bib=None):
        """Load physical properties from `Lowell Observatory`_ 

        The Lowell database will provide a database of physical
        properties which is a compilation of a number of different sources.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        bib : SBPy Bibliography instance, optional, default None
            Bibliography instance that will be populated

        Returns
        -------
        Astropy Table

        Examples
        --------
        >>> from sbpy.data import Phys
        >>> phys = Phys.from_astorb('ceres'(

        not yet implemented

        .. _Lowell Observatory: http://asteroid.lowell.edu/ 

        """

    def derive_absmag(self):
        """Derive absolute magnitude from diameter and geometric albedo"""

    def derive_diam(self):
        """Derive diameter from absolute magnitude and geometric albedo"""
        
    def derive_pv(self):
        """Derive geometric albedo from diameter and absolute magnitude"""

    def derive_bondalbedo(self):
        """Derive Bond albedo from geometric albedo and photometric phase slope"""

