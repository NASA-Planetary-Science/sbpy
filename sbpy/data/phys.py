# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
SBPy data.Phys Module
=====================

Class for storing and querying physical properties

created on June 04, 2017
"""


from .core import DataClass

__all__ = ['Phys']


class Phys(DataClass):
    """Class for storing and querying physical properties"""

    @classmethod
    def from_horizons(cls, targetid, bib=None):
        """Load physical properties from JPL Horizons
        (https://ssd.jpl.nasa.gov/horizons.cgi)

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
        #>>> from sbpy.data import Phys
        #>>> phys = Phys.from_horizons('Ceres')

        not yet implemented

        """

    @classmethod
    def from_lowell(cls, targetid, bib=None):
        """Load physical properties from Lowell Observatory
        (http://asteroid.lowell.edu/).

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
        #>>> from sbpy.data import Phys
        #>>> phys = Phys.from_astorb('Ceres')

        not yet implemented

        """

    def derive_absmag(self):
        """Derive absolute magnitude from diameter and geometric albedo"""

    def derive_diam(self):
        """Derive diameter from absolute magnitude and geometric albedo"""
        
    def derive_pv(self):
        """Derive geometric albedo from diameter and absolute magnitude"""

    def derive_bondalbedo(self):
        """Derive Bond albedo from geometric albedo and photometric phase slope"""
