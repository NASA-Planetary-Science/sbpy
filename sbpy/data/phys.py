# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Phys Module
=====================

Class for storing and querying physical properties

created on June 04, 2017
"""

from collections import OrderedDict

from numpy import ndarray, isnan
from astroquery.jplsbdb import SBDB

from .core import DataClass, conf

__all__ = ['Phys']


class Phys(DataClass):
    """Class for storing and querying physical properties"""

    @classmethod
    def from_sbdb(cls, targetids, references=False, notes=False):
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

        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        alldata = []
        for targetid in targetids:

            sbdb = SBDB.query(targetid, phys=True)

            print(sbdb['phys_par'])
            # print(type(sbdb['phys_par']['diameter_sig']))

            data = {}
            for key, val in sbdb['phys_par'].items():
                if val is None or val == 'None':
                    val = -99
                if '_note' in key:
                    if notes:
                        data[key] = val
                elif '_ref' in key:
                    if references:
                        data[key] = val
                else:
                    try:
                        if isnan(val):
                            val = None
                    except TypeError:
                        pass
                data[key] = val
            alldata.append(cls.from_dict(data))

            alldata[-1].add_column([sbdb['object']['fullname']],
                                   name='targetname', index=0)
            # print(alldata[-1].table['H_sig'][0])

        if len(alldata) <= 1:
            return alldata[0]
        else:
            out = alldata[0]
            for tab in alldata[1:]:
                print(type(out.table['diameter_sig'][0]),
                      type(tab['diameter_sig'][0]
                           ), out.table['diameter_sig'][0],
                      tab['diameter_sig'][0])
                out.add_rows(tab)
            return out

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
