# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Phys Module
=====================

Class for storing and querying physical properties

created on June 04, 2017
"""

from collections import OrderedDict

from numpy import ndarray, array, isnan, nan
import astropy.units as u
from astroquery.jplsbdb import SBDB

from .core import DataClass
from .. import bib

__all__ = ['Phys']


class Phys(DataClass):
    """Class for storing and querying physical properties"""

    @classmethod
    def from_sbdb(cls, targetids, references=False, notes=False):
        """Load physical properties from `JPL Small-Body Database (SBDB)
        <https://ssd.jpl.nasa.gov/sbdb.cgi>`_ using
        `~astroquery.jplsbdb` for one or more targets. Builds a
        `~Phys` object from the output of `'phys_par'` from
        SBDB. Units are applied, where available. Missing data are
        filled up as nan values. Note that SBDB only serves physical
        properties data for a limited number of objects.

        Parameters
        ----------
        targetids : str, int or iterable thereof
            Target identifier(s) to be queried; use object numbers, names,
            or designations as unambiguous as possible.

        Returns
        -------
        `~Phys` object

        Examples
        --------
        >>> from sbpy.data import Phys
        >>> phys = Phys.from_sbdb(['Ceres', '12893', '3552'])
        >>> print(phys['targetname', 'H', 'diameter'])  # doctest: +SKIP
                targetname                 H          diameter
                                          mag            km
        -------------------------- ------------------ --------
                           1 Ceres               3.34    939.4
         12893 Mommert (1998 QS55)               13.9    5.214
        3552 Don Quixote (1983 SA) 12.800000000000002     19.0

        """

        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        alldata = []
        columnnames = ['targetname']
        columnunits = OrderedDict([('targetname', set())])
        for targetid in targetids:

            sbdb = SBDB.query(str(targetid), phys=True)

            # assemble data from sbdb output
            data = OrderedDict([('targetname', sbdb['object']['fullname'])])
            for key, val in sbdb['phys_par'].items():
                if val is None or val == 'None':
                    val = nan
                if '_note' in key:
                    if notes:
                        data[key] = val
                elif '_ref' in key:
                    if references:
                        data[key] = val
                else:
                    try:
                        if isnan(val):
                            val = nan
                    except TypeError:
                        pass
                data[key] = val

                # add to columnnames if not yet there
                if key not in columnnames:
                    columnnames.append(key)
                    columnunits[key] = set()

                # identify units
                if isinstance(val, u.Quantity):
                    columnunits[key].add(val.unit)
                elif isinstance(val, u.CompositeUnit):
                    for unit in val.bases:
                        columnunits[key].add(unit)

            alldata.append(data)

        # re-assemble data on a per-column basis
        coldata = []
        for col in columnnames:
            data = []

            for obj in alldata:
                try:
                    data.append(obj[col])
                except KeyError:
                    data.append(nan)

            # identify common unit (or at least any unit)
            try:
                unit = list(columnunits[col])[0]
                # transform data to this unit
                newdata = []
                for dat in data:
                    if isinstance(dat, (u.Quantity, u.CompositeUnit)):
                        try:
                            newdata.append(dat.to(unit))
                        except u.UnitConversionError:
                            # keep data untouched if conversion fails
                            unit = 1
                            newdata = data
                            break
                    else:
                        newdata.append(dat)
            except IndexError:
                # data has no unit assigned
                unit = 1
                newdata = data

            # convert lists of strings to floats, where possible
            try:
                data = array(newdata).astype(float)
            except (ValueError, TypeError):
                data = newdata

            # apply unit, if available
            if unit != 1:
                coldata.append(data*unit)
            else:
                coldata.append(data)

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Phys.from_sbdb',
                         {'data service url':
                          'https://ssd.jpl.nasa.gov/sbdb.cgi'})

        # assemble data as Phys object
        return cls.from_array(coldata, names=columnnames)

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
