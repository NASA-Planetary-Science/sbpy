# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Phys Module
=====================

Class for storing and querying physical properties

created on June 04, 2017
"""

from collections import OrderedDict

from numpy import ndarray, array, isnan, nan, interp, log, exp
import astropy.units as u
from astroquery.jplsbdb import SBDB
from astroquery.jplspec import JPLSpec

from .core import DataClass
from ..bib import cite
from ..exceptions import SbpyException

__all__ = ['Phys']


class JPLSpecQueryFailed(SbpyException):
    '''
    Raise warning if molecular data query fails
    '''


class Phys(DataClass):
    """Class for storing and querying physical properties"""

    @classmethod
    @cite({'software: astroquery': '2019AJ....157...98G'})
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
        >>> phys = Phys.from_sbdb(['Ceres', '12893', '3552']) # doctest: +REMOTE_DATA
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

        # assemble data as Phys object
        return cls.from_columns(coldata, names=columnnames)

    @classmethod
    @cite({'software: astroquery': '2019AJ....157...98G'})
    def from_jplspec(cls, temp_estimate, transition_freq, mol_tag):
        """Returns relevant constants from JPLSpec catalog and energy
        calculations

        Parameters
        ----------
        temp_estimate : `~astropy.units.Quantity`
            Estimated temperature in Kelvins

        transition_freq : `~astropy.units.Quantity`
            Transition frequency in MHz

        mol_tag : int or str
            Molecule identifier. Make sure it is an exclusive identifier,
            although this function can take a regex as your molecule tag,
            it will return an error if there is ambiguity on what the
            molecule of interest is. The function
            `~astroquery.jplspec.JPLSpec.query_lines_async`
            with the option `parse_name_locally=True` can be used to parse
            for the exclusive identifier of a molecule you might be
            interested in. For more information, visit
            `astroquery.jplspec` documentation.

        Returns
        -------
        Molecular data : `~sbpy.data.Phys` instance
            Quantities in the following order from JPL Spectral Molecular Catalog:
                | Transition frequency
                | Temperature
                | Integrated line intensity at 300 K
                | Partition function at 300 K
                | Partition function at designated temperature
                | Upper state degeneracy
                | Upper level energy in Joules
                | Lower level energy in Joules
                | Degrees of freedom

        """

        if isinstance(mol_tag, str):
            query = JPLSpec.query_lines_async(
                min_frequency=(transition_freq - (1 * u.GHz)),
                max_frequency=(transition_freq + (1 * u.GHz)),
                molecule=mol_tag,
                parse_name_locally=True,
                get_query_payload=True)

            res = dict(query)
            # python request payloads aren't stable (could be
            # dictionary or list)
            # depending on the version, so make
            # sure to check back from time to time
            if len(res['Mol']) > 1:
                raise JPLSpecQueryFailed(
                    ("Ambiguious choice for molecule,\
                    more than one molecule was found for \
                    the given mol_tag. Please refine \
                    your search to one of the following tags\
                    {} by using JPLSpec.get_species_table()\
                    (as shown in JPLSpec documentation)\
                    to parse their names and choose your \
                    molecule of interest, or refine your\
                    regex to be more specific (hint '^name$'\
                    will match 'name' exactly with no\
                    ambiguity).").format(res['Mol']))
            else:
                mol_tag = res['Mol'][0]

        query = JPLSpec.query_lines(
            min_frequency=(transition_freq - (1 * u.GHz)),
            max_frequency=(transition_freq + (1 * u.GHz)),
            molecule=mol_tag)

        freq_list = query['FREQ']

        if freq_list[0] == 'Zero lines we':
            raise JPLSpecQueryFailed(
                ("Zero lines were found by JPLSpec in a +/- 1 GHz "
                 "range from your provided transition frequency for "
                 "molecule tag {}.").format(mol_tag))

        t_freq = min(list(freq_list.quantity),
                     key=lambda x: abs(x-transition_freq))

        data = query[query['FREQ'] == t_freq.value]

        df = int(data['DR'].data)

        lgint = float(data['LGINT'].data)

        lgint = 10**lgint * u.nm * u.nm * u.MHz

        elo = float(data['ELO'].data) / u.cm

        gu = float(data['GUP'].data)

        cat = JPLSpec.get_species_table()

        mol = cat[cat['TAG'] == mol_tag]

        temp_list = cat.meta['Temperature (K)'] * u.K

        part = list(mol['QLOG1', 'QLOG2', 'QLOG3', 'QLOG4', 'QLOG5', 'QLOG6',
                        'QLOG7'][0])

        temp = temp_estimate

        f = interp(log(temp.value), log(
            temp_list.value[::-1]), log(part[::-1]))

        f = exp(f)

        partition = 10**(f)

        part300 = 10 ** (float(mol['QLOG1'].data))

        # yields in 1/cm
        energy = elo + (t_freq.to(1/u.cm, equivalencies=u.spectral()))

        energy_J = energy.to(u.J, equivalencies=u.spectral())
        elo_J = elo.to(u.J, equivalencies=u.spectral())

        quantities = [t_freq, temp, lgint, part300, partition, gu, energy_J,
                      elo_J, df, mol_tag]

        names = ['t_freq', 'temp', 'lgint300', 'partfn300', 'partfn',
                 'dgup', 'eup_J', 'elo_J', 'degfreedom', 'mol_tag']

        # names = ('Transition frequency',
        #          'Temperature',
        #          'Integrated line intensity at 300 K',
        #          'Partition function at 300 K',
        #          'Partition function at designated temperature',
        #          'Upper state degeneracy',
        #          'Upper level energy in Joules',
        #          'Lower level energy in Joules',
        #          'Degrees of freedom', 'Molecule Identifier')

        result = cls.from_dict(dict(zip(names, quantities)))

        return result
