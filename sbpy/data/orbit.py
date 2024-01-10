# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Orbit Module
======================

Class for querying, manipulating, integrating, and fitting orbital elements.

created on June 04, 2017

"""

__all__ = ['Orbit', 'OrbitError', 'OpenOrbError']

__doctest_requires__ = {
    ("Orbit.oo_propagate", "Orbit.oo_transform"): ["pyoorb", "astroquery"],
    ("Orbit.D_criterion", "Orbit.tisserand", "Orbit.from_horizons", "Orbit.from_mpc"): ["astroquery"],
}

import os
import itertools
from warnings import warn
import numpy as np
from numpy import array, ndarray, double, arange
from astropy.time import Time
from astropy.table import vstack, QTable
import astropy.units as u

# optional imports
try:
    from astroquery.jplhorizons import Horizons
    from astroquery.mpc import MPC
except ImportError:
    pass

try:
    import pyoorb
except ImportError:
    pass

from ..bib import cite, register
from ..exceptions import RequiredPackageUnavailable, SbpyException
from . import Conf, DataClass, QueryError, TimeScaleWarning
from ..utils.decorators import requires


class OrbitError(SbpyException):
    """Generic Error used in sbpy.data.orbit"""
    pass


class OpenOrbError(SbpyException):
    """To be raised by issues with OpenOrb"""
    pass


class Orbit(DataClass):
    """Class for querying, manipulating, integrating, and fitting orbital
    elements"""

    @classmethod
    @requires("astroquery")
    @cite({'data source': '1996DPS....28.2504G',
           'software: astroquery': '2019AJ....157...98G'})
    def from_horizons(cls, targetids, id_type='smallbody',
                      epochs=None, center='500@10',
                      **kwargs):
        """Load target orbital elements from
        `JPL Horizons <https://ssd.jpl.nasa.gov/horizons/>`_ using
        `astroquery.jplhorizons.HorizonsClass.elements`

        Parameters
        ----------
        targetids : str or iterable of str
            Target identifier, i.e., a number, name, designation, or JPL
            Horizons record number, for one or more targets.
        id_type : str, optional
            The nature of ``targetids`` provided; possible values are
            ``'smallbody'`` (asteroid or comet), ``'majorbody'`` (planet or
            satellite), ``'designation'`` (asteroid or comet designation),
            ``'name'`` (asteroid or comet name), ``'asteroid_name'``,
            ``'comet_name'``, ``'id'`` (Horizons id).
            Default: ``'smallbody'``
        epochs : `~astropy.time.Time` or dict, optional
            Epochs of elements to be queried; requires a
            `~astropy.time.Time` object with a single or multiple epochs. A
            dictionary including keywords ``start`` and ``stop``, as well
            as either ``step`` or ``number``, can be used to generate a range
            of epochs. ``start`` and ``stop`` have to be
            `~astropy.time.Time` objects (see :ref:`epochs`).
            If ``step`` is provided, a range
            of epochs will be queries starting at ``start`` and ending at
            ``stop`` in steps of ``step``; ``step`` has to be provided as
            a `~astropy.units.Quantity` object with integer value and a
            unit of either minutes, hours, days, or years. If
            ``number`` is provided as an integer, the
            interval defined by
            ``start`` and ``stop`` is split into ``number`` equidistant
            intervals. If ``None`` is
            provided, current date and time are
            used. All epochs should be provided in TDB; if not, they will be
            converted to TDB and a `~sbpy.data.TimeScaleWarning` will be
            raised.  Default: ``None``
        center : str, optional, default ``'500@10'`` (center of the Sun)
            Elements will be provided relative to this position.
        **kwargs : optional
            Arguments that will be provided to
            `astroquery.jplhorizons.HorizonsClass.elements`.

        Notes
        -----
        * For detailed explanations of the queried fields, refer to
          `astroquery.jplhorizons.HorizonsClass.elements` and the
          `JPL Horizons documentation <https://ssd.jpl.nasa.gov/horizons/manual.html>`_.
        * By default, elements are provided in the J2000.0 reference
          system and relative to the ecliptic and mean equinox of the
          reference epoch. Different settings can be chosen using
          additional keyword arguments as used by
          `astroquery.jplhorizons.HorizonsClass.elements`.

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='tdb')  # doctest: +REMOTE_DATA
        >>> eph = Orbit.from_horizons('Ceres', epochs=epoch)  # doctest: +REMOTE_DATA
        """

        # modify epoch input to make it work with astroquery.jplhorizons
        # maybe this stuff should really go into that module....
        if epochs is None:
            epochs = [Time.now().tdb.jd]
        elif isinstance(epochs, Time):
            if epochs.scale != 'tdb':
                warn(('converting {} epochs to tdb for use in '
                      'astroquery.jplhorizons').format(epochs.scale),
                     TimeScaleWarning)
            epochs = epochs.tdb.jd
        elif isinstance(epochs, dict):
            if 'start' in epochs and 'stop' in epochs and 'number' in epochs:
                epochs['step'] = epochs['number']*u.dimensionless_unscaled
            # convert to tdb and iso for astroquery.jplhorizons
            epochs['start'] = epochs['start'].tdb.iso
            epochs['stop'] = epochs['stop'].tdb.iso
            if 'step' in epochs:
                if epochs['step'].unit is not u.dimensionless_unscaled:
                    epochs['step'] = '{:d}{:s}'.format(
                        int(epochs['step'].value),
                        {u.minute: 'm', u.hour: 'h', u.d: 'd',
                         u.year: 'y'}[epochs['step'].unit])
                else:
                    epochs['step'] = '{:d}'.format(
                        int(epochs['step'].value-1))

        # if targetids is a list, run separate Horizons queries and append
        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        # append elements table for each targetid
        all_elem = None
        for targetid in targetids:

            # load elements using astroquery.jplhorizons
            obj = Horizons(id=targetid, id_type=id_type,
                           location=center, epochs=epochs)
            try:
                elem = obj.elements(**kwargs)
            except ValueError as e:
                raise QueryError(
                    ('Error raised by astroquery.jplhorizons: {:s}\n'
                     'The following query was attempted: {:s}').format(
                         str(e), obj.uri))

            # workaround for current version of astroquery to make
            # column units compatible with astropy.table.QTable
            # should really change '---' units to None in
            # astroquery.jplhorizons.__init__.py
            for column_name in elem.columns:
                if elem[column_name].unit == '---':
                    elem[column_name].unit = None

            if all_elem is None:
                all_elem = elem
            else:
                all_elem = vstack([all_elem, elem])

        # turn epochs into astropy.time.Time and apply timescale
        # https://ssd.jpl.nasa.gov/?horizons_doc
        all_elem['epoch'] = Time(all_elem['datetime_jd'], format='jd',
                                 scale='tdb')
        all_elem['Tp'] = Time(all_elem['Tp_jd'], format='jd',
                              scale='tdb')

        all_elem.remove_column('datetime_jd')
        all_elem.remove_column('datetime_str')
        all_elem.remove_column('Tp_jd')

        return cls.from_table(all_elem)

    @classmethod
    @requires("astroquery")
    @cite({'data source': 'https://minorplanetcenter.net/iau/MPEph/MPEph.html',
           'software: astroquery': '2019AJ....157...98G'})
    def from_mpc(cls, targetids, id_type=None, target_type=None, **kwargs):
        """Load latest orbital elements from the
        `Minor Planet Center <https://minorplanetcenter.net>`_.

        Parameters
        ----------
        targetids : str or iterable of str
            Target identifier, resolvable by the Minor Planet
            Ephemeris Service. If multiple targetids are provided in a list,
            the same format (number, name, or designation) must be used.

        id_type : str, optional
            ``'name'``, ``'number'``, ``'designation'``, or ``None`` to
            indicate
            type of identifiers provided in ``targetids``. If ``None``,
            automatic identification is attempted using
            `~sbpy.data.names`. Default: ``None``

        target_type : str, optional
            ``'asteroid'``, ``'comet'``, or ``None`` to indicate
            target type. If ``None``, automatic identification is
            attempted using
            `~sbpy.data.names`. Default: ``None``

        **kwargs : optional
            Additional keyword arguments are passed to
            `~astroquery.mpc.MPC.query_object`

        Returns
        -------
        `~Orbit` object
            The resulting object will be populated with most columns as
            defined in
            `~astroquery.mpc.query_object`; refer
            to that document on information on how to modify the list
            of queried parameters.


        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> orb = Orbit.from_mpc('Ceres') # doctest: +REMOTE_DATA
        >>> orb  # doctest: +SKIP
        <QTable length=1>
         absmag    Q      arc      w     ...     a        Tj   moid_uranus moid_venus
          mag      AU      d      deg    ...     AU                 AU         AU
        float64 float64 float64 float64  ...  float64  float64   float64    float64
        ------- ------- ------- -------- ... --------- ------- ----------- ----------
           3.34    2.98 79653.0 73.59764 ... 2.7691652     3.3     15.6642    1.84632
        """

        from ..data import Names

        # if targetids is a list, run separate Horizons queries and append
        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        for targetid in targetids:
            if target_type is None:
                target_type = Names.asteroid_or_comet(targetid)
            if id_type is None:
                if target_type == 'asteroid':
                    ident = Names.parse_asteroid(targetid)
                elif target_type == 'comet':
                    ident = Names.parse_comet(targetid)
                if 'name' in ident:
                    id_type = 'name'
                elif 'desig' in ident:
                    id_type = 'designation'
                elif 'number' in ident:
                    id_type = 'number'

        # append ephemerides table for each targetid
        all_elem = None
        for targetid in targetids:

            # get elements
            kwargs[id_type] = targetid
            e = MPC.query_object(target_type, **kwargs)

            # parse results from MPC.query_object
            results = {}
            for key, val in e[0].items():
                # skip if key not in Conf.mpc_orbit_fields
                if key not in Conf.mpc_orbit_fields:
                    continue

                fieldname, fieldunit = Conf.mpc_orbit_fields[key]
                # try to convert to float
                try:
                    val = float(val)
                except (ValueError, TypeError):
                    pass

                if fieldname == 'mpc_orbit_type':
                    results[fieldname] = [{
                        0: 'Unclassified', 1: 'Atira',
                        2: 'Aten', 3: 'Apollo', 4: 'Amor',
                        5: 'Mars Crosser', 6: 'Hungaria',
                        7: 'Phocaeas', 8: 'Hilda', 9: 'Jupiter Trojan',
                        10: 'Distant Object'}[int(val)]]
                elif fieldunit is None:
                    results[fieldname] = [val]
                elif fieldunit == 'time_jd_utc':
                    results[fieldname] = Time([val],
                                              scale='utc', format='jd')
                else:
                    if val is None:
                        continue
                    results[fieldname] = [val]*u.Unit(fieldunit)

            if all_elem is None:
                all_elem = QTable(results)
            else:
                all_elem = vstack([all_elem, QTable(results)])

        return cls.from_table(all_elem)

    # functions using pyoorb

    def _to_oo(self):
        """Converts this orbit object to a openorb-compatible orbit array

        Notes
        -----
        * Epochs must be provided as `~astropy.time.Time` objects with
          time scales that are compatible with `pyoorb`. If epochs are
          not provided appropriately, they will be adjusted and a
          `TimeScaleWarning` will be raised.
        """

        # identify orbit type based on available table columns
        if 'orbtype' in self.field_names:
            orbittype = self.table['orbtype'][0]
        else:
            orbittype = None

            for testtype in ['KEP', 'COM', 'CART']:
                field_names = [
                    field[0] for field in Conf.oorb_orbit_fields[testtype]
                ]
                try:
                    for field in field_names[1:6]:
                        self.__getitem__(field)
                    orbittype = testtype
                    break
                except KeyError:
                    pass

        if orbittype is None:
            raise OrbitError(
                'orbit type cannot be determined from elements')

        # implant ``targetname`` field information, if not available
        if 'targetname' not in self.field_names:
            self.table['targetname'] = ['orbit_'+str(i) for i in
                                        range(len(self.table))]

        # check that epochs are astropy.time.Time
        if not isinstance(self['epoch'][0], Time):
            raise OrbitError(
                'epochs have to be provided as astropy.time.Time objects')

        # check that pyoorb can deal with time scale
        if self['epoch'][0].scale.upper() not in ('UTC', 'UT1', 'TT', 'TAI'):
            warn(('epochs time scale is {} which is incompatible with '
                  'pyoorb; converting time scale to TT.').format(
                self['epoch'][0].scale.upper()))
            self['epoch'] = self['epoch'].tt

        # assemble orbit array for oorb_ephemeris
        if orbittype == 'COM':
            # cometary orbit: id q e i node argperi t_p otype epoch t H G
            orbits = array(array([arange(0, len(self.table), 1),
                                  self['q'].to('au').value,
                                  self['e'].data,
                                  self['i'].to_value('radian'),
                                  self['Omega'].to_value('radian'),
                                  self['w'].to_value('radian'),
                                  self['Tp'].mjd,
                                  [Conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  self['epoch'].mjd,
                                  ([Conf.oorb_timeScales
                                    [self['epoch'][0].scale.upper()]] *
                                   len(self.table)),
                                  self['H'].value,
                                  self['G'].data]).transpose(),
                           dtype=double, order='F')
        elif orbittype == 'KEP':
            # keplerian orbit: id a e i node argperi M otype epoch ttype H G
            orbits = array(array([arange(0, len(self.table), 1),
                                  self['a'].to('au').value,
                                  self['e'].data,
                                  self['incl'].to_value('radian'),
                                  self['Omega'].to_value('radian'),
                                  self['w'].to_value('radian'),
                                  self['M'].to_value('radian'),
                                  [Conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  self['epoch'].mjd,
                                  ([Conf.oorb_timeScales
                                    [self['epoch'][0].scale.upper()]] *
                                   len(self.table)),
                                  self['H'].value,
                                  self['G'].data]).transpose(),
                           dtype=double, order='F')
        elif orbittype == 'CART':
            # cartesian orbit: id x y z dx dy dz otype epoch ttype H G
            orbits = array(array([arange(0, len(self.table), 1),
                                  self['x'].to('au').value,
                                  self['y'].to('au').value,
                                  self['z'].to('au').value,
                                  self['dx'].to('au/d').value,
                                  self['dy'].to('au/d').value,
                                  self['dz'].to('au/d').value,
                                  [Conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  self['epoch'].mjd,
                                  ([Conf.oorb_timeScales
                                    [self['epoch'][0].scale.upper()]] *
                                   len(self.table)),
                                  self['H'].data,
                                  self['G'].data]).transpose(),
                           dtype=double, order='F')

        return orbits

    @staticmethod
    def _from_oo(oo_orbits, orbittype, timescale):
        """Convert openorb-compatible orbit array to ``Orbit``."""

        # reorder data in Orbit object
        fields = Conf.oorb_orbit_fields[orbittype]
        columns = []
        field_names = []
        for i, col in enumerate(oo_orbits.transpose()):
            field_name, field_unit = fields[i]
            column = Orbit._unit_apply(col, field_unit)

            # return units as degrees, not radians
            if field_unit == 'rad':
                column = np.degrees(column)

            # convert epoch and Tp to Time object, and convert back to user's
            # original time scale
            if field_name in ['epoch', 'Tp']:
                column = Time(
                    Time(column, format='mjd', scale='tt'),
                    scale=timescale.lower()
                )

            columns.append(column)
            field_names.append(field_name)

        return Orbit.from_columns(columns, names=field_names)

    @staticmethod
    def _from_oo_propagatation(oo_orbits, orbittype, timescale):
        """Convert openorb orbit array from oorb_propagation to ``Orbit``.

        pyoorb.oorb_propagation returns degrees, but elsewhere pyoorb returns
        radians.

        """

        # convert columns of degrees to radians, then _from_oo will convert
        # back to degrees
        fields = Conf.oorb_orbit_fields[orbittype]
        for i, (field_name, field_unit) in enumerate(fields):
            if field_unit == 'rad':
                oo_orbits[:, i] = np.radians(oo_orbits[:, i])

        return Orbit._from_oo(oo_orbits, orbittype, timescale)

    @requires("pyoorb")
    @cite({'method': '2009M&PS...44.1853G',
           'software': 'https://github.com/oorb/oorb'})
    def oo_transform(self, orbittype, ephfile='de430'):
        """Uses pyoorb to transform this orbit object to a different
        orbit type definition. Required fields are:

        * target identifier (``'targetname'``)
        * semi-major axis (``'a'``, for Keplerian orbit) or perihelion
          distance (``'q'``, for cometary orbit), typically in au or
          or x-component of state vector (``'x'``, for cartesian orbit),
          typically in au
        * eccentricity (``'e'``, for Keplerian or cometary orbit) or
          y-component of state vector (``'y'``, for cartesian orbit) in
          au
        * inclination (``'i'``, for Keplerian or cometary orbit) in
          degrees or z-component of state vector (``'z'``, for cartesian
          orbit) in au
        * longitude of the ascending node (``'Omega'``, for Keplerian or
          cometary orbit) in degrees or x-component of velocity vector
          (``'vx'``, for cartesian orbit), au/day
        * argument of the periapsis (``'w'``, for Keplerian or cometary
          orbit) in degrees or y-component of velocity vector (``'vy'``,
          for cartesian orbit) in au/day
        * mean anomaly (``'M'``, for Keplerian orbits) in degrees or
          perihelion epoch (``'Tp'``, for cometary orbits) as astropy Time
          object or z-component of velocity vector (``'vz'``, for cartesian
          orbit) in au/day
        * epoch (``'epoch'``) as `~astropy.time.Time` object
        * absolute magnitude (``'H'``) in mag
        * photometric phase slope (``'G'``)

        Parameters
        ----------
        orbittype : str
            Orbit definition to be transformed to; available orbit
            definitions are ``KEP`` (Keplerian elements), ``CART``
            (cartesian elements), ``COM`` (cometary elements).
        ephfile : str, optional
            Planet and Lunar ephemeris file version as provided by JPL
            to be used in the propagation. Default: ``'de430'``

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        Obtain the current state vector (cartesian definition, ``CART``) for
        asteroid Ceres.

        >>> from sbpy.data import Orbit
        >>> ceres = Orbit.from_horizons('Ceres')  # doctest: +REMOTE_DATA
        >>> statevec = ceres.oo_transform('CART') # doctest: +SKIP
        >>> print(statevec)  # doctest: +SKIP
        <QTable length=1>
           id           x                   y          ...    H       G    timescale
                        AU                  AU         ...   mag
          str7       float64             float64       ... float64 float64    str2
        ------- ------------------ ------------------- ... ------- ------- ---------
        1 Ceres -1.967176101061908 -1.7891189971612211 ...    3.34    0.12        TT
        """

        if pyoorb is None:
            raise RequiredPackageUnavailable('pyoorb')

        # initialize pyoorb
        if os.getenv('OORB_DATA') is None:
            # oorb installed using conda
            pyoorb.pyoorb.oorb_init()
        else:
            ephfile = os.path.join(os.getenv('OORB_DATA'), ephfile+'.dat')
            pyoorb.pyoorb.oorb_init(ephfile)

        # extract time scale
        timescale = self.table['epoch'].scale.upper()

        # derive and apply default units
        default_units = {}
        for field_name, field_unit in Conf.oorb_orbit_fields[orbittype]:
            try:
                # use the primary field name
                primary_field_name = self._translate_columns(field_name)[0]
            except KeyError:
                continue
            default_units[primary_field_name] = field_unit

        for colname in self.field_names:
            if (colname in default_units.keys() and
                not isinstance(self[colname],
                               (u.Quantity, u.CompositeUnit))):
                self[colname].unit = default_units[colname]

        # convert epochs to TT and MJD
        in_orbits = Orbit.from_table(self.table)
        in_orbits['epoch'] = in_orbits['epoch'].tt

        oo_orbits, err = pyoorb.pyoorb.oorb_element_transformation(
            in_orbits=in_orbits._to_oo(),
            in_element_type={'CART': 1, 'COM': 2, 'KEP': 3,
                             'DEL': 4, 'EQX': 5}[orbittype])

        if err != 0:
            OpenOrbError('pyoorb failed with error code {:d}'.format(err))

        orbits = Orbit._from_oo(oo_orbits, orbittype, timescale)

        # replace id column with actual target names from original orbits
        orbits.table.replace_column('id', self['targetname'])

        # replace orbtype column
        orbits.table.replace_column('orbtype', [orbittype] * len(orbits.table))

        return orbits

    @requires("pyoorb")
    @cite({'method': '2009M&PS...44.1853G',
           'software': 'https://github.com/oorb/oorb'})
    def oo_propagate(self, epochs, dynmodel='N', ephfile='de430'):
        """Uses pyoorb to propagate this `~Orbit` object. Required fields are:

        * target identifier (``'targetname'``)
        * semi-major axis (``'a'``, for Keplerian orbit) or perihelion
          distance (``'q'``, for cometary orbit), typically in au or
          or x-component of state vector (``'x'``, for cartesian orbit),
          typically in au
        * eccentricity (``'e'``, for Keplerian or cometary orbit) or
          y-component of state vector (``'y'``, for cartesian orbit) in
          au
        * inclination (``'i'``, for Keplerian or cometary orbit) in
          degrees or z-component of state vector (``'z'``, for cartesian
          orbit) in au
        * longitude of the ascending node (``'Omega'``, for Keplerian or
          cometary orbit) in degrees or x-component of velocity vector
          (``'vx'``, for cartesian orbit), au/day
        * argument of the periapsis (``'w'``, for Keplerian or cometary
          orbit) in degrees or y-component of velocity vector (``'vy'``,
          for cartesian orbit) in au/day
        * mean anomaly (``'M'``, for Keplerian orbits) in degrees or
          perihelion epoch (``'Tp'``, for cometary orbits) as astropy Time
          object or z-component of velocity vector (``'vz'``, for cartesian
          orbit) in au/day
        * epoch (``'epoch'``) as `~astropy.time.Time` object
        * absolute magnitude (``'H'``) in mag
        * photometric phase slope (``'G'``)

        Parameters
        ----------
        epochs : `~astropy.time.Time` object
            Epoch to which the orbit will be propagated to. Must be a
            `~astropy.time.Time` object holding a single epoch or
            multiple epochs (see :ref:`epochs`).
            The resulting `~sbpy.data.Orbit` object will
            have the same time scale as ``epochs``.
        dynmodel : str, optional
            The dynamical model to be used in the propagation: ``'N'``
            for n-body simulation or ``'2'`` for a 2-body
            simulation. Default: ``'N'``
        ephfile : str, optional
            Planet and Lunar ephemeris file version as provided by JPL
            to be used in the propagation. Default: ``'de430'``

        Returns
        -------
        `~Orbit` object

        Examples
        --------
        Propagate the orbit of Ceres 100 days into the future:

        >>> from sbpy.data import Orbit
        >>> from astropy.time import Time
        >>> epoch = Time(Time.now().jd + 100, format='jd')
        >>> ceres = Orbit.from_horizons('Ceres')      # doctest: +REMOTE_DATA
        >>> future_ceres = ceres.oo_propagate(epoch)  # doctest: +SKIP
        >>> print(future_ceres)  # doctest: +SKIP
        <QTable length=1>
           id           a                  e          ...    H       G    timescale
                        AU                            ...   mag
          str7       float64            float64       ... float64 float64    str3
        ------- ----------------- ------------------- ... ------- ------- ---------
        1 Ceres 2.769331727251861 0.07605371361208543 ...    3.34    0.12       UTC

        """

        if pyoorb is None:
            raise RequiredPackageUnavailable('pyoorb')

        # initialize pyoorb
        if os.getenv('OORB_DATA') is None:
            # oorb installed using conda
            pyoorb.pyoorb.oorb_init()
        else:
            ephfile = os.path.join(os.getenv('OORB_DATA'), ephfile+'.dat')
            pyoorb.pyoorb.oorb_init(ephfile)

        # extract time scale
        timescale = epochs.scale.upper()

        # identify orbit type based on available table columns
        orbittype = None
        for testtype in ['KEP', 'COM', 'CART']:
            field_names = [
                field[0] for field in Conf.oorb_orbit_fields[testtype]
            ]
            try:
                self._translate_columns(field_names[1:6])
                orbittype = testtype
                break
            except KeyError:
                pass

        if orbittype is None:
            raise OrbitError(
                'orbit type cannot be determined from elements')

        # derive and apply default units
        default_units = {}
        for field_name, field_unit in Conf.oorb_orbit_fields[orbittype]:
            try:
                # use the primary field name
                primary_field_name = self._translate_columns(field_name)[0]
            except KeyError:
                continue

            default_units[primary_field_name] = field_unit

        for colname in self.field_names:
            if (colname in default_units.keys() and
                not isinstance(self[colname],
                               (u.Quantity, u.CompositeUnit))):
                self[colname].unit = default_units[colname]

        # epochs may be a scalar value, but we need an interable
        mjd = epochs.tt.mjd.reshape(epochs.size)
        ooepoch = np.array(
            list(itertools.zip_longest(mjd, [3])), dtype=np.double, order="F"
        )

        # convert epochs to TT and MJD
        in_orbits = Orbit.from_table(self.table)
        in_orbits['epoch'] = in_orbits['epoch'].tt

        oo_orbits, err = pyoorb.pyoorb.oorb_propagation(
            in_orbits=in_orbits._to_oo(),
            in_epoch=ooepoch,
            in_dynmodel=dynmodel)

        if err != 0:
            OpenOrbError('pyoorb failed with error code {:d}'.format(err))

        orbits = Orbit._from_oo_propagatation(oo_orbits, orbittype, timescale)

        # replace id column with actual target names from original orbits
        orbits.table.replace_column('id', self.table['targetname'])

        orbits.meta['orbit_type'] = orbittype
        orbits.table.remove_column('epoch_scale'),

        return orbits

    @cite({'method': '1997Icar..127...13L'})
    def tisserand(self, planet='599', epoch=None):
        """Tisserand parameter with respect to a planet


        Parameters
        ----------
        planet : str, `~Orbit` object, or sequence of str, optional
            Planet(s) against which the Tisserand parameter is calculated.
            If `str`, then the orbital elements of the planet will be
            automaticallyl pulled from JPL Horizons. If `self` and/or
            `planet` contains more than one object, then `numpy` broadcasting
            rules apply.  Default is Jupiter.
        epoch : `~astropy.time.Time`, optional
            The epoch of planet orbit if pulled from JPL Horizons.  This
            parameter will be passed to `~sbpy.data.Orbit.from_horizons`.
            If the planet orbit is passed directly, then this parameter
            has no effect.

        Returns
        -------
        u.Quantity
            The Tisserand parameter(s)


        References
        ----------
        Levison, H. F., Duncan, M. J. 1997, Icarus 127, 13


        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> comets = Orbit.from_horizons(['252P', 'P/2016 BA14'],
        ...     id_type='designation', closest_apparition=True)
        ...     # doctest: +REMOTE_DATA
        >>> T_J = comets.tisserand() # doctest: +REMOTE_DATA
        >>>
        >>> halley = Orbit.from_horizons('1P', id_type='designation',
        ...     closest_apparition=True) # doctest: +REMOTE_DATA
        >>> T = halley.tisserand(['599', '699', '799', '899']) # doctest: +REMOTE_DATA
        """

        if isinstance(planet, str) \
            or (hasattr(planet, '__iter__')
                and np.all([isinstance(x, str) for x in planet])):
            planet = Orbit.from_horizons(planet, id_type=None, epochs=epoch)

        a_p = planet['a']
        t = a_p / self['a'] + 2 * np.cos(self['i']) * \
            np.sqrt((1 - self['e']**2) * self['a'] / a_p)
        t = u.Quantity(t)
        if len(t) == 1:
            t = t[0]
        return t

    def D_criterion(self, obj, version='sh'):
        """Evaluate orbit similarity D-criterion

        Three different versions of D-criterion are defined and often compared
        to each other, includingthe Southworth-Hawkins function [SH63]_,
        Drummond function [D81]_, and the hybrid version [J93]_.  See review by
        [W19]_.

        Parameters
        ----------
        obj : `~Orbit` object
            Object(s) against which to calculate D-criterion
        version : ['sh', 'd', 'h'], optional
            Select the versions of D-criterion formula.  Case insensitive.
            'sh' : Southworth-Hawkins function
            'd' : Drummond function
            'h' : Hybrid function
            See references for the details of each version.


        Returns
        -------
        float or numpy.ndarray


        References
        ----------
        .. [SH63] `Southwarth, R. B., & Hawkins, G. S. 1963, SCoA, 7, 261
           <https://ui.adsabs.harvard.edu/abs/1963SCoA....7..261S/abstract>`_

        .. [D81] `Drummond, J. D. 1981, Icarus 45, 545
           <https://ui.adsabs.harvard.edu/abs/1981Icar...45..545D/abstract>`_

        .. [J93] `Jopek, T. J. 1993, Icarus 106, 603
           <https://ui.adsabs.harvard.edu/abs/1993Icar..106..603J/abstract>`_

        .. [W19] `Williams, I. P., Jopek, T. J., Rudawska, R., Tóth, J.,
           Kornoš, L. 2019, In: Ryabova, G. O., Asher, D. J., Compbell-Brown
           M. D., eds.), Cambridge, UK: Cambridge University Press, 210-234.
           <https://ui.adsabs.harvard.edu/abs/2019msme.book..210W/abstract>`_


        Examples
        --------
        >>> from sbpy.data import Orbit
        >>> comets = Orbit.from_horizons(['252P', 'P/2016 BA14'],
        ...     id_type='designation', closest_apparition=True)
        ...     # doctest: +REMOTE_DATA
        >>> # Southworth & Hawkins function
        >>> D_SH = comets[0].D_criterion(comets[1]) # doctest: +REMOTE_DATA
        >>> # Drummond function
        >>> D_D = comets[0].D_criterion(comets[1], version='d') # doctest: +REMOTE_DATA
        >>> # hybrid function
        >>> D_H = comets[0].D_criterion(comets[1], version='h') # doctest: +REMOTE_DATA
        """

        if version.lower() not in ['sh', 'd', 'h']:
            raise ValueError("version should be one of ['sh', 'd', 'h'] (case "
                             "insensitive, {} received".format(version))

        diff_e = obj['e'] - self['e']
        sum_e = obj['e'] + self['e']
        diff_q = obj['q'].to_value(u.au) - self['q'].to_value(u.au)
        sum_q = obj['q'].to_value(u.au) + self['q'].to_value(u.au)
        diff_i = obj['i'] - self['i']
        sum_i = obj['i'] + self['i']
        diff_omega = obj['Omega'] - self['Omega']
        diff_w = obj['w'] - self['w']

        # sin(I_AB / 2)**2
        sin_i2 = np.sin(diff_i / 2)**2 \
            + np.sin(self['i']) * np.sin(obj['i']) * np.sin(diff_omega / 2)**2

        if version.lower() == 'd':
            # Drummond function
            register(self.D_criterion, {'method': '1981Icar...45..545D'})
            i_ba = np.arcsin(np.sqrt(sin_i2)) * 2
            beta = [np.arcsin(np.sin(o['i']) * np.sin(o['w']))
                    for o in [obj, self]]
            lamb = [o['Omega'] + np.arctan(np.cos(o['i']) * np.tan(o['w']))
                    + (np.cos(o['w']) < 0).astype(int) * np.pi * u.rad
                    for o in [obj, self]]
            theta_ba = np.arccos(np.sin(beta[0]) * np.sin(beta[1])
                                 + np.cos(beta[0]) * np.cos(beta[1])
                                 * np.cos(lamb[1] - lamb[0]))
            d2 = (diff_e / sum_e)**2 + (diff_q / sum_q)**2 \
                + (i_ba / (np.pi * u.rad))**2 \
                + (sum_e * theta_ba / (2 * np.pi * u.rad))**2
        else:
            cos_i2 = np.sqrt(1 - sin_i2)
            sign = (abs(diff_omega) <= 180 * u.deg).astype(int) * 2 - 1
            pi_ba = diff_w + 2 * sign * np.arcsin(
                np.cos(sum_i / 2) * np.sin(diff_omega / 2) / cos_i2)

            if version.lower() == 'sh':
                # Southworth-Hawkins function
                register(self.D_criterion, {'method': '1963SCoA....7..261S'})
                d2 = diff_e**2 + diff_q**2 + 4 * sin_i2 \
                    + (sum_e * np.sin(pi_ba / 2))**2
            else:
                # hybrid function
                register(self.D_criterion, {'method': '1993Icar..106..603J'})
                d2 = diff_e**2 + (diff_q / sum_q)**2 + 4 * sin_i2 \
                    + (sum_e * np.sin(pi_ba / 2))**2
        d = np.sqrt(u.Quantity(d2))
        if len(d) == 1:
            d = d[0]
        return d
