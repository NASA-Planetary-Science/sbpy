# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Orbit Module
======================

Class for querying, manipulating, integrating, and fitting orbital elements.

created on June 04, 2017
"""
import os
from numpy import array, ndarray, double, arange, rad2deg
from astropy.time import Time
from astropy.table import vstack, QTable
from astroquery.jplhorizons import Horizons
from astroquery.mpc import MPC
import astropy.units as u
from warnings import warn

from ..bib import cite
from ..exceptions import SbpyException
from . import conf, DataClass, QueryError, TimeScaleWarning

__all__ = ['Orbit', 'OrbitError', 'OpenOrbError']


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
    @cite({'data source': '1996DPS....28.2504G'})
    @cite({'software: astroquery': '2019AJ....157...98G'})
    def from_horizons(cls, targetids, id_type='smallbody',
                      epochs=None, center='500@10',
                      **kwargs):
        """Load target orbital elements from
        `JPL Horizons <https://ssd.jpl.nasa.gov/horizons.cgi>`_ using
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
          `JPL Horizons documentation <https://ssd.jpl.nasa.gov/?horizons_doc>`_.
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
        >>> epoch = Time('2018-05-14', scale='tdb')
        >>> eph = Orbit.from_horizons('Ceres', epochs=epoch)  # doctest: +REMOTE_DATA
        """

        # modify epoch input to make it work with astroquery.jplhorizons
        # maybe this stuff should really go into that module....
        if epochs is None:
            epochs = [Time.now().tdb.jd]
        elif isinstance(epochs, Time):
            if epochs.scale is not 'tdb':
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
    @cite({'data source':
           'https://minorplanetcenter.net/iau/MPEph/MPEph.html'})
    @cite({'software: astroquery': '2019AJ....157...98G'})
    def from_mpc(cls, targetids, id_type=None, target_type=None, **kwargs):
        """Load latest orbital elements from the
        `Minor Planet Center <http://minorplanetcenter.net>`_.

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
                # skip if key not in conf.mpc_orbit_fields
                if key not in conf.mpc_orbit_fields:
                    continue

                fieldname, fieldunit = conf.mpc_orbit_fields[key]
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
                try:
                    for field in conf.oorb_orbit_fields[testtype][1:6]:
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
                                  self['i'].to('radian').value,
                                  self['Omega'].to('radian').value,
                                  self['w'].to('radian').value,
                                  (self['Tp_jd'].to('d').value -
                                   2400000.5),
                                  [conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  self['epoch'].mjd,
                                  ([conf.oorb_timeScales
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
                                  self['incl'].to('radian').value,
                                  self['Omega'].to('radian').value,
                                  self['w'].to('radian').value,
                                  self['M'].to('radian').value,
                                  [conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  self['epoch'].mjd,
                                  ([conf.oorb_timeScales
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
                                  [conf.oorb_elemType[orbittype]] *
                                  len(self.table),
                                  self['epoch'].mjd,
                                  ([conf.oorb_timeScales
                                    [self['epoch'][0].scale.upper()]] *
                                   len(self.table)),
                                  self['H'].data,
                                  self['G'].data]).transpose(),
                           dtype=double, order='F')

        return orbits

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
          perihelion epoch (``'Tp_jd'``, for cometary orbits) in JD or
          z-component of velocity vector (``'vz'``, for cartesian orbit)
          in au/day
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
        import pyoorb

        # initialize pyoorb
        if os.getenv('OORB_DATA') is None:
            # oorb installed using conda
            pyoorb.pyoorb.oorb_init()
        else:
            ephfile = os.path.join(os.getenv('OORB_DATA'), ephfile+'.dat')
            pyoorb.pyoorb.oorb_init(ephfile)

        # extract time scale
        timescale = self.table['epoch'][0].scale.upper()

        # derive and apply default units
        default_units = {}
        for idx, field in enumerate(conf.oorb_orbit_fields[orbittype]):
            try:
                default_units[self._translate_columns(
                    field)[0]] = conf.oorb_orbit_units[orbittype][idx]
            except KeyError:
                pass
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

        # reorder data in Orbit object
        field_names = conf.oorb_orbit_fields[orbittype]

        columns = []
        for i, col in enumerate(oo_orbits.transpose()):
            columns.append(Orbit._unit_apply(
                col, conf.oorb_orbit_units[orbittype][i]))
        orbits = self.from_columns(columns, names=field_names)

        for i, col in enumerate(orbits.field_names):
            # convert from radians to degrees where unit == deg
            if conf.oorb_orbit_units[orbittype][i] == 'deg':
                orbits._table[col] = rad2deg(orbits[col])

        # replace id column with actual target names from original orbits
        orbits.table.replace_column('id', self['targetname'])

        # epoch column into astropy.time.Time, original time scale, and JD
        orbits.table['epoch'] = Time(
            orbits['epoch'].to('d').value,
            format='mjd', scale='tt').__getattr__(timescale.lower())

        # replace orbtype column
        orbits.table.replace_column('orbtype',
                                    [orbittype] * len(orbits.table))

        return orbits

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
          perihelion epoch (``'Tp_jd'``, for cometary orbits) in JD or
          z-component of velocity vector (``'vz'``, for cartesian orbit)
          in au/day
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
        1 Ceres 2.769331727251861 0.07605371361208543 ...    3.34    0.12       UTC        """

        import pyoorb

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
            try:
                self._translate_columns(
                    conf.oorb_orbit_fields[testtype][1:6])
                orbittype = testtype
                break
            except KeyError:
                pass

        if orbittype is None:
            raise OrbitError(
                'orbit type cannot be determined from elements')

        # derive and apply default units
        default_units = {}
        for idx, field in enumerate(conf.oorb_orbit_fields[orbittype]):
            try:
                default_units[self._translate_columns(
                    field)[0]] = conf.oorb_orbit_units[orbittype][idx]
            except KeyError:
                pass
        for colname in self.field_names:
            if (colname in default_units.keys() and
                not isinstance(self[colname],
                               (u.Quantity, u.CompositeUnit))):
                self[colname].unit = default_units[colname]

        ooepoch = [epochs.tt.mjd, conf.oorb_timeScales['TT']]

        # convert epochs to TT and MJD
        in_orbits = Orbit.from_table(self.table)
        in_orbits['epoch'] = in_orbits['epoch'].tt

        oo_orbits, err = pyoorb.pyoorb.oorb_propagation(
            in_orbits=in_orbits._to_oo(),
            in_epoch=ooepoch,
            in_dynmodel=dynmodel)

        if err != 0:
            OpenOrbError('pyoorb failed with error code {:d}'.format(err))

        # reorder data in Orbit object
        field_names = conf.oorb_orbit_fields[orbittype]

        columns = []
        for i, col in enumerate(oo_orbits.transpose()):
            columns.append(Orbit._unit_apply(
                col, conf.oorb_orbit_units[orbittype][i]))
        orbits = self.from_columns(columns, names=field_names)

        for i, col in enumerate(orbits.field_names):
            # convert from radians to degrees where unit == deg
            if conf.oorb_orbit_units[orbittype][i] == 'deg':
                orbits._table[col] = rad2deg(orbits[col])

        # replace id column with actual target names from original orbits
        orbits.table.replace_column('id', self.table['targetname'])

        orbits.meta['orbit_type'] = orbittype
        orbits.table.remove_column('epoch_scale'),

        # adjust epochs to standard jd
        orbits.table['epoch'] = Time(
            Time(orbits.table['epoch'], format='mjd',
                 scale='tt').__getattr__(timescale.lower()),
            format='jd')

        return orbits
