# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Ephem Module
======================

Class for storing and querying ephemerides

created on June 04, 2017
"""
import os
from warnings import warn

from numpy import ndarray, hstack, iterable
from astropy.time import Time
from astropy.table import vstack, Column
import astropy.units as u
from astroquery.jplhorizons import Horizons
from astroquery.mpc import MPC
from astroquery.imcce import Miriade
from astroquery.exceptions import InvalidQueryError
from astropy.coordinates import EarthLocation


from ..bib import cite
from .core import DataClass, conf, QueryError, TimeScaleWarning
from ..exceptions import SbpyException
from .orbit import OpenOrbError

__all__ = ['Ephem']


class Ephem(DataClass):
    """Class for querying, manipulating, and calculating ephemerides"""

    @classmethod
    @cite({'data source': '1996DPS....28.2504G'})
    @cite({'software: astroquery': '2019AJ....157...98G'})
    def from_horizons(cls, targetids, id_type='smallbody',
                      epochs=None, location='500', **kwargs):
        """Load target ephemerides from
        `JPL Horizons <https://ssd.jpl.nasa.gov/horizons.cgi>`_ using
        `astroquery.jplhorizons.HorizonsClass.ephemerides`

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
        epochs : `~astropy.time.Time` object, or dictionary, optional
            Epochs of elements to be queried; `~astropy.time.Time` objects
            support single and multiple epochs;
            a dictionary including keywords ``start`` and ``stop``, as well
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
            used. All epochs should be provided in UTC; if not, they will be
            converted to UTC and a `~sbpy.data.TimeScaleWarning` will be
            raised.  Default: ``None``
        location : str or `~astropy.coordinates.EarthLocation`, optional
            Location of the observer using IAU observatory codes
            (see `IAU observatory codes
            <https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html>`__)
            or as `~astropy.coordinates.EarthLocation`.
            Default: ``'500'`` (geocentric)
        **kwargs : optional
            Arguments that will be provided to
            `astroquery.jplhorizons.HorizonsClass.ephemerides`.

        Notes
        -----
        * For detailed explanations of the queried fields, refer to
          `astroquery.jplhorizons.HorizonsClass.ephemerides` and the
          `JPL Horizons documentation <https://ssd.jpl.nasa.gov/?horizons_doc>`_.
        * By default, all properties are provided in the J2000.0 reference
          system. Different settings can be chosen using
          additional keyword arguments as used by
          `astroquery.jplhorizons.HorizonsClass.ephemerides`.


        Returns
        -------
        `~Ephem` object
            The resulting object will be populated with columns as
            defined in
            `~astroquery.jplhorizons.HorizonsClass.ephemerides`; refer
            to that document on information on how to modify the list
            of queried parameters.

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_horizons('ceres', epochs=epoch) # doctest: +SKIP

        """

        # modify epoch input to make it work with astroquery.jplhorizons
        # maybe this stuff should really go into that module....
        if epochs is None:
            epochs = [Time.now().utc.jd]
        elif isinstance(epochs, Time):
            if epochs.scale is not 'utc':
                warn(('converting {} epochs to utc for use in '
                      'astroquery.jplhorizons').format(epochs.scale),
                     TimeScaleWarning)
            epochs = epochs.utc.jd
        elif isinstance(epochs, dict):
            if 'start' in epochs and 'stop' in epochs and 'number' in epochs:
                epochs['step'] = epochs['number']*u.dimensionless_unscaled
            # convert to utc and iso for astroquery.jplhorizons
            epochs['start'] = epochs['start'].utc.iso
            epochs['stop'] = epochs['stop'].utc.iso
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

        # turn EarthLocation into dictionary of strings as used by
        # astroquery.jplhorizons
        if isinstance(location, EarthLocation):
            location = {'lon': location.lon.deg,
                        'lat': location.lat.deg,
                        'elevation': location.height.to('km')}

        # append ephemerides table for each targetid
        all_eph = None
        for targetid in targetids:

            # load ephemerides using astroquery.jplhorizons
            obj = Horizons(id=targetid, id_type=id_type,
                           location=location, epochs=epochs)
            try:
                eph = obj.ephemerides(**kwargs)
            except ValueError as e:
                raise QueryError(
                    ('Error raised by astroquery.jplhorizons: {:s}\n'
                     'The following query was attempted: {:s}').format(
                         str(e), obj.uri))

            # workaround for current version of astroquery to make
            # column units compatible with astropy.table.QTable
            # should really change '---' units to None in
            # astroquery.jplhorizons.__init__.py
            for column_name in eph.columns:
                if eph[column_name].unit == '---':
                    eph[column_name].unit = None

            # workaround for astroquery 0.3.9.dev5056 and earlier,
            # Horizons column named RA_rate always includes the
            # cos(Dec) term:
            if 'RA_rate' in eph.colnames:
                eph['RA_rate'].name = 'RA*cos(Dec)_rate'

            if all_eph is None:
                all_eph = eph
            else:
                all_eph = vstack([all_eph, eph])

        # turn epochs into astropy.time.Time and apply timescale
        # convert ut1 epochs to utc
        # https://ssd.jpl.nasa.gov/?horizons_doc
        if any(all_eph['datetime_jd'] < 2437665.5):
            all_eph['datetime_jd'][all_eph['datetime_jd'] <
                                   2437665.5] = Time(
                all_eph['datetime_jd'][all_eph['datetime_jd'] <
                                       2437665.5],
                                       scale='ut1', format='jd').utc.jd
        all_eph['epoch'] = Time(all_eph['datetime_jd'], format='jd',
                                scale='utc')
        all_eph['siderealtime'].unit = u.Unit('hour')

        all_eph.remove_column('datetime_jd')
        all_eph.remove_column('datetime_str')

        return cls.from_table(all_eph)

    @classmethod
    @cite({'data source':
           'https://minorplanetcenter.net/iau/MPEph/MPEph.html'})
    @cite({'software: astroquery': '2019AJ....157...98G'})
    def from_mpc(cls, targetids, epochs=None, location='500', **kwargs):
        """Load ephemerides from the
        `Minor Planet Center <http://minorplanetcenter.net>`_.

        Parameters
        ----------
        targetids : str or iterable of str
            Target identifier, resolvable by the Minor Planet
            Ephemeris Service [MPES]_, e.g., 2P, C/1995 O1, P/Encke,
            (1), 3200, Ceres, and packed designations, for one or more
            targets.

        epochs : `~astropy.time.Time` object, or dictionary, optional
            Request ephemerides at these epochs.  May be a single
            epoch or multiple epochs as `~astropy.time.Time`
            (see :ref:`epochs`)or a
            dictionary describing a linearly-spaced array
            of epochs. All epochs should be provided in UTC; if not,
            they will be converted to UTC and a
            `~sbpy.data.TimeScaleWarning` will be raised.
            If ``None`` (default), the current date and
            time will be used.

            For the dictionary format, the keys ``start`` (start
            epoch), ``step`` (step size), ``stop`` (end epoch), and/or
            ``number`` (number of epochs total) are used.  Only one of
            ``stop`` and ``number`` may be specified at a time.
            ``step``, ``stop``, and ``number`` are optional. The values of
            of ``start`` and ``stop`` must be `~astropy.time.Time` objects.
            ``number`` should be an integer value; ``step`` should be a
            `~astropy.units.Quantity` with an integer value and units of
            seconds, minutes, hours, or days.

            All epochs should be provided in UTC; if not, they will be
            converted to UTC and a `~sbpy.data.TimeScaleWarning` will be
            raised.

        location : various, optional
            Location of the observer as an IAU observatory code
            [OBSCODES]_ (string), a 3-element array of Earth
            longitude, latitude, altitude, or an
            `~astropy.coordinates.EarthLocation`.  Longitude and
            latitude should be parseable by
            `~astropy.coordinates.Angle`, and altitude should be
            parsable by `~astropy.units.Quantity` (with units of
            length).  If ``None``, then the geocenter (code 500) is
            used.

        **kwargs
            Additional keyword arguments are passed to
            `~astroquery.mpc.MPC.get_ephemerides`: ``eph_type``,
            ``ra_format``, ``dec_format``, ``proper_motion``,
            ``proper_motion_unit``, ``suppress_daytime``,
            ``suppress_set``, ``perturbed``, ``unc_links``, ``cache``.

        Returns
        -------
        `~Ephem` object
            The resulting object will be populated with columns as
            defined in
            `~astroquery.mpc.get_ephemerides`; refer
            to that document on information on how to modify the list
            of queried parameters.


        Examples
        --------
        Query a single set of ephemerides of Ceres as observed from
        Maunakea:
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_mpc('ceres', epoch, 568) # doctest: +REMOTE_DATA

        Query a range of ephemerides of comet 2P/Encke as observed from
        Maunakea:
        >>> epochs = {'start': Time('2019-01-01'),
        ...           'step': 1*u.d, 'number': 365}
        >>> eph = Ephem.from_mpc('2P', epochs, 568) # doctest: +REMOTE_DATA

        Notes
        -----
        * All properties are provided in the J2000.0 reference system.
        * See `astroquery.mpc.MPC.get_ephemerides` and the Minor
          Planet Ephemeris Service user's guide [MPES]_ for details,
          including accetable target names.


        References
        ----------
        .. [MPES] Wiliams, G. The Minor Planet Ephemeris Service.
           https://minorplanetcenter.org/iau/info/MPES.pdf

        .. [OBSCODES] IAU Minor Planet Center.  List of observatory
           codes. https://minorplanetcenter.org/iau/lists/ObsCodesF.html

        """

        # parameter check

        # if targetids is a list, run separate Horizons queries and append
        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        if isinstance(epochs, Time):
            if epochs.scale is not 'utc':
                warn(('converting {} epochs to utc for use in '
                      'astroquery.mpc').format(epochs.scale),
                     TimeScaleWarning)
                epochs = epochs.utc
            start = None
        elif isinstance(epochs, dict):
            start = epochs['start']  # required
            if start.scale is not 'utc':
                warn(('converting {} start epoch to utc for use in '
                      'astroquery.mpc').format(start.scale),
                     TimeScaleWarning)
                start = start.utc
            step = epochs.get('step')
            stop = epochs.get('stop')
            if stop is not None and stop.scale is not 'utc':
                warn(('converting {} stop epoch to utc for use in '
                      'astroquery.mpc').format(stop.scale),
                     TimeScaleWarning)
                stop = stop.utc
            number = epochs.get('number')

            if step is not None and stop is None:
                step = u.Quantity(step)
                if step.unit not in (u.d, u.h, u.min, u.s):
                    raise QueryError(
                        'step must have units of days, hours, minutes,'
                        ' or seconds')
            if stop is not None:
                if step is not None and number is None:
                    # start and stop both defined, estimate number of steps
                    dt = (Time(stop).jd - Time(start).jd) * u.d
                    number = int((dt / step).decompose()) + 1
                elif step is None and number is not None:
                    step = int((stop-start).jd*1440/(number-1))*u.minute
                else:
                    raise QueryError(
                        ('epoch definition unclear; step xor number '
                         'must be provided with start and stop'))
        else:
            start = None

        if epochs is None:
            epochs = Time.now()

        if not iterable(epochs):
            epochs = [epochs]

        # append ephemerides table for each targetid
        all_eph = None
        for targetid in targetids:

            try:
                # get ephemeris
                if start is None:
                    eph = []
                    for i in range(len(epochs)):
                        e = MPC.get_ephemeris(targetid, location=location,
                                              start=Time(epochs[i],
                                                         scale='utc'),
                                              number=1, **kwargs)
                        e['Date'] = e['Date'].iso  # for vstack to work
                        eph.append(e)
                    eph = vstack(eph)
                    eph['Date'] = Time(eph['Date'], scale='utc')
                else:
                    eph = MPC.get_ephemeris(targetid, location=location,
                                            start=start, step=step,
                                            number=number, **kwargs)
            except InvalidQueryError as e:
                raise QueryError(
                    'Error raised by astroquery.mpc: {:s}'.format(str(e)))

            # add targetname column
            eph.add_column(Column([targetid]*len(eph),
                                  name='Targetname'), index=0)

            if all_eph is None:
                all_eph = eph
            else:
                all_eph = vstack([all_eph, eph])

        # if ra_format or dec_format is defined, then units must be
        # dropped or else QTable will raise an exception because
        # strings cannot have units
        if 'ra_format' in kwargs:
            all_eph['RA'].unit = None
        if 'dec_format' in kwargs:
            all_eph['Dec'].unit = None

        return cls.from_table(all_eph)

    @classmethod
    @cite({'data source': 'http://vo.imcce.fr/webservices/miriade/'})
    @cite({'software: astroquery': '2019AJ....157...98G'})
    def from_miriade(cls, targetids, objtype='asteroid',
                     epochs=None, location='500', **kwargs):
        """Load target ephemerides from
        `IMCCE Miriade <http://vo.imcce.fr/webservices/miriade/>`_ using
        `astroquery.imcce.MiriadeClass.get_ephemerides`

        Parameters
        ----------
        targetids : str or iterable of str
            Target identifier, i.e., a number, name, designation, or JPL
            Horizons record number, for one or more targets.
        objtype : str, optional
            The nature of ``targetids`` provided; possible values are
            ``'asteroid'``, ``'comet'``, ``'dwarf planet'``,
            ``'planet'``, or ``'satellite'``. Default: ``'asteroid'``
        epochs : `~astropy.time.Time` object, or dictionary, optional
            Epochs of elements to be queried; `~astropy.time.Time` objects
            support single and multiple epochs;
            a dictionary including keywords ``start`` and ``stop``, as well
            as either ``step`` or ``number``, can be used to generate a range
            of epochs. ``start`` and ``stop`` have to be
            `~astropy.time.Time` objects (see :ref:`epochs`).
            If ``step`` is provided, a range
            of epochs will be queried starting at ``start`` and ending at
            ``stop`` in steps of ``step``; ``step`` has to be provided as
            a `~astropy.units.Quantity` object with integer value and a
            unit of either seconds, minutes, hours, or days. If
            ``number`` is provided as an integer, the
            interval defined by
            ``start`` and ``stop`` is split into ``number`` equidistant
            intervals. If ``None`` is
            provided, current date and time are
            used. All epochs should be provided in UTC; if not, they will be
            converted to UTC and a `~sbpy.data.TimeScaleWarning` will be
            raised. Default: ``None``
        location : str or `~astropy.coordinates.EarthLocation`, optional
            Location of the observer using IAU observatory codes
            (see `IAU observatory codes
            <https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html>`__)
            or as `~astropy.coordinates.EarthLocation`.
            Default: ``'500'`` (geocentric)
        **kwargs : optional
            Arguments that will be provided to
            `astroquery.imcce.MiriadeClass.get_ephemerides`.

        Notes
        -----
        * For detailed explanations of the queried fields, refer to
          `astroquery.imcce.MiriadeClass.get_ephemerides` and the
          `Miriade documentation
          <http://vo.imcce.fr/webservices/miriade/?documentation>`_.
        * By default, all properties are provided in the J2000.0 reference
          system. Different settings can be chosen using
          additional keyword arguments as used by
          `astroquery.imcce.MiriadeClass.get_ephemerides`.

        Returns
        -------
        `~Ephem` object
            The resulting object will be populated with columns as
            defined in
            `~astroquery.imcce.MiriadeClass.get_ephemerides`; refer
            to that document on information on how to modify the list
            of queried parameters.

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_horizons('ceres', epochs=epoch) # doctest: +SKIP
        """

        # modify epoch input to make it work with astroquery.imcce.Miriade
        if epochs is None:
            epochs = {'start': Time.now().utc.jd}
        elif isinstance(epochs, Time):
            if epochs.scale is not 'utc':
                warn(('converting {} epochs to utc for use in '
                      'astroquery.imcce').format(epochs.scale),
                     TimeScaleWarning)
                epochs = epochs.utc
            epochs = {'start': epochs}
        elif isinstance(epochs, dict):
            if epochs['start'].scale is not 'utc':
                warn(('converting {} start epoch to utc for use in '
                      'astroquery.imcce').format(epochs['start'].scale),
                     TimeScaleWarning)
                epochs['start'] = epochs['start'].utc
            if 'stop' in epochs and epochs['stop'].scale is not 'utc':
                warn(('converting {} stop epoch to utc for use in '
                      'astroquery.imcce').format(epochs['stop'].scale),
                     TimeScaleWarning)
                epochs['stop'] = epochs['stop'].utc
            if 'number' in epochs:
                # turn interval/number into step size based on full minutes
                epochs['step'] = int((Time(epochs['stop']) -
                                      Time(epochs['start'])).jd *
                                     86400/(epochs['number']-1))*u.s
            elif 'step' in epochs:
                epochs['number'] = ((Time(epochs['stop']) -
                                     Time(epochs['start'])).jd *
                                    86400/epochs['step'].to('s').value)+1
            if 'step' in epochs:
                epochs['step'] = '{:f}{:s}'.format(
                    epochs['step'].value,
                    {u.s: 's', u.minute: 'm', u.hour: 'h',
                     u.day: 'd'}[epochs['step'].unit])

        # if targetids is a list, run separate Horizons queries and append
        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        # turn EarthLocation into dictionary of strings as used by
        # astroquery.jplhorizons

        if isinstance(location, EarthLocation):
            location = '{:+f} {:+f} {:.1f}'.format(
                location.lon.deg,
                location.lat.deg,
                location.height.to('m').value)

        # append ephemerides table for each targetid
        all_eph = None
        for targetid in targetids:
            query = Miriade()
            try:
                if 'step' not in epochs and 'number' not in epochs:
                    if not iterable(epochs['start']):
                        # single epoch
                        eph = query.get_ephemerides(targetname=targetid,
                                                    objtype=objtype,
                                                    location=location,
                                                    epoch=epochs['start'],
                                                    **kwargs)
                    else:
                        # multiple epochs
                        eph = []
                        for i in range(len(epochs['start'])):
                            e = query.get_ephemerides(targetname=targetid,
                                                      objtype=objtype,
                                                      location=location,
                                                      epoch=epochs['start'][i],
                                                      **kwargs)
                            e['epoch'] = Time(e['epoch'], format='jd',
                                              scale='utc').iso
                            eph.append(e)
                        eph = vstack(eph)
                        eph['epoch'] = Time(eph['epoch'], scale='utc',
                                            format='iso')
                else:
                    # dictionary
                    eph = query.get_ephemerides(
                        targetname=targetid, objtype=objtype,
                        location=location, epoch=epochs['start'],
                        epoch_step=epochs['step'],
                        epoch_nsteps=epochs['number'],
                        **kwargs)
            except RuntimeError as e:
                raise QueryError(
                    ('Error raised by astroquery.imcce: {:s}\n'
                     'The following query was attempted: {:s}').format(
                         str(e), query.uri))

            if all_eph is None:
                all_eph = eph
            else:
                all_eph = vstack([all_eph, eph])

        self = cls.from_table(all_eph)

        # turn epochs into astropy.time.Time and apply timescale
        self.table['epoch'] = Time(self.table['epoch'],
                                   format='jd', scale='utc')

        return self

    @classmethod
    @cite({'method': '2009M&PS...44.1853G',
           'software': 'https://github.com/oorb/oorb'})
    def from_oo(self, orbit, epochs=None, location='500', scope='full',
                dynmodel='N', ephfile='de430'):
        """Uses pyoorb to derive ephemerides from an `~Orbit` object. For a
        list of output parameters, please read the `pyoorb documentation
        <https://github.com/oorb/oorb/tree/master/python>`_.

        Parameters
        ----------
        orbit : `~Orbit` object
            Can contain any number of orbits, ephemerides will be calculated
            for each orbit. Required fields are:

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
            * epoch (``'epoch'``) as `~astropy.time.Time`
            * absolute magnitude (``'H'``) in mag
            * photometric phase slope (``'G'``)

        epochs : `~astropy.time.Time` object, optional
            Epochs of elements to be queried; must be a
            `~astropy.time.Time` object holding a single or multiple epochs
            (see :ref:`epochs`).
            If ``None`` is provided, current date and
            time are used. The same time scale that is used in ``epochs``
            will be applied to the results. Default: ``None``
        location : str, optional, default ``'500'`` (geocentric)
            Location of the observer.
        scope : str
            Scope of data to be determined: ``'full'`` obtains all
            available properties, ``'basic'`` obtains only a limited
            amount of data. Default: ``'full'``
        dynmodel : str, optional
            The dynamical model to be used in the propagation: ``'N'``
            for n-body simulation or ``'2'`` for a 2-body
            simulation. Default: ``'N'``
        ephfile : str, optional
            Planet and Lunar ephemeris file version as provided by JPL
            to be used in the propagation. Default: ``'de430'``

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        Compute ephemerides for Ceres as seen from the Discovery Channel
        Telescope for the next 10 days at 1hr intervals:

        >>> import numpy as np
        >>> from sbpy.data import Orbit, Ephem
        >>> from astropy.time import Time
        >>> epochs = Time(Time.now().jd + np.arange(0, 10, 1/24), format='jd')
        >>> ceres = Orbit.from_horizons('1')    # doctest: +REMOTE_DATA
        >>> eph = Ephem.from_oo(ceres, epochs, 'G37') # doctest: +REMOTE_DATA
        >>> eph  # doctest: +SKIP
        <QTable length=240>
        targetname       epoch        ...           obsz               trueanom
                           d          ...            AU                  deg
           str7         float64       ...         float64              float64
        ---------- ------------------ ... ----------------------- -----------------
           1 Ceres  2458519.316966272 ...  3.2083678848104924e-06  68.0863831954328
           1 Ceres 2458519.3586329385 ...  2.7022422510736277e-07 68.09589266358881
           1 Ceres 2458519.4002996054 ...  -3.111046209036683e-06 68.10540191585879
           1 Ceres  2458519.441966272 ...  -6.700369254264427e-06 68.11491095202307
           1 Ceres 2458519.4836329385 ... -1.0248419404668141e-05 68.12441977218093
           1 Ceres 2458519.5252996054 ... -1.3508703580356052e-05 68.13392837643161
               ...                ... ...                     ...               ...
           1 Ceres  2458529.066966272 ...  1.2522500440509399e-05 70.30569661787204
           1 Ceres 2458529.1086329385 ...  1.4101698473351076e-05 70.31515536712485
           1 Ceres 2458529.1502996054 ...  1.4771304981564537e-05  70.3246138990413
           1 Ceres  2458529.191966272 ...   1.448582020449618e-05 70.33407221340468
           1 Ceres 2458529.2336329385 ...   1.326517587380005e-05 70.34353031031534
           1 Ceres 2458529.2752996054 ...  1.1193369555934085e-05 70.35298818987367        """

        import pyoorb

        # create a copy of orbit
        from . import Orbit
        orb = Orbit.from_table(orbit.table)

        if epochs is None:
            epochs = Time.now()

        # extract time scale
        timescale = epochs.scale.upper()

        # initialize pyoorb
        if os.getenv('OORB_DATA') is None:
            # oorb installed using conda
            pyoorb.pyoorb.oorb_init()
        else:
            ephfile = os.path.join(os.getenv('OORB_DATA'), ephfile+'.dat')
            pyoorb.pyoorb.oorb_init(ephfile)

        # identify orbit type based on available table columns
        orbittype = None
        for testtype in ['KEP', 'COM', 'CART']:
            try:
                orb._translate_columns(
                    conf.oorb_orbit_fields[testtype][1:6])
                orbittype = testtype
                break
            except KeyError:
                pass

        if orbittype is None:
            raise OpenOrbError(
                'orbit type cannot be determined from elements')

        # add/update orbittype column
        orb['orbittype'] = [orbittype] * len(orb)

        # derive and apply default units
        default_units = {}
        for idx, field in enumerate(conf.oorb_orbit_fields[orbittype]):
            try:
                default_units[orb._translate_columns(
                    field)[0]] = conf.oorb_orbit_units[orbittype][idx]
            except KeyError:
                pass
        for colname in orb.field_names:
            if (colname in default_units.keys() and
                not isinstance(orb[colname],
                               (u.Quantity, u.CompositeUnit, Time))):
                orb[colname].unit = default_units[colname]

        # convert epochs to TT
        orb['epoch'] = orb['epoch'].tt
        epochs = epochs.tt

        try:
            epochs = list(zip(epochs.mjd,
                              [conf.oorb_timeScales['TT']]*len(epochs)))
        except TypeError:
            epochs = [(epochs.mjd, conf.oorb_timeScales['TT'])]

        if scope == 'full':
            oo_eph, err = pyoorb.pyoorb.oorb_ephemeris_full(
                orb._to_oo(),
                location,
                epochs,
                dynmodel)
        elif scope == 'basic':
            oo_eph, err = pyoorb.pyoorb.oorb_ephemeris_basic(
                orb._to_oo(),
                location,
                epochs,
                dynmodel)

        if err != 0:
            OpenOrbError('pyoorb failed with error code {:d}'.format(err))

        # reorder data on per-column basis and apply units
        oo_eph_col = hstack([oo_eph.transpose()[:, :, i]
                             for i in range(oo_eph.shape[0])]).tolist()
        oo_eph_col_u = []
        if scope == 'full':
            for i, col in enumerate(oo_eph_col):
                oo_eph_col_u.append(Ephem._unit_apply(
                    col, conf.oorb_ephem_full_units[i]))
            ephem = self.from_columns(oo_eph_col_u,
                                      names=conf.oorb_ephem_full_fields)
        elif scope == 'basic':
            for i, col in enumerate(oo_eph_col):
                oo_eph_col_u.append(Ephem._unit_apply(
                    col, conf.oorb_ephem_basic_units[i]))
            ephem = self.from_columns(oo_eph_col_u,
                                      names=conf.oorb_ephem_basic_fields)

        # add targetname column
        ephem.table.add_column(Column(data=sum([[orb['targetname'][i]] *
                                                len(epochs) for i in
                                                range(len(orb.table))],
                                               []),
                                      name='targetname'),
                               index=0)

        # convert MJD to astropy.time.TimeJulian Date
        ephem.table['epoch'] = Time(Time(ephem['MJD'], format='mjd',
                                         scale=timescale.lower()),
                                    format='jd')
        ephem.table.remove_column('MJD')

        return ephem
