# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Ephem Module
======================

Class for storing and querying ephemerides

created on June 04, 2017
"""
import os
from copy import deepcopy

from numpy import ndarray, array, hstack, iterable
from astropy.time import Time
from astropy.table import vstack, Column
import astropy.units as u
from astroquery.jplhorizons import Horizons
from astroquery.mpc import MPC
from astroquery.imcce import Miriade
from astropy.coordinates import EarthLocation

from ..exceptions import SbpyException
from .. import bib
from .core import DataClass, conf

__all__ = ['Ephem']


class QueryError(SbpyException):
    pass


class Ephem(DataClass):
    """Class for querying, manipulating, and calculating ephemerides"""

    @classmethod
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
            support iterables within the object so an `~astropy.time.Time`
            object should still be used for a number of discrete epochs;
            a dictionary including keywords ``start`` and ``stop``, as well
            as either ``step`` or ``number``, can be used to generate a range
            of epochs. ``start`` and ``stop`` have to be
            `~astropy.time.Time` objects. If ``step`` is provided as a string,
            the interval
            defined by ``start``
            and ``stop`` is split into time steps defined by ``step`` (see
            `~astroquery.jplhorizons.HorizonsClass.Horizons.ephemerides`
            for details). If ``number`` is provided as an integer, the
            interval defined by
            ``start`` and ``stop`` is split into ``number`` equidistant
            intervals by calculating ``step`` from the interval and ``number``
            in units of full minutes. If ``None`` is
            provided, current date and time are
            used. Default: ``None``
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
            epochs = [Time.now().jd]
        elif isinstance(epochs, Time):
            epochs = epochs.jd
            if isinstance(epochs, float):
                epochs = [epochs]
            new_epochs = [None] * len(epochs)
            for i in range(len(epochs)):
                new_epochs[i] = epochs[i]
            epochs = new_epochs
        elif isinstance(epochs, dict):
            if 'start' in epochs and 'stop' in epochs and 'number' in epochs:
                # turn interval/number into step size based on full minutes
                epochs['step'] = '{:d}m'.format(
                    int((epochs['stop']-epochs['start']).jd*1440 /
                        (epochs['number']-1)))
            if 'start' in epochs and 'stop' in epochs and 'step' in epochs:
                # default format used in astroquery.jplhorizons
                for key, val in epochs.items():
                    if isinstance(val, Time):
                        val.format = 'iso'
                        val.out_subfmt = 'date_hms'
                        epochs[key] = val.value
                    else:
                        epochs[key] = epochs[key]

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
                raise RuntimeError(
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
        # https://ssd.jpl.nasa.gov/?horizons_doc
        all_eph['epoch'] = Time(all_eph['datetime_jd'], format='jd',
                                scale='utc')
        if any(all_eph['datetime_jd'] < 2437665.5):
            print('pimmel')
            all_eph['epoch'][all_eph['datetime_jd'] < 2437665.5] = (
                Time(all_eph['datetime_jd'][all_eph['datetime_jd'] <
                                            2437665.5],
                     format='jd', scale='ut1'))
        all_eph.remove_column('datetime_jd')

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem.from_horizons',
                         {'data service': '1996DPS....28.2504G'})

        return cls.from_table(all_eph)

    @classmethod
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

        epochs : various, optional
            Request ephemerides at these epochs.  May be a single
            epoch as a string or `~astropy.time.Time`, an array of
            epochs, or a dictionary describing a linearly-spaced array
            of epochs.  If ``None`` (default), the current date and
            time will be used.

            For the dictionary format, the keys ``start`` (start
            epoch), ``step`` (step size), ``stop`` (end epoch), and/or
            ``number`` (number of epochs total) are used.  Only one of
            ``stop`` and ``number`` may be specified at a time.
            ``step``, ``stop``, and ``number`` are optional.  See
            `~astroquery.mpc.MPC.get_ephemeris` for defaults.

            Epochs, including ``start`` and ``stop``, are
            `~astropy.time.Time` objects, anything that can initialize
            a ``Time`` object, or a float indicating a Julian date.
            Unless otherwise specified, the UTC scale is assumed.

            ``step`` should be an integer in units of seconds,
            minutes, hours, or days.  Anything that can initialize an
            `~astropy.units.Quantity` object is allowed.

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
        >>> epochs = {'start': '2019-01-01', 'step': '1d', 'number': 365}
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

        if isinstance(epochs, dict):
            start = epochs['start']  # required
            step = epochs.get('step')
            stop = epochs.get('stop')
            number = epochs.get('number')

            if isinstance(start, (float, int)):
                start = Time(start, format='jd', scale='utc')

            if isinstance(stop, (float, int)):
                stop = Time(stop, format='jd', scale='utc')

            if step is not None:
                step = u.Quantity(step)
                if step.unit not in (u.d, u.h, u.min, u.s):
                    raise ValueError(
                        'step must have units of days, hours, minutes,'
                        ' or seconds')

            if stop is not None:
                if step is None:
                    raise ValueError(
                        'step is required when start and stop are provided')

                # start and stop both defined, estimate number of steps
                dt = (Time(stop).jd - Time(start).jd) * u.d
                number = int((dt / step).decompose()) + 1
        else:
            start = None

            if epochs is None:
                epochs = Time.now()

            if not iterable(epochs) or isinstance(epochs, str):
                epochs = [epochs]

            # check for Julian dates
            for i in range(len(epochs)):
                if isinstance(epochs[i], (float, int)):
                    epochs[i] = Time(epochs[i], format='jd', scale='utc')

        # append ephemerides table for each targetid
        all_eph = None
        for targetid in targetids:

            # get ephemeris
            if start is None:
                eph = []
                for i in range(len(epochs)):
                    e = MPC.get_ephemeris(targetid, location=location,
                                          start=Time(epochs[i], scale='utc'),
                                          number=1, **kwargs)
                    e['Date'] = e['Date'].iso  # for vstack to work
                    eph.append(e)
                eph = vstack(eph)
                eph['Date'] = Time(eph['Date'], scale='utc')
            else:
                eph = MPC.get_ephemeris(targetid, location=location,
                                        start=start, step=step,
                                        number=number, **kwargs)

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
            only support single epochs;
            a dictionary including keywords ``start`` and ``stop``, as well
            as either ``step`` or ``number``, can be used to generate a range
            of epochs. ``start`` and ``stop`` have to be
            `~astropy.time.Time` objects or strings in ISO format
            (``'YYYY-MM-DD HH:MM:SS'``). If ``step`` is provided as a
            string, the interval
            defined by ``start``
            and ``stop`` is split into time steps defined by ``step``. If
            ``number`` is provided as an integer, the
            interval defined by
            ``start`` and ``stop`` is split into ``number`` equidistant
            intervals. If ``None`` is
            provided, current date and time are
            used. Default: ``None``
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

        Examples
        --------
        >>> from sbpy.data import Ephem
        >>> from astropy.time import Time
        >>> epoch = Time('2018-05-14', scale='utc')
        >>> eph = Ephem.from_horizons('ceres', epochs=epoch) # doctest: +SKIP
        """

        # modify epoch input to make it work with astroquery.imcce.Miriade
        if epochs is None:
            epochs = {'start': Time.now().jd}
        elif isinstance(epochs, Time):
            epochs = {'start': epochs}
        elif isinstance(epochs, dict):
            if 'number' in epochs:
                # turn interval/number into step size based on full minutes
                epochs['step'] = '{:d}m'.format(
                    int((Time(epochs['stop'])-Time(epochs['start'])).jd *
                        1440 / (epochs['number']-1)))
            elif 'step' in epochs:
                # parse 'step'
                stepsize_s = (float(epochs['step'][:-1]) *
                              {'s': 1, 'm': 60, 'h': 3600, 'd': 86400}[
                                  epochs['step'][-1]])
                epochs['number'] = ((Time(epochs['stop']) -
                                     Time(epochs['start'])).jd *
                                    86400/stepsize_s)+1

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
            if 'step' not in epochs and 'number' not in epochs:
                # single epoch
                try:
                    # load ephemerides using astroquery.imcce
                    eph = query.get_ephemerides(targetname=targetid,
                                                objtype=objtype,
                                                location=location,
                                                epoch=epochs['start'],
                                                **kwargs)
                except ValueError as e:
                    raise RuntimeError(
                        ('Error raised by astroquery.imcce: {:s}\n'
                         'The following query was attempted: {:s}').format(
                             str(e), query.uri))
            else:
                # range of epochs
                try:
                    # load ephemerides using astroquery.imcce
                    eph = query.get_ephemerides(
                        targetname=targetid, objtype=objtype,
                        location=location, epoch=epochs['start'],
                        epoch_step=epochs['step'],
                        epoch_nsteps=epochs['number'],
                        **kwargs)
                except ValueError as e:
                    raise RuntimeError(
                        ('Error raised by astroquery.imcce: {:s}\n'
                         'The following query was attempted: {:s}').format(
                             str(e), query.uri))

            if all_eph is None:
                all_eph = eph
            else:
                all_eph = vstack([all_eph, eph])

        self = cls.from_table(all_eph)

        # turn epochs into astropy.time.Time and apply timescale
        timescale = 'UTC'
        if 'timescale' in kwargs:
            timescale = kwargs['timescale']
        self.table['epoch'] = Time(self.table['epoch'],
                                   format='jd', scale=timescale.lower())

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem.from_miriade',
                         {'data service':
                          'http://vo.imcce.fr/webservices/miriade/'})

        return self

    @classmethod
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
            * epoch (``'epoch'``) in JD
            * absolute magnitude (``'H'``) in mag
            * photometric phase slope (``'G'``)

        epochs : `~astropy.time.Time` object or iterable thereof, optional
            Epochs of elements to be queried; a list, tuple or
            `~numpy.ndarray` of `~astropy.time.Time` objects or Julian
            Dates as floats should be used for a number of discrete
            epochs; if ``None`` is provided, current date and
            time are used. Default: ``None``
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

        # extract time scale
        timescale = orb.table['epoch'][0].scale.upper()

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
            raise ValueError(
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

        # modify epochs input to make it work with pyoorb
        if epochs is None:
            epochs = Time.now()
        # elif isinstance(epochs, Time):
        #     epochs = [Time(epochs)]
        # elif isinstance(epochs, (float, int)):
        #     epochs = [Time(epochs, format='jd')]
        # elif isinstance(epochs, str):
        #     epochs = [Time(epochs, format='iso')]
        # elif isinstance(epochs, (list, tuple, ndarray)):
        #     new_epochs = [None] * len(epochs)
        #     for i in range(len(epochs)):
        #         if isinstance(epochs[i], Time):
        #             new_epochs[i] = epochs[i]
        #         else:
        #             new_epochs[i] = Time(epochs[i], format='jd')
        #     epochs = new_epochs

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
        else:
            raise ValueError('only \'full\' or \'basic\' allowed for scope')

        if err != 0:
            RuntimeError('pyoorb failed with error code {:d}'.format(err))

        # reorder data on per-column basis and apply units
        oo_eph_col = hstack([oo_eph.transpose()[:, :, i]
                             for i in range(oo_eph.shape[0])]).tolist()
        oo_eph_col_u = []
        for i, col in enumerate(oo_eph_col):
            oo_eph_col_u.append(Ephem._unit_apply(col,
                                                  conf.oorb_ephem_units[i]))
        ephem = self.from_columns(oo_eph_col_u,
                                  names=conf.oorb_ephem_fields)

        # add targetname column
        ephem.table.add_column(Column(data=sum([[orb['targetname'][i]] *
                                                len(epochs) for i in
                                                range(len(orb.table))],
                                               []),
                                      name='targetname'),
                               index=0)

        # convert MJD to astropy.time.TimeJulian Date
        ephem.table['ephem'] = Time(ephem['MJD'], format='mjd',
                                    scale=timescale.lower())
        ephem.table.remove_column('MJD')

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem.from_oo',
                         {'method': '2009M&PS...44.1853G',
                          'implementation': 'https://github.com/oorb/oorb'})

        return ephem
