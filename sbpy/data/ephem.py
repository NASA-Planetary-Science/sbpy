# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Ephem Module
======================

Class for storing and querying ephemerides

created on June 04, 2017
"""
import os

from numpy import ndarray, array, hstack
from astropy.time import Time
from astropy.table import vstack, Column
from astroquery.jplhorizons import Horizons

from .. import bib
from .core import DataClass, conf

__all__ = ['Ephem']


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
        epochs : `~astropy.time.Time` object or iterable thereof, or dictionary, optional
            Epochs of elements to be queried; a list, tuple or
            `~numpy.ndarray` of `~astropy.time.Time` objects or Julian
            Dates as floats should be used for a number of discrete
            epochs; a dictionary including keywords ``start``,
            ``step``, and ``stop`` can be used to generate a range of
            epochs (see
            `~astroquery.jplhorizons.HorizonsClass.Horizons.ephemerides`
            for details); if ``None`` is provided, current date and
            time are used. Default: ``None``
        location : str, optional, default ``'500'`` (geocentric)
            Location of the observer.
        **kwargs : optional
            Arguments that will be provided to
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
            epochs = [Time(epochs).jd]
        elif isinstance(epochs, dict):
            for key, val in epochs.items():
                if isinstance(val, Time):
                    val.format = 'iso'
                    val.out_subfmt = 'date_hm'
                    epochs[key] = val.value
                else:
                    epochs[key] = epochs[key]
        elif isinstance(epochs, (list, tuple, ndarray)):
            new_epochs = [None] * len(epochs)
            for i in range(len(epochs)):
                if isinstance(epochs[i], Time):
                    new_epochs[i] = epochs[i].jd
                else:
                    new_epochs[i] = epochs[i]
            epochs = new_epochs

        # if targetids is a list, run separate Horizons queries and append
        if not isinstance(targetids, (list, ndarray, tuple)):
            targetids = [targetids]

        # append ephemerides table for each targetid
        all_eph = None
        for targetid in targetids:

            # load ephemerides using astroquery.jplhorizons
            obj = Horizons(id=targetid, id_type=id_type,
                           location=location, epochs=epochs)
            eph = obj.ephemerides(**kwargs)

            # workaround for current version of astroquery to make
            # column units compatible with astropy.table.QTable
            # should really change '---' units to None in
            # astroquery.jplhorizons.__init__.py
            for column_name in eph.columns:
                if eph[column_name].unit == '---':
                    eph[column_name].unit = None

            if all_eph is None:
                all_eph = eph
            else:
                all_eph = vstack([all_eph, eph])

        # identify time scales returned by Horizons query
        timescales = array(['UTC'] * len(all_eph))
        timescales[all_eph['datetime_jd'] < 2437665.5] = 'UT1'
        # according to Horizons documentation
        all_eph.add_column(Column(timescales, name='timescale'))

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem.from_horizons',
                         {'data service': '1996DPS....28.2504G'})

        return cls.from_table(all_eph)

    @classmethod
    def from_mpc(cls, targetid, epoch, observatory='500'):
        """
        Load ephemerides from the
        `Minor Planet Center <http://minorplanetcenter.net>`_.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default ``None``
            epoch of elements; if ``None`` is provided, current date is used
        observatory : str, optional, default ``'500'`` (geocentric)
            location of observer

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> from astropy.time import Time  # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc')  # doctest: +SKIP
        >>> eph = Ephem.from_mpc('ceres', '568', epoch)  # doctest: +SKIP

        not yet implemented

        """

    def report_to_mpc():
        """
        Format ephemerides as a report to the
        `Minor Planet Center <http://minorplanetcenter.net>`_.

        Returns
        -------
        list of strings

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> eph = Ephem.from_array...  # doctest: +SKIP
        >>> report = eph.report_to_mpc()  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_imcce(cls, targetid, epoch, observatory='500'):
        """
        Load orbital elements from
        `IMCCE <http://vo.imcce.fr/webservices/miriade/>`_.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default ``None``
            epoch of elements; if ``None`` is provided, current date is used
        observatory : str, optional, default '500' (geocentric)
            location of observer

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> from astropy.time import Time  # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc')  # doctest: +SKIP
        >>> eph = Ephem.from_imcce('ceres', '568', epoch)  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_lowell(cls, targetid, epoch, observatory='500'):
        """
        Load orbital elements from
        Lowell Observatory's `astorb <http://asteroid.lowell.edu/>`_.

        Parameters
        ----------
        targetid : str, mandatory
            target identifier
        epochs : astropy Time instance or iterable, optional, default ``None``
            epoch of elements; if ``None`` is provided, current date is used
        observatory : str, optional, default '500' (geocentric)
            location of observer

        Returns
        -------
        `~Ephem` object

        Examples
        --------
        >>> from sbpy.data import Ephem  # doctest: +SKIP
        >>> from astropy.time import Time  # doctest: +SKIP
        >>> epoch = Time('2018-05-14', scale='utc')  # doctest: +SKIP
        >>> eph = Ephem.from_lowell('ceres', '568', epoch)  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_pyephem(cls, orb, location, epoch):
        """
        Derives ephemerides using
        `PyEphem <http://rhodesmill.org/pyephem/>`_.

        Parameters
        ----------
        orb : ``~Orbit`` object, mandatory
            complete set of orbital elements
        location : str or dictionary, mandatory
            name of location or a dictionary fully describing the location
        epoch : astropy ``Time`` object

        Examples
        --------
        >>> from sbpy.data import Ephem, Orbit  # doctest: +SKIP
        >>> orb = Orbit.from_...  # doctest: +SKIP
        >>> eph = Ephem.from_pyephem(orb,  # doctest: +SKIP
        >>>                          location={'name': 'Flagstaff',  # doctest: +SKIP
        >>>                                    'geolon': 35.199167,  # doctest: +SKIP
        >>>                                    'geolat': -111.631111,  # doctest: +SKIP
        >>>                                    'altitude': '2106'},  # doctest: +SKIP
        >>>                          epoch=epoch)  # doctest: +SKIP

        not yet implemented

        """

    @classmethod
    def from_oo(self, orbit, epochs=None, location='500', scope='full',
                timescale='UTC', dynmodel='N', ephfile='de430'):
        """Uses pyoorb to derive ephemerides from an `~Orbit` object. For a
        list of output parameters, please read the `pyoorb documentation
        <https://github.com/oorb/oorb/tree/master/python>`_.

        Parameters
        ----------
        orbit : `~Orbit` object
            orbit can contain any number of orbits. Required fields are XXX.
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
        timescale : str, optional
            Timescale to be used in the computation; the following
            values are allowed: ``'UTC'``, ``'UT1'``, ``'TT'``,
            ``'TAI'``. Default: ``'UTC'``
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
        >>> epochs = Time.now().jd + np.arange(0, 10, 1/24)
        >>> ceres = Orbit.from_horizons('1')
        >>> eph = Ephem.from_oo(ceres, epochs=epochs, location='G37')  # doctest: +SKIP
        >>> print(eph.table)  # doctest: +SKIP
        targetname      MJD [1]       ...        obsy [1]              obsz [1]
                           d          ...           AU                    AU
        ---------- ------------------ ... --------------------- ----------------------
           1 Ceres 58374.720415079966 ...   -0.1640418731222332 1.3660753531152814e-05
           1 Ceres  58374.76208174648 ...  -0.16334416599555382 1.6732994041007698e-05
           1 Ceres 58374.803748413455 ...  -0.16264729902661218 2.0200328928084155e-05
           1 Ceres 58374.845415079966 ...  -0.16195072092478624 2.3823231905778508e-05
           1 Ceres  58374.88708174648 ...  -0.16125385509757997  2.735153478080482e-05
           1 Ceres 58374.928748413455 ...  -0.16055613920683476 3.0541568772989025e-05
                ...                ... ...                   ...                    ...
           1 Ceres 58384.428748413455 ... 0.0016096754330388497  9.924120661052244e-06
           1 Ceres 58384.470415079966 ... 0.0023287044344341605   7.69766111133525e-06
           1 Ceres  58384.51208174648 ... 0.0030458232636104473  6.300640241761616e-06
           1 Ceres 58384.553748413455 ...  0.003760809893911351 5.8280310798125914e-06
           1 Ceres 58384.595415079966 ...  0.004473588211662766  6.311456253324348e-06
           1 Ceres  58384.63708174648 ...  0.005184233254950517  7.717021060406424e-06
           1 Ceres 58384.678748413455 ...  0.005892966025131429  9.947635868821306e-06
        Length = 240 rows
        """

        import pyoorb

        # initialize pyoorb
        ephfile = os.path.join(os.getenv('OORB_DATA'), ephfile+'.dat')
        pyoorb.pyoorb.oorb_init(ephfile)

        # identify orbit type based on available table columns
        orbittype = None
        for testtype in ['KEP', 'COM', 'CART']:
            try:
                orbit._translate_columns(
                    conf.oorb_orbit_fields[testtype][1:6])
                orbittype = testtype
                break
            except KeyError:
                pass

        if orbittype is None:
            raise ValueError(
                'orbit type cannot be determined from elements')

        # modify epochs input to make it work with pyoorb
        if epochs is None:
            epochs = [Time.now()]
        elif isinstance(epochs, Time):
            epochs = [Time(epochs)]
        elif isinstance(epochs, (list, tuple, ndarray)):
            new_epochs = [None] * len(epochs)
            for i in range(len(epochs)):
                if isinstance(epochs[i], Time):
                    new_epochs[i] = epochs[i]
                else:
                    new_epochs[i] = Time(epochs[i], format='jd')
            epochs = new_epochs
        epochs = list(zip([epoch.jd-2400000.5 for epoch in epochs],
                          [conf.oorb_timeScales[timescale]]*len(epochs)))

        if scope == 'full':
            oo_eph, err = pyoorb.pyoorb.oorb_ephemeris_full(
                orbit._to_oo(),
                location,
                epochs,
                dynmodel)
        elif scope == 'basic':
            oo_eph, err = pyoorb.pyoorb.oorb_ephemeris_basic(
                orbit._to_oo(),
                location,
                epochs,
                dynmodel)
        else:
            raise ValueError('only \'full\' or \'basic\' allowed for scope')

        if err != 0:
            RuntimeError('pyoorb failed with error code {:d}'.format(err))

        # reorder data in Orbit object
        ephem = self.from_array(hstack([oo_eph.transpose()[:, :, i]
                                        for i in range(oo_eph.shape[0])]),
                                names=conf.oorb_ephem_fields)

        # apply units
        for i, col in enumerate(ephem.column_names):
            ephem[col].unit = conf.oorb_ephem_units[i]

        # add targetname column
        ephem.table.add_column(Column(data=sum([[orbit['targetname'][i]] *
                                                len(epochs) for i in
                                                range(len(orbit.table))],
                                               []),
                                      name='targetname'),
                               index=0)

        # remove trueanom column for now as it only holds a dummy value
        ephem.table.remove_column('trueanom')

        if bib.status() is None or bib.status():
            bib.register('sbpy.data.Ephem.from_oo',
                         {'method': '2009M&PS...44.1853G',
                          'implementation': 'https://github.com/oorb/oorb'})

        return ephem
