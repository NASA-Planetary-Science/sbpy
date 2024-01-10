=====================================
Ephemeris Objects (`sbpy.data.Ephem`)
=====================================

As described in :ref:`How to use Data Containers`, `~sbpy.data.Ephem` objects
can be created on the fly. However, `~sbpy.data.Ephem` can also be used to
access ephemeris information from other sources with a largely uniform API.


JPL Horizons (`~sbpy.data.Ephem.from_horizons`)
-----------------------------------------------

`~sbpy.data.Ephem.from_horizons` uses one or more target names, an observer
location, and the ephemeris epoch.  For instance, the following few lines will
query for ephemerides of asteroid Ceres on a given date and for the position of
Mauna Kea Observatory (IAU observatory code 568) from the `JPL Horizons service
<https://ssd.jpl.nasa.gov/horizons/>`_:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> from sbpy.data import Ephem
    >>> from astropy.time import Time
    >>> epoch = Time('2018-08-03 14:20', scale='utc') # time in UT
    >>> eph = Ephem.from_horizons('Ceres',
    ...                           location='568',
    ...                           epochs=epoch)
    >>> eph['epoch', 'ra', 'dec', 'rh', 'delta', 'phase']  # doctest: +SKIP
    <QTable length=1>
          epoch           RA      DEC          r             delta        alpha 
                         deg      deg          AU              AU          deg  
           Time        float64  float64     float64         float64      float64
    ----------------- --------- -------- -------------- ---------------- -------
    2458334.097222222 169.23471 13.40412 2.573882458571 3.34213404869399  12.978

The full column name list in the data table can be retrieved with the
`~sbpy.data.DataClass.field_names` property:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> eph.field_names
    ['targetname', 'H', 'G', 'solar_presence', 'flags', 'RA', 'DEC', 'RA_app', 'DEC_app', 'RA*cos(Dec)_rate', 'DEC_rate', 'AZ', 'EL', 'AZ_rate', 'EL_rate', 'sat_X', 'sat_Y', 'sat_PANG', 'siderealtime', 'airmass', 'magextinct', 'V', 'surfbright', 'illumination', 'illum_defect', 'sat_sep', 'sat_vis', 'ang_width', 'PDObsLon', 'PDObsLat', 'PDSunLon', 'PDSunLat', 'SubSol_ang', 'SubSol_dist', 'NPole_ang', 'NPole_dist', 'EclLon', 'EclLat', 'r', 'r_rate', 'delta', 'delta_rate', 'lighttime', 'vel_sun', 'vel_obs', 'elong', 'elongFlag', 'alpha', 'lunar_elong', 'lunar_illum', 'sat_alpha', 'sunTargetPA', 'velocityPA', 'OrbPlaneAng', 'constellation', 'TDB-UT', 'ObsEclLon', 'ObsEclLat', 'NPole_RA', 'NPole_DEC', 'GlxLon', 'GlxLat', 'solartime', 'earth_lighttime', 'RA_3sigma', 'DEC_3sigma', 'SMAA_3sigma', 'SMIA_3sigma', 'Theta_3sigma', 'Area_3sigma', 'RSS_3sigma', 'r_3sigma', 'r_rate_3sigma', 'SBand_3sigma', 'XBand_3sigma', 'DoppDelay_3sigma', 'true_anom', 'hour_angle', 'alpha_true', 'PABLon', 'PABLat', 'epoch']


Time options
^^^^^^^^^^^^

In the above example a specific epoch was specified, but multiple epochs may be
requested, or even a range of epochs.  All time references are
`~astropy.time.Time` objects.

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> epochs = Time(['2022-06-04', '2023-06-04'])
    >>> eph = Ephem.from_horizons('Ceres', location='568', epochs=epochs)
    >>> eph['epoch', 'rh', 'delta']  # doctest: +SKIP
    <QTable length=2>
      epoch         r             delta      
                    AU              AU       
       Time      float64         float64     
    --------- -------------- ----------------
    2459734.5  2.60719514844 3.49018226961962
    2460099.5 2.603599947555 2.18175266069806

Note that the total number of epochs queried using this option should be less
than a few hundred to prevent corruption of the query (see
`~astroquery.jplhorizons.HorizonsClass.ephemerides` for details).

To specify a range of epochs:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> import astropy.units as u
    >>>
    >>> epochs = {'start': Time('2022-01-01'),
    ...           'stop': Time('2023-01-01'),
    ...           'step': 30 * u.day}
    >>> eph = Ephem.from_horizons('Ceres', location='568', epochs=epochs)
    >>> eph['epoch', 'rh', 'delta']  # doctest: +SKIP
    <QTable length=13>
      epoch         r             delta      
                    AU              AU       
       Time      float64         float64     
    --------- -------------- ----------------
    2459580.5 2.718302197006  1.9062817959135
    2459610.5 2.694384101844 2.22554330199232
    ...
    2459940.5 2.549833352462 2.30003114258051

As an alternative to ``'step'`` one could specify the number of epochs with
``''number'``.


Mulitple targets
^^^^^^^^^^^^^^^^

An additional feature of `~sbpy.data.Ephem.from_horizons` is that you can
automatically concatenate queries for a number of objects:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> epoch1 = Time('2018-08-03 14:20')
    >>> eph = Ephem.from_horizons(['Ceres', 'Pallas', 12893, '1983 SA'],
    ...                           location='568',
    ...                           epochs=epoch1)
    >>> eph  # doctest: +SKIP
    <QTable masked=True length=4>
            targetname            H       G    ...  PABLat        epoch
                                 mag           ...   deg
              str26            float64 float64 ... float64        object
    -------------------------- ------- ------- ... -------- -----------------
                       1 Ceres    3.34    0.12 ...   9.3473 2458334.097222222
                      2 Pallas    4.13    0.11 ... -20.1396 2458334.097222222
     12893 Mommert (1998 QS55)    13.9    0.15 ...  -2.0567 2458334.097222222
    3552 Don Quixote (1983 SA)    12.9    0.15 ...  13.3365 2458334.097222222

    
Please be aware that these queries are not simultaneous. The more targets you
query, the longer the query will take. Furthermore, keep in mind that asteroids
and comets have slightly different table layouts (e.g., different magnitude
systems: ``T-mag`` and ``N-mag`` instead of ``V-mag``), which will complicate
the interpretation of the data. It might be safest to query asteroids and comets
separately.


Observer locations
^^^^^^^^^^^^^^^^^^

Observer locations can be defined as strings using official `IAU observatory
codes <https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html>`__ as above,
or by using `~astropy.coordinates.EarthLocation` as shown in the following
example:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> from astropy.coordinates import EarthLocation
    >>> lowell = EarthLocation.of_site('Lowell Observatory')
    >>> eph = Ephem.from_horizons(1, epochs=Time('2018-01-01'),
    ...                           location=lowell)
    >>> eph # doctest: +SKIP
    <QTable masked=True length=1>
    targetname    H       G    solar_presence ...  PABLon   PABLat   epoch  
                 mag                          ...   deg      deg            
       str7    float64 float64      str1      ... float64  float64   object 
    ---------- ------- ------- -------------- ... -------- ------- ---------
       1 Ceres    3.34    0.12              * ... 130.4303  9.2004 2458119.5


Optional parameters
^^^^^^^^^^^^^^^^^^^

`~sbpy.data.Ephem.from_horizons` is actually a wrapper around
`astroquery.jplhorizons.HorizonsClass`.  Additional optional parameters provided
to `~sbpy.data.Ephem.from_horizons` are directly passed on to
`astroquery.jplhorizons.HorizonsClass.ephemerides`, maintaining the full
flexibility of the latter function.  For example one may use the
``skip_daylight`` keyword argument:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> epoch1 = Time('2018-08-03 14:20', scale='utc')
    >>> epoch2 = Time('2018-08-04 07:30', scale='utc')
    >>> eph = Ephem.from_horizons('Ceres',
    ...                           location='568',
    ...                           epochs={'start': epoch1,
    ...                                   'stop': epoch2,
    ...                                   'step': 10 * u.minute},
    ...                           skip_daylight=True)

Or, a common option for periodic cometary targets is to limit orbit look-ups to
the apparition closest to the epochs being queried (requires
``id_type='designation'``):

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> eph = Ephem.from_horizons('2P')   # doctest: +SKIP
    Traceback (most recent call last):
    ...
    ValueError: Ambiguous target name; provide unique id:
    Record #  Epoch-yr  >MATCH DESIG<  Primary Desig  Name  
    --------  --------  -------------  -------------  -------------------------
    90000034    1786    2P             2P              Encke
    90000035    1796    2P             2P              Encke
    90000036    1805    2P             2P              Encke
    ...
    >>> eph = Ephem.from_horizons('2P', id_type='designation', closest_apparition=True)
    >>> print(eph['targetname'])                                                                    
    targetname
    ----------
      2P/Encke


Minor Planet Center's Ephemeris Service (`~sbpy.data.Ephem.from_mpc`)
---------------------------------------------------------------------

Offering similar functionality, the `~sbpy.data.Ephem.from_mpc` method will
retrieve ephemerides from the `Minor Planet Center's Ephemeris Service
<https://minorplanetcenter.net/iau/MPEph/MPEph.html>`_:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> eph = Ephem.from_mpc('2P', location='568',
    ...                      epochs={'start': Time('2018-10-22'),
    ...                              'stop': Time('2018-10-26'),
    ...                              'step': 1*u.day})
    >>> eph  # doctest: +SKIP
    <QTable length=5>
    Targetname           Date          ... Moon distance Moon altitude
                                       ...      deg           deg
       str2             object         ...    float64       float64
    ---------- ----------------------- ... ------------- -------------
            2P 2018-10-22 00:00:00.000 ...          28.0         -33.0
            2P 2018-10-24 00:00:00.000 ...          54.0         -48.0
            2P 2018-10-25 00:00:00.000 ...          67.0         -53.0
            2P 2018-10-26 00:00:00.000 ...          81.0         -56.0
            2P 2018-10-23 00:00:00.000 ...          41.0         -41.0


IMCCE's Miriade (`~sbpy.data.Ephem.from_miriade`)
-------------------------------------------------

Finally, `~sbpy.data.Ephem.from_miriade` will retrieve ephemerides from the
`Miriade ephemeris generator <http://vo.imcce.fr/webservices/miriade/>`_ at
`Institut de Mécanique Céleste et de Calcul des Éphémérides
<https://www.imcce.fr/>`_:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> eph = Ephem.from_miriade('2P', objtype='comet', location='568',
    ...                          epochs={'start': Time('2018-10-22'),
    ...                                  'stop': Time('2018-10-26'),
    ...                                  'step': 1*u.day})
    >>> eph  # doctest: +SKIP
    <QTable masked=True length=5>
     target   epoch           RA         ...   DEC_rate    delta_rate 
                             deg         ... arcsec / min    km / s   
    bytes20   object       float64       ...   float64      float64   
    ------- --------- ------------------ ... ------------ ------------
         2P 2458413.5 329.99213124999994 ...    -0.063365   24.7933113
         2P 2458414.5 329.91132124999996 ...    -0.059361   25.0280603
         2P 2458415.5 329.83517041666664 ...    -0.055369    25.253586
         2P 2458416.5 329.76366666666667 ...    -0.051392   25.4700287
         2P 2458417.5  329.6967958333333 ...     -0.04743    25.677518


Using an orbit and OpenOrb (`~sbpy.data.Ephem.from_oo`)
-------------------------------------------------------

Ephemerides can also be derived from `~sbpy.data.Orbit` objects using `sbpy`'s
interface to `OpenOrb <https://github.com/oorb/oorb>`_ with the function
`~sbpy.data.Ephem.from_oo`. The following example computes ephemerides for the
next ten days in steps of 1 hr for Ceres as seen from the Discovery Channel
Telescope:

.. doctest-requires:: oorb

    >>> import numpy as np
    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from sbpy.data import Orbit, Ephem
    >>>
    >>> ceres = Orbit.from_dict({'targetname': 'Ceres',
    ...                          'orbtype': 'KEP',
    ...                          'a': 2.77 * u.au,
    ...                          'e': 0.0786,
    ...                          'i': 10.6 * u.deg,
    ...                          'w': 73.6 * u.deg,
    ...                          'Omega': 80.3 * u.deg,
    ...                          'M': 320.3 * u.deg,
    ...                          'epoch': Time(2459735.0, format='jd'),
    ...                          'H': 3.3 * u.mag,
    ...                          'G': 0.15})
    >>> epochs = Time('2022-06-01') + np.arange(31) * u.day
    >>> eph = Ephem.from_oo(ceres, epochs, 'G37')
    >>> print(eph)
    <QTable length=31>
    targetname         RA                DEC         ...      trueanom            epoch      
                      deg                deg         ...        deg                          
       str5         float64            float64       ...      float64              Time      
    ---------- ------------------ ------------------ ... ------------------ -----------------
         Ceres  97.53508969190534 26.840028524616123 ... 313.25073865769383 2459731.500800741
         Ceres  98.00890363845723 26.840660100535846 ...  313.4904467480157 2459732.500800741
         Ceres  98.48352873396654  26.83985349294406 ...  313.7302633914505 2459733.500800741
         Ceres    98.958932708787  26.83760580100297 ... 313.97018822933694 2459734.500800741
         Ceres  99.43508364177217 26.833914483074352 ...  314.2102209001088 2459735.500800741
         Ceres  99.91195006004624 26.828777353292317 ...  314.4503610392964 2459736.500800741
         Ceres 100.38950104500788 26.822192575950346 ...  314.6906082795292 2459737.500800741
         Ceres 100.86770636534342 26.814158657288353 ... 314.93096225053836 2459738.500800741
           ...                ...                ... ...                ...               ...
         Ceres 108.10282016659154 26.519545221745105 ... 318.54881889008544 2459753.500800741
         Ceres 108.58838133652219 26.488326782806954 ... 318.79082784334094 2459754.500800741
         Ceres 109.07422128038874 26.455671250331037 ... 319.03293679237225 2459755.500800741
         Ceres 109.56031412985443  26.42158188064386 ...  319.2751453153177 2459756.500800741
         Ceres 110.04663355942071 26.386062295091307 ...  319.5174529874914 2459757.500800741
         Ceres 110.53315292361725  26.34911647975649 ...  319.7598593813884 2459758.500800741
         Ceres 111.01984538555946 26.310748782938763 ...  320.0023640666907 2459759.500800741
         Ceres 111.50668404363827 26.270963909983408 ... 320.24496661027257 2459760.500800741
         Ceres 111.99364205884406 26.229766915255638 ... 320.48766657620695 2459761.500800741

The properties computed by pyoorb and listed in the resulting table are defined
in the `pyoorb documentation
<https://github.com/oorb/oorb/tree/master/python>`_. Note that this function
requires pyoorb to be installed, which is not a requirement for `sbpy`.