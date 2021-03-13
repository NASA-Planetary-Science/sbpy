===========
Using Ephem
===========

As shown above (:ref:`How to use Data Containers`),
`~sbpy.data.Ephem` objects can be created on the fly. However,
`~sbpy.data.Ephem` can also be used to access ephemerides information
from remote services with a largely uniform API.

For instance, the following few lines will query
ephemerides for asteroid Ceres on a given date and for the position of
Mauna Kea Observatory (IAU observatory code ``568``) from the `JPL Horizons service <https://ssd.jpl.nasa.gov/horizons.cgi>`_:

    >>> from sbpy.data import Ephem
    >>> from astropy.time import Time
    >>> epoch = Time('2018-08-03 14:20', scale='utc') # time in UT
    >>> eph = Ephem.from_horizons('Ceres',
    ...                           location='568',
    ...                           epochs=epoch)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA +SKIP
    <QTable masked=True length=1>
    targetname    H       G    solar_presence ...  PABLon  PABLat       epoch
                 mag                          ...   deg     deg
       str7    float64 float64      str1      ... float64 float64       object
    ---------- ------- ------- -------------- ... ------- ------- -----------------
       1 Ceres    3.34    0.12                ... 171.275  9.3473 2458334.097222222
    >>> eph.field_names  # doctest: +REMOTE_DATA
    ['targetname', 'H', 'G', 'solar_presence', 'flags', 'RA', 'DEC', 'RA_app', 'DEC_app', 'RA*cos(Dec)_rate', 'DEC_rate', 'AZ', 'EL', 'AZ_rate', 'EL_rate', 'sat_X', 'sat_Y', 'sat_PANG', 'siderealtime', 'airmass', 'magextinct', 'V', 'surfbright', 'illumination', 'illum_defect', 'sat_sep', 'sat_vis', 'ang_width', 'PDObsLon', 'PDObsLat', 'PDSunLon', 'PDSunLat', 'SubSol_ang', 'SubSol_dist', 'NPole_ang', 'NPole_dist', 'EclLon', 'EclLat', 'r', 'r_rate', 'delta', 'delta_rate', 'lighttime', 'vel_sun', 'vel_obs', 'elong', 'elongFlag', 'alpha', 'lunar_elong', 'lunar_illum', 'sat_alpha', 'sunTargetPA', 'velocityPA', 'OrbPlaneAng', 'constellation', 'TDB-UT', 'ObsEclLon', 'ObsEclLat', 'NPole_RA', 'NPole_DEC', 'GlxLon', 'GlxLat', 'solartime', 'earth_lighttime', 'RA_3sigma', 'DEC_3sigma', 'SMAA_3sigma', 'SMIA_3sigma', 'Theta_3sigma', 'Area_3sigma', 'RSS_3sigma', 'r_3sigma', 'r_rate_3sigma', 'SBand_3sigma', 'XBand_3sigma', 'DoppDelay_3sigma', 'true_anom', 'hour_angle', 'alpha_true', 'PABLon', 'PABLat', 'epoch']

`~sbpy.data.Ephem.from_horizons` uses one or more target names, an
observer location in the form of an IAU observatory code, and a list
of discrete epochs or a range of epochs defined in a dictionary (see
`~sbpy.data.Ephem.from_horizons`) to query the JPL Horizons
service. Epochs have to be provided in the form of `~astropy.time.Time`
objects. The column names in the data table can be inquired using
`~sbpy.data.DataClass.field_names`.

`~sbpy.data.Ephem.from_horizons` is actually a wrapper around
`~astroquery.jplhorizons.HorizonsClass.ephemerides`. This function
conveniently combines the creation of a
`~astroquery.jplhorizons.HorizonsClass` query and the actual
ephemerides information retrieval into a single function. Additional
optional parameters provided to `~sbpy.data.Ephem.from_horizons` are
directly passed on to
`~astroquery.jplhorizons.HorizonsClass.ephemerides`, maintaining the
full flexibility of the latter function:

    >>> import astropy.units as u
    >>> epoch1 = Time('2018-08-03 14:20', scale='utc')
    >>> epoch2 = Time('2018-08-04 07:30', scale='utc')
    >>> eph = Ephem.from_horizons('Ceres',
    ...                           location='568',
    ...                           epochs={'start': epoch1,
    ...                                   'stop': epoch2,
    ...                                   'step': 10*u.minute},
    ...                           skip_daylight=True)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA +SKIP
    <QTable masked=True length=26>
    targetname    H       G    ...  PABLon   PABLat       epoch
                 mag           ...   deg      deg
       str7    float64 float64 ... float64  float64       object
    ---------- ------- ------- ... -------- ------- -----------------
       1 Ceres    3.34    0.12 ...  171.275  9.3473 2458334.097222222
       1 Ceres    3.34    0.12 ... 171.2774  9.3472 2458334.104166667
       1 Ceres    3.34    0.12 ... 171.2798  9.3471 2458334.111111111
       1 Ceres    3.34    0.12 ... 171.2822   9.347 2458334.118055556
       1 Ceres    3.34    0.12 ... 171.2846  9.3469       2458334.125
       1 Ceres    3.34    0.12 ... 171.2869  9.3468 2458334.131944444
           ...     ...     ... ...      ...     ...               ...
       1 Ceres    3.34    0.12 ... 171.5076  9.3369 2458334.777777778
       1 Ceres    3.34    0.12 ... 171.5099  9.3368 2458334.784722222
       1 Ceres    3.34    0.12 ... 171.5123  9.3367 2458334.791666667
       1 Ceres    3.34    0.12 ... 171.5147  9.3366 2458334.798611111
       1 Ceres    3.34    0.12 ... 171.5171  9.3365 2458334.805555556
       1 Ceres    3.34    0.12 ... 171.5195  9.3364      2458334.8125

Note that ``skip_daylight`` is an optional parameter of
`~astroquery.jplhorizons.HorizonsClass.ephemerides` and it can be used
here as well. An additional feature of
`~sbpy.data.Ephem.from_horizons` is that you can automatically
concatenate queries for a number of objects:

    >>> eph = Ephem.from_horizons(['Ceres', 'Pallas', 12893, '1983 SA'],
    ...                           location='568',
    ...                           epochs=epoch1)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA +SKIP
    <QTable masked=True length=4>
            targetname            H       G    ...  PABLat        epoch
                                 mag           ...   deg
              str26            float64 float64 ... float64        object
    -------------------------- ------- ------- ... -------- -----------------
                       1 Ceres    3.34    0.12 ...   9.3473 2458334.097222222
                      2 Pallas    4.13    0.11 ... -20.1396 2458334.097222222
     12893 Mommert (1998 QS55)    13.9    0.15 ...  -2.0567 2458334.097222222
    3552 Don Quixote (1983 SA)    12.9    0.15 ...  13.3365 2458334.097222222

    
Please be aware that these queries are not simultaneous. The more
targets you query, the longer the query will take. Furthermore, keep
in mind that asteroids and comets have slightly different table
layouts (e.g., different magnitude systems: ``T-mag`` and ``N-mag``
instead of ``V-mag``), which will complicate the interpretation of the
data. It might be safest to query asteroids and comets separately.

Also note that the two examples shown above use different ways to
define epochs. The first example uses a dictionary that defines a
``start`` and ``stop`` epoch, as well as a ``step`` size (see
`~sbpy.data.Ephem.from_horizons` for
details). Instead of a ``step`` size, the user can also provide a
``number`` of steps as an integer; ``step`` is then automatically
derived from the interval and ``number`` in units of full minutes. The
second example uses a specific epoch as input. The
`~astropy.time.Time` object provided to ``epochs`` can also be
initialized with a list of epochs, querying multiple epochs at the
same time. Note that the total number of epochs queried using this
option should be less than a few hundred to prevent corruption of the
query (see `~astroquery.jplhorizons.HorizonsClass.ephemerides` for
details).

Observer locations can be defined as strings using official `IAU
observatory codes
<https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html>`__ or
using `~astropy.coordinates.EarthLocation` as shown in the following
example:

    >>> from astropy.coordinates import EarthLocation
    >>> lowell = EarthLocation.of_site('Lowell Observatory')  # doctest: +SKIP
    >>> eph = Ephem.from_horizons(1, epochs=Time('2018-01-01', format='iso'),
    ... 			  location=lowell) # doctest: +SKIP
    >>> eph # doctest: +REMOTE_DATA +SKIP
    <QTable masked=True length=1>
    targetname    H       G    solar_presence ...  PABLon   PABLat   epoch  
                 mag                          ...   deg      deg            
       str7    float64 float64      str1      ... float64  float64   object 
    ---------- ------- ------- -------------- ... -------- ------- ---------
       1 Ceres    3.34    0.12              * ... 130.4303  9.2004 2458119.5

Offering almost identical functionality, the
`~sbpy.data.Ephem.from_mpc` method will retrieve ephemerides from the
`Minor Planet Center <https://minorplanetcenter.net/>`_:

    >>> eph = Ephem.from_mpc('2P', location='568',
    ...                      epochs={'start': Time('2018-10-22'),
    ...                              'stop': Time('2018-10-26'),
    ...                              'step': 1*u.day})  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA +SKIP
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

Finally, `~sbpy.data.Ephem.from_miriade` will retrieve ephemerides
from the `Miriade ephemeris generator
<http://vo.imcce.fr/webservices/miriade/>`_ at `IMCCE
<https://www.imcce.fr/>`_:

    >>> eph = Ephem.from_miriade('2P', objtype='comet', location='568',
    ...                          epochs={'start': Time('2018-10-22'),
    ...                                  'stop': Time('2018-10-26'),
    ...                                  'step': 1*u.day})  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA +SKIP
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
    
Ephemerides can also be derived from `~sbpy.data.Orbit` objects using
`sbpy`'s interface to `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_ with the function
`~sbpy.data.Ephem.from_oo`. The following example computes
ephemerides for the next ten days in steps of 1 hr for Ceres as seen
from the Discovery Channel Telescope:

    >>> import numpy as np
    >>> from sbpy.data import Orbit, Ephem
    >>> from astropy.time import Time
    >>> epochs = Time(Time.now().jd + np.arange(0, 10, 1/24), format='jd')
    >>> ceres = Orbit.from_horizons('1')  # doctest: +REMOTE_DATA
    >>> eph = Ephem.from_oo(ceres, epochs, 'G37') # doctest: +SKIP 
    >>> eph # doctest: +SKIP 
    <QTable length=240>
    targetname         RA         ...      trueanom            epoch       
                      deg         ...        deg                           
       str7         float64       ...      float64             object      
    ---------- ------------------ ... ------------------ ------------------
       1 Ceres 238.56187075007446 ...  105.8270438687299 2458694.6423231447
       1 Ceres   238.564318627966 ... 105.83566067245822  2458694.683989811
       1 Ceres 238.56680284927273 ...  105.8442772820886  2458694.725656478
       1 Ceres 238.56933812666867 ...  105.8528936974433 2458694.7673231447
       1 Ceres 238.57193638137088 ...  105.8615099186335  2458694.808989811
       1 Ceres 238.57460592776462 ... 105.87012594577034  2458694.850656478
           ...                ... ...                ...                ...
       1 Ceres  239.4677754274348 ... 107.83811369526742 2458704.3923231447
       1 Ceres  239.4726928414698 ...   107.846685468736  2458704.433989811
       1 Ceres 239.47756694312102 ... 107.85525705166283  2458704.475656478
       1 Ceres 239.48240809475683 ...  107.8638284438719 2458704.5173231447
       1 Ceres 239.48722955376766 ... 107.87239964547449  2458704.558989811
       1 Ceres 239.49204656314026 ... 107.88097065658197  2458704.600656478

     
The properties computed by pyoorb and listed in the resulting table
are defined in the `pyoorb documentation
<https://github.com/oorb/oorb/tree/master/python>`_. Note that this function requires pyoorb to be installed, which is not a requirement for `sbpy`.


