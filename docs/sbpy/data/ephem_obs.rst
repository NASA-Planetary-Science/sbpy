===================
Using Ephem and Obs
===================

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
    >>> eph  # doctest: +REMOTE_DATA
    <QTable masked=True length=1>
    targetname       datetime_str          datetime_jd    ...  PABLat timescale
						d         ...   deg
       str7             str24                float64      ... float64    str3
    ---------- ------------------------ ----------------- ... ------- ---------
       1 Ceres 2018-Aug-03 14:20:00.000 2458334.097222222 ...  9.3473       UTC

    >>> eph.field_names  # doctest: +REMOTE_DATA
    <TableColumns names=('targetname','datetime_str','datetime_jd','H','G','solar_presence','flags','RA','DEC','RA_app','DEC_app','RA*cos(Dec)_rate','DEC_rate','AZ','EL','AZ_rate','EL_rate','sat_X','sat_Y','sat_PANG','siderealtime','airmass','magextinct','V','surfbright','illumination','illum_defect','sat_sep','sat_vis','ang_width','PDObsLon','PDObsLat','PDSunLon','PDSunLat','SubSol_ang','SubSol_dist','NPole_ang','NPole_dist','EclLon','EclLat','r','r_rate','delta','delta_rate','lighttime','vel_sun','vel_obs','elong','elongFlag','alpha','lunar_elong','lunar_illum','sat_alpha','sunTargetPA','velocityPA','OrbPlaneAng','constellation','TDB-UT','ObsEclLon','ObsEclLat','NPole_RA','NPole_DEC','GlxLon','GlxLat','solartime','earth_lighttime','RA_3sigma','DEC_3sigma','SMAA_3sigma','SMIA_3sigma','Theta_3sigma','Area_3sigma','RSS_3sigma','r_3sigma','r_rate_3sigma','SBand_3sigma','XBand_3sigma','DoppDelay_3sigma','true_anom','hour_angle','alpha_true','PABLon','PABLat','timescale')>

`~sbpy.data.Ephem.from_horizons` uses one or more target names, an
observer location in the form of an IAU observatory code, and a list
of discrete epochs or a range of epochs defined in a dictionary (see
`~sbpy.data.Ephem.from_horizons`) to query the JPL Horizons
service. Due to different requirements of the JPL Horizons service for
the epoch format, we recommend to use `~astropy.time.Time`
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

    >>> epoch1 = Time('2018-08-03 14:20', scale='utc')
    >>> epoch2 = Time('2018-08-04 07:30', scale='utc')
    >>> eph = Ephem.from_horizons('Ceres',
    ...                           location='568',
    ...                           epochs={'start': epoch1,
    ...                                   'stop': epoch2,
    ...                                   'step': '10m'},
    ...                           skip_daylight=True)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable masked=True length=26>
    targetname    datetime_str      datetime_jd    ...  PABLon   PABLat timescale
					 d         ...   deg      deg
       str7          str17            float64      ... float64  float64    str3
    ---------- ----------------- ----------------- ... -------- ------- ---------
       1 Ceres 2018-Aug-03 14:20 2458334.097222222 ...  171.275  9.3473       UTC
       1 Ceres 2018-Aug-03 14:30 2458334.104166667 ... 171.2774  9.3472       UTC
       1 Ceres 2018-Aug-03 14:40 2458334.111111111 ... 171.2798  9.3471       UTC
       1 Ceres 2018-Aug-03 14:50 2458334.118055556 ... 171.2822   9.347       UTC
       1 Ceres 2018-Aug-03 15:00       2458334.125 ... 171.2846  9.3469       UTC
       1 Ceres 2018-Aug-03 15:10 2458334.131944444 ... 171.2869  9.3468       UTC
	   ...               ...               ... ...      ...     ...       ...
       1 Ceres 2018-Aug-04 06:40 2458334.777777778 ... 171.5076  9.3369       UTC
       1 Ceres 2018-Aug-04 06:50 2458334.784722222 ... 171.5099  9.3368       UTC
       1 Ceres 2018-Aug-04 07:00 2458334.791666667 ... 171.5123  9.3367       UTC
       1 Ceres 2018-Aug-04 07:10 2458334.798611111 ... 171.5147  9.3366       UTC
       1 Ceres 2018-Aug-04 07:20 2458334.805555556 ... 171.5171  9.3365       UTC
       1 Ceres 2018-Aug-04 07:30      2458334.8125 ... 171.5195  9.3364       UTC

Note that ``skip_daylight`` is an optional parameter of
`~astroquery.jplhorizons.HorizonsClass.ephemerides` and it can be used
here as well. An additional feature of
`~sbpy.data.Ephem.from_horizons` is that you can automatically
concatenate queries for a number of objects:

    >>> eph = Ephem.from_horizons(['Ceres', 'Pallas', 12893, '1983 SA'],
    ...                           location='568',
    ...                           epochs=epoch1)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable masked=True length=4>
	    targetname               datetime_str       ...  PABLat  timescale
							...   deg
	      str26                     str24           ... float64     str3
    -------------------------- ------------------------ ... -------- ---------
		       1 Ceres 2018-Aug-03 14:20:00.000 ...   9.3473       UTC
		      2 Pallas 2018-Aug-03 14:20:00.000 ... -20.1396       UTC
     12893 Mommert (1998 QS55) 2018-Aug-03 14:20:00.000 ...  -2.0567       UTC
    3552 Don Quixote (1983 SA) 2018-Aug-03 14:20:00.000 ...  13.3365       UTC

Please be aware that these queries are not simultaneous. The more
targets you query, the longer the query will take. Furthermore, keep
in mind that asteroids and comets have slightly different table
layouts (e.g., different magnitude systems: ``T-mag`` and ``N-mag``
instead of ``V-mag``), which will complicate the interpretation of the
data. It might be safest to query asteroids and comets separately.

Also note that the two examples shown above use different ways to
define epochs. The first example uses a dictionary that defines a
``start`` and ``stop`` epoch, as well as a ``step`` size (see
`~astroquery.jplhorizons.HorizonsClass.ephemerides` for
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

Observer locations can be defined as strings using offical `IAU
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
    targetname       datetime_str       datetime_jd ...  PABLon   PABLat timescale
					     d      ...   deg      deg            
       str7             str24             float64   ... float64  float64    str3  
    ---------- ------------------------ ----------- ... -------- ------- ---------
    1 Ceres 2018-Jan-01 00:00:00.000   2458119.5 ... 130.4303  9.2004       UTC

Offering almost identical functionality, the
`~sbpy.data.Ephem.from_mpc` method will retrieve ephemerides from the
`Minor Planet Center <https://minorplanetcenter.net/>`_:

    >>> eph = Ephem.from_mpc('2P', location='568',
    ...                      epochs={'start': '2018-10-22',
    ...                              'stop': '2018-10-26',
    ...                              'step': '1d'})  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable length=5>
	      Date          timescale ... Moon distance Moon altitude
				      ...      deg           deg
	     object            str3   ...    float64       float64
    ----------------------- --------- ... ------------- -------------
    2018-10-22 00:00:00.000       UTC ...          28.0         -33.0
    2018-10-23 00:00:00.000       UTC ...          41.0         -41.0
    2018-10-24 00:00:00.000       UTC ...          54.0         -48.0
    2018-10-25 00:00:00.000       UTC ...          67.0         -53.0
    2018-10-26 00:00:00.000       UTC ...          81.0         -56.0

Finally, `~sbpy.data.Ephem.from_miriade` will retrieve ephemerides
from the `Miriade ephemeris generator
<http://vo.imcce.fr/webservices/miriade/>`_ at `IMCCE
<https://www.imcce.fr/>`_:

    >>> eph = Ephem.from_miriade('2P', objtype='comet', location='568',
    ...                          epochs={'start': '2018-10-22',
    ...                                  'stop': '2018-10-26',
    ...                                  'step': '1d'})  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable masked=True length=5>
     target        epoch                 RA         ...  delta_rate  timescale
		     d                  deg         ...    km / s             
    bytes20       float64             float64       ...   float64       str3  
    ------- -------------------- ------------------ ... ------------ ---------
	 2P            2458413.5 329.99213124999994 ...   24.7933113       UTC
	 2P            2458414.5 329.91132124999996 ...   25.0280603       UTC
	 2P            2458415.5 329.83517041666664 ...    25.253586       UTC
	 2P            2458416.5 329.76366666666667 ...   25.4700287       UTC
	 2P            2458417.5  329.6967958333333 ...    25.677518       UTC    

    
Ephemerides can also be derived from `~Orbit` objects using `sbpy`'s
interface to `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_ with the function
`~sbpy.data.Ephem.from_oorb`. The following example computes
ephemerides for the next ten days in steps of 1 hr for Ceres as seen
from the Discovery Channel Telescope:

     >>> import numpy as np
     >>> from sbpy.data import Orbit, Ephem
     >>> from astropy.time import Time
     >>> epochs = Time.now().jd + np.arange(0, 10, 1/24)
     >>> ceres = Orbit.from_horizons('1')  # doctest: +REMOTE_DATA
     >>> eph = Ephem.from_oo(ceres, epochs, 'G37') # doctest: +SKIP 
     >>> print(eph) # doctest: +SKIP 
     <QTable length=240>
     targetname       epoch        ...           obsz               trueanom
			d          ...            AU                  deg
	str7         float64       ...         float64              float64
     ---------- ------------------ ... ----------------------- -----------------
	1 Ceres 2458519.2878717002 ...   4.886414464166933e-06 68.07980642088688
	1 Ceres 2458519.3295383668 ...  2.3814767035612583e-06  68.0893160393968
	1 Ceres 2458519.3712050337 ...  -7.136200919632962e-07 68.09882544202566
	1 Ceres 2458519.4128717002 ...   -4.18340743346679e-06 68.10833462855386
	1 Ceres 2458519.4545383668 ...  -7.786747377891423e-06 68.11784359908062
	1 Ceres 2458519.4962050337 ... -1.1273355301266719e-05 68.12735235370518
	    ...                ... ...                     ...               ...
	1 Ceres 2458529.0378717002 ...   1.093565783852335e-05 70.29915515170745
	1 Ceres 2458529.0795383668 ...  1.3089531693877277e-05  70.3086140523456
	1 Ceres 2458529.1212050337 ...  1.4402894355114437e-05 70.31807273565124
	1 Ceres 2458529.1628717002 ...  1.4786143903738891e-05 70.32753120140761
	1 Ceres 2458529.2045383668 ...  1.4213398342149963e-05 70.33698944971509
	1 Ceres 2458529.2462050337 ...  1.2724269065650384e-05 70.34644748067402

The properties computed by pyoorb and listed in the resulting table
are defined in the `pyoorb documentation
<https://github.com/oorb/oorb/tree/master/python>`_. Note that this function requires pyoorb to be installed, which is not a requirement for `sbpy`.

`~sbpy.data.Obs` works exactly like `~sbpy.data.Ephem`, but this class
will feature in the future some convenience functions to be able to
better deal with observational data.


Obs Functionality
=================

`~sbpy.data.Obs` objects have the same functionality as
`~sbpy.data.Ephem` as well as some unique functions.

For instance, this class allows you to query observations reported to
the Minor Planet Center for a given target:

    >>> from sbpy.data import Obs
    >>> data = Obs.from_mpc('2019 AA', id_type='asteroid designation') # doctest: +REMOTE_DATA
    >>> data
    <QTable masked=True length=33>
    number  desig  discovery note1 ...        DEC           mag   band observatory
				   ...        deg           mag
    int64    str7     str1    str1 ...      float64       float64 str1     str3
    ------ ------- --------- ----- ... ------------------ ------- ---- -----------
	-- 2019 AA        --    -- ...  42.32416944444445    20.2    G         F51
	-- 2019 AA        --    -- ...  42.32879722222223    20.3    G         F51
	-- 2019 AA        --    -- ... 42.333225000000006    20.3    G         F51
	-- 2019 AA         *    -- ...  46.52321666666666    20.0    w         F51
	-- 2019 AA        --    -- ...  46.52748611111111    20.0    w         F51
	-- 2019 AA        --    -- ... 46.531755555555556    20.0    w         F51
       ...     ...       ...   ... ...                ...     ...  ...         ...
	-- 2019 AA        --    -- ... 46.706500000000005    20.2    V         033
	-- 2019 AA        --    -- ...  46.70652777777778    20.2    V         033
	-- 2019 AA        --    -- ...  49.73566111111111    20.1    i         F52
	-- 2019 AA        --    -- ...  49.73788888888889    20.1    i         F52
	-- 2019 AA        --    -- ...  49.74008611111111    20.1    i         F52
	-- 2019 AA        --    -- ... 49.742266666666666    20.2    i         F52


For a given `~sbpy.data.Obs` object, `~sbpy.data.Obs.supplement`
allows you to supplement the information content of this object by
adding ephemeris data for the target(s) and epochs provided. This
function makes use of the query functions that are part of
`~sbpy.data.Ephem` and allows you to pick a service from which you
would like to obtain the data.

    >>> data.field_names
    <TableColumns names=('number','desig','discovery','note1','note2','epoch','RA','DEC','mag','band','observatory')>
    >>> data_sup = data.supplement()
    <TableColumns names=('number','desig','discovery','note1','note2','epoch','RA_obs','DEC_obs','mag','band','observatory','targetname','datetime_str','H','G','solar_presence','flags','RA','DEC','RA_app','DEC_app','RA*cos(Dec)_rate','DEC_rate','AZ','EL','AZ_rate','EL_rate','sat_X','sat_Y','sat_PANG','siderealtime','airmass','magextinct','V','illumination','illum_defect','sat_sep','sat_vis','ang_width','PDObsLon','PDObsLat','PDSunLon','PDSunLat','SubSol_ang','SubSol_dist','NPole_ang','NPole_dist','EclLon','EclLat','r','r_rate','delta','delta_rate','lighttime','vel_sun','vel_obs','elong','elongFlag','alpha','lunar_elong','lunar_illum','sat_alpha','sunTargetPA','velocityPA','OrbPlaneAng','constellation','TDB-UT','ObsEclLon','ObsEclLat','NPole_RA','NPole_DEC','GlxLon','GlxLat','solartime','earth_lighttime','RA_3sigma','DEC_3sigma','SMAA_3sigma','SMIA_3sigma','Theta_3sigma','Area_3sigma','RSS_3sigma','r_3sigma','r_rate_3sigma','SBand_3sigma','XBand_3sigma','DoppDelay_3sigma','true_anom','hour_angle','alpha_true','PABLon','PABLat')>
    >>> data_sup['r']  # doctest: +SKIP
    [1.87637455 1.87638696 1.87639935 1.8891401  1.88915504 1.88916999
 1.88918492 1.88963885 1.88964435 1.88965078 1.8896515  1.88965233
 1.88965306 1.88970198 1.88970271 1.88970353 1.88970426 1.88970508
 1.88970581 1.88978031 1.88978103 1.88978186 1.8897826  1.88978341
 1.88978415 1.88978497 1.88978569 1.88978653 1.88978725 1.90342431
 1.903438   1.90345164 1.90346529] AU

`~sbpy.data.Obs.supplement` queries in this case ephemerides from the
JPL Horizons system (the default ``service`` to be used) and appends
the additional data as new columns to the `~sbpy.data.Obs`
object. Targetnames and epochs are taken from the underlying
`~sbpy.data.Obs` object and must hence be present there; keyword
arguments ``id_field`` and ``epoch_field`` control the fieldnames from
which these information are taken. In order to prevent duplicate
column names, the keyword argument ``modify_fieldnames`` defines
whether field names in the `~sbpy.data.Obs` object or in the new
`~sbpy.data.Ephem` data are to be modified.

Note that services using `~sbpy.data.Ephem.from_mpc` and
`~sbpy.data.Ephem.from_miriade` are not optimized for queries of
multiple epochs. Hence, using these service will result in long
loading times. `~sbpy.data.Ephem.from_horizons`, however, is optimized
for this type of query.


