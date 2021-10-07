===========
 Using Obs
===========

`~sbpy.data.Obs` objects mostly share their functionality with
`~sbpy.data.Ephem`, but there are some unique features tailored to observational data.

For instance, this class allows you to query observations reported to
the Minor Planet Center for a given target:

    >>> from sbpy.data import Obs
    >>> data = Obs.from_mpc('2019 AA', id_type='asteroid designation') # doctest: +REMOTE_DATA
    >>> data # doctest: +REMOTE_DATA +SKIP
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

    >>> data.field_names # doctest: +REMOTE_DATA
    ['number', 'desig', 'discovery', 'note1', 'note2', 'epoch', 'RA', 'DEC', 'mag', 'band', 'observatory']
    >>> data_sup = data.supplement(id_field='desig') # doctest: +REMOTE_DATA
    >>> data_sup.field_names # doctest: +REMOTE_DATA
    ['number', 'desig', 'discovery', 'note1', 'note2', 'epoch', 'RA_obs', 'DEC_obs', 'mag', 'band', 'observatory', 'targetname', 'H', 'G', 'solar_presence', 'flags', 'RA', 'DEC', 'RA_app', 'DEC_app', 'RA*cos(Dec)_rate', 'DEC_rate', 'AZ', 'EL', 'AZ_rate', 'EL_rate', 'sat_X', 'sat_Y', 'sat_PANG', 'siderealtime', 'airmass', 'magextinct', 'V', 'surfbright', 'illumination', 'illum_defect', 'sat_sep', 'sat_vis', 'ang_width', 'PDObsLon', 'PDObsLat', 'PDSunLon', 'PDSunLat', 'SubSol_ang', 'SubSol_dist', 'NPole_ang', 'NPole_dist', 'EclLon', 'EclLat', 'r', 'r_rate', 'delta', 'delta_rate', 'lighttime', 'vel_sun', 'vel_obs', 'elong', 'elongFlag', 'alpha', 'lunar_elong', 'lunar_illum', 'sat_alpha', 'sunTargetPA', 'velocityPA', 'OrbPlaneAng', 'constellation', 'TDB-UT', 'ObsEclLon', 'ObsEclLat', 'NPole_RA', 'NPole_DEC', 'GlxLon', 'GlxLat', 'solartime', 'earth_lighttime', 'RA_3sigma', 'DEC_3sigma', 'SMAA_3sigma', 'SMIA_3sigma', 'Theta_3sigma', 'Area_3sigma', 'RSS_3sigma', 'r_3sigma', 'r_rate_3sigma', 'SBand_3sigma', 'XBand_3sigma', 'DoppDelay_3sigma', 'true_anom', 'hour_angle', 'alpha_true', 'PABLon', 'PABLat']
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


