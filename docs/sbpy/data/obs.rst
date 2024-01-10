============================================
Observational Data Objects (`sbpy.data.Obs`)
============================================

`~sbpy.data.Obs` objects mostly share their functionality with
`~sbpy.data.Ephem`, but there are some unique features tailored to observational data.

For instance, this class allows you to query observations reported to the Minor
Planet Center for a given target via `astroquery.mpc.MPCClass.get_observations`:

... .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> from sbpy.data import Obs
    >>> data = Obs.from_mpc('2019 AA', id_type='asteroid designation')
    >>> data = data[:10]  # limit the number of rows for this example
    >>> data # doctest: +SKIP
    <QTable length=10>
    number  desig  discovery note1 ...         DEC           mag   band observatory
                                   ...         deg           mag                   
    int64    str7     str1    str1 ...       float64       float64 str1     str3   
    ------ ------- --------- ----- ... ------------------- ------- ---- -----------
        -- 2019 AA        --    -- ... -2.6917250000000004    22.2    w         F51
        -- 2019 AA        --    -- ... -2.6893833333333337    22.1    w         F51
        -- 2019 AA        --    -- ... -2.6872222222222226    22.2    w         F51
        -- 2019 AA        --    -- ...   18.10141666666667    20.3    V         703
        -- 2019 AA        --    -- ...  18.100916666666667    19.4    V         703
        -- 2019 AA        --    -- ...  18.100555555555555    20.2    V         703
        -- 2019 AA        --    -- ...   18.10086111111111    20.8    V         703
        -- 2019 AA        --    -- ...   22.21022222222222    19.6    V         G96
        -- 2019 AA        --    -- ...  22.210083333333333    19.8    V         G96
        -- 2019 AA        --    -- ...   22.21022222222222    20.2    V         G96

For a given `~sbpy.data.Obs` object, `~sbpy.data.Obs.supplement`
allows you to supplement the information content of this object by
adding ephemeris data for the target(s) and epochs provided. This
function makes use of the query functions that are part of
`~sbpy.data.Ephem` and allows you to pick a service from which you
would like to obtain the data.

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> data.field_names
    ['number', 'desig', 'discovery', 'note1', 'note2', 'epoch', 'RA', 'DEC', 'mag', 'band', 'observatory']
    >>> data_sup = data.supplement(id_field='desig')
    >>> data_sup.field_names
    ['number', 'desig', 'discovery', 'note1', 'note2', 'epoch', 'RA_obs', 'DEC_obs', 'mag', 'band', 'observatory', 'targetname', 'H', 'G', 'solar_presence', 'flags', 'RA', 'DEC', 'RA_app', 'DEC_app', 'RA*cos(Dec)_rate', 'DEC_rate', 'AZ', 'EL', 'AZ_rate', 'EL_rate', 'sat_X', 'sat_Y', 'sat_PANG', 'siderealtime', 'airmass', 'magextinct', 'V', 'surfbright', 'illumination', 'illum_defect', 'sat_sep', 'sat_vis', 'ang_width', 'PDObsLon', 'PDObsLat', 'PDSunLon', 'PDSunLat', 'SubSol_ang', 'SubSol_dist', 'NPole_ang', 'NPole_dist', 'EclLon', 'EclLat', 'r', 'r_rate', 'delta', 'delta_rate', 'lighttime', 'vel_sun', 'vel_obs', 'elong', 'elongFlag', 'alpha', 'lunar_elong', 'lunar_illum', 'sat_alpha', 'sunTargetPA', 'velocityPA', 'OrbPlaneAng', 'constellation', 'TDB-UT', 'ObsEclLon', 'ObsEclLat', 'NPole_RA', 'NPole_DEC', 'GlxLon', 'GlxLat', 'solartime', 'earth_lighttime', 'RA_3sigma', 'DEC_3sigma', 'SMAA_3sigma', 'SMIA_3sigma', 'Theta_3sigma', 'Area_3sigma', 'RSS_3sigma', 'r_3sigma', 'r_rate_3sigma', 'SBand_3sigma', 'XBand_3sigma', 'DoppDelay_3sigma', 'true_anom', 'hour_angle', 'alpha_true', 'PABLon', 'PABLat']
    >>> data_sup['r']  # doctest: +SKIP
    <MaskedQuantity [1.74110682, 1.74111208, 1.74111672, 1.97466336,
                     1.97465973, 1.97465603, 1.97465236, 1.98676694,
                     1.98676467, 1.9867624 ] AU>

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


