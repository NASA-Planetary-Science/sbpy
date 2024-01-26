# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
sbpy.data
---------

:author: Michael Mommert (mommermiscience@gmail.com)

``Conf`` contains metadata for ``sbpy`` `~sbpy.DataClass` field names.

``Conf.fieldnames_info`` is a list of dictionaries, one per field, with the
keys: `'description'`, `'fieldnames'`, `'provenance'`, `'dimension'`, and
`'equivalencies'`:

    * description: text description of the field
    * provenance: list of `~sbpy.DataClass` objects which use the field
    * fieldnames: list of field names as strings, the first is considered the
      primary, the remaining strings, if any, are alternates
    * dimension: English description of the field's dimension (e.g., length)
    * equivalencies: `~astropy.units` list of equivalencies for unit
      conversion (optional)

"""

__all__ = [
    "DataClass",
    "Ephem",
    "Obs",
    "Orbit",
    "Phys",
    "Names",
    "Conf",
    "DataClassError",
    "quantity_to_dataclass",
    "natural_sort_key",
    "dataclass_input",
    "QueryError",
    "TimeScaleWarning",
]

import astropy.units as u
from astropy.time import Time
from . import dimensions


class Conf:
    # acceptable field names for DataClass
    fieldnames_info = [
        # General
        {
            "description": "Target Identifier",
            "fieldnames": ["targetname", "id", "Object"],
            "provenance": ["orbit", "ephem", "obs", "phys"],
            "dimension": None,
        },
        {
            "description": "Target Designation",
            "fieldnames": ["desig", "designation"],
            "provenance": ["orbit", "ephem", "obs", "phys"],
            "dimension": None,
        },
        {
            "description": "Target Number",
            "fieldnames": ["number"],
            "provenance": ["orbit", "ephem", "obs", "phys"],
            "dimension": None,
        },
        {
            "description": "Target Name",
            "fieldnames": ["name"],
            "provenance": ["orbit", "ephem", "obs", "phys"],
            "dimension": None,
        },
        {
            "description": "Epoch",
            "fieldnames": ["epoch", "datetime", "Date", "date", "Time", "time"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.time_object,
        },
        {
            "description": "relative time",
            "fieldnames": ["t_relative", "t_rel", "dt"],
            "provenance": ["ephem"],
            "dimension": dimensions.time,
        },
        # Orbital Elements
        {
            "description": "Semi-Major Axis",
            "fieldnames": ["a", "sma"],
            "provenance": ["orbit"],
            "dimension": dimensions.length,
        },
        {
            "description": "Eccentricity",
            "fieldnames": ["e", "ecc"],
            "provenance": ["orbit"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Inclination",
            "fieldnames": ["i", "inc", "incl"],
            "provenance": ["orbit"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Perihelion Distance",
            "fieldnames": ["q", "periheldist"],
            "provenance": ["orbit"],
            "dimension": dimensions.length,
        },
        {
            "description": "Aphelion Distance",
            "fieldnames": ["Q", "apheldist"],
            "provenance": ["orbit"],
            "dimension": dimensions.length,
        },
        {
            "description": "Longitude of the Ascending Node",
            "fieldnames": ["Omega", "longnode", "node"],
            "provenance": ["orbit"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Argument of the Periapsis",
            "fieldnames": ["w", "argper"],
            "provenance": ["orbit"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Mean Anomaly",
            "fieldnames": ["M", "mean_anom"],
            "provenance": ["orbit"],
            "dimension": dimensions.angle,
        },
        {
            "description": "True Anomaly",
            "fieldnames": ["v", "true_anom", "true_anomaly"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Arc Length",
            "fieldnames": ["arc", "arc_length"],
            "provenance": ["orbit", "ephem"],
            "dimension": dimensions.time,
        },
        {
            "description": "Delta-v",
            "fieldnames": ["delta_v", "delta-v"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Mercury",
            "fieldnames": ["moid_mercury"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Earth",
            "fieldnames": ["moid_earth"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Venus",
            "fieldnames": ["moid_venus"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Mars",
            "fieldnames": ["moid_mars"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Jupiter",
            "fieldnames": ["moid_jupiter"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Saturn",
            "fieldnames": ["moid_saturn"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Uranus",
            "fieldnames": ["moid_uranus"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Minimum Orbit Intersection Distance wrt Neptune",
            "fieldnames": ["moid_neptune"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Tisserand Parameter wrt Jupiter",
            "fieldnames": ["Tj", "tj"],
            "provenance": ["orbit", "phys"],
            "dimension": None,
        },
        {
            "description": "MPC Orbit Type",
            "fieldnames": ["mpc_orb_type"],
            "provenance": ["orbit", "phys"],
            "dimension": None,
        },
        {
            "description": "Epoch of Perihelion Passage",
            "fieldnames": ["Tp"],
            "provenance": ["orbit"],
            "dimension": dimensions.time_object,
        },
        {
            "description": "Orbital Period",
            "fieldnames": ["P", "period"],
            "provenance": ["orbit", "phys"],
            "dimension": dimensions.time,
        },
        # Ephemerides properties
        {
            "description": "Heliocentric Distance",
            "fieldnames": ["r", "rh", "r_hel", "heldist"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "Heliocentric Radial Velocity",
            "fieldnames": ["r_rate", "rh_rate", "rdot", "r-dot", "rhdot", "rh-dot"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "Distance to the Observer",
            "fieldnames": ["delta", "Delta", "obsdist"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "Observer-Target Radial Velocity",
            "fieldnames": ["delta_rate", "deltadot", "delta-dot", "deldot", "del-dot"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "Right Ascension",
            "fieldnames": ["ra", "RA"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Declination",
            "fieldnames": ["dec", "DEC", "Dec"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Right Ascension Rate",
            "fieldnames": ["ra_rate", "RA_rate", "ra_rates", "RA_rates", "dRA", "dra"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle_per_time,
        },
        {
            "description": "RA*cos(Dec) Rate",
            "fieldnames": [
                "RA*cos(Dec)_rate",
                "dra cos(dec)",
                "dRA cos(Dec)",
                "dra*cos(dec)",
                "dRA*cos(Dec)",
            ],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle_per_time,
        },
        {
            "description": "Declination Rate",
            "fieldnames": [
                "dec_rate",
                "DEC_rate",
                "Dec_rate",
                "dec_rates",
                "DEC_rates",
                "Dec_rates",
                "dDec",
                "dDEC",
                "ddec",
            ],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle_per_time,
        },
        {
            "description": "Proper Motion",
            "fieldnames": ["mu", "Proper motion"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle_per_time,
        },
        {
            "description": "Proper Motion Direction",
            "fieldnames": ["Direction", "direction"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Solar Phase Angle",
            "fieldnames": ["alpha", "phaseangle", "Phase", "phase"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Solar Elongation Angle",
            "fieldnames": [
                "elong",
                "solarelong",
                "solarelongation",
                "elongation",
                "Elongation",
            ],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "V-band Magnitude",
            "fieldnames": ["V", "Vmag"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.magnitude,
        },
        {
            "description": "Heliocentric Ecliptic Longitude",
            "fieldnames": ["hlon", "EclLon", "ecllon", "HelEclLon", "helecllon"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Heliocentric Ecliptic Latitude",
            "fieldnames": ["hlat", "EclLat", "ecllat", "HelEclLat", "helecllat"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Horizontal Elevation",
            "fieldnames": ["el", "EL", "elevation", "alt", "altitude", "Altitude"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Horizontal Azimuth",
            "fieldnames": ["az", "AZ", "azimuth"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Lunar Elongation",
            "fieldnames": [
                "lunar_elong",
                "elong_moon",
                "elongation_moon",
                "lunar_elongation",
                "lunarelong",
            ],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "X State Vector Component",
            "fieldnames": ["x", "X", "x_vec"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "Y State Vector Component",
            "fieldnames": ["y", "Y", "y_vec"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "Z State Vector Component",
            "fieldnames": ["z", "Z", "z_vec"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "X Velocity Vector Component",
            "fieldnames": ["vx", "dx", "dx/dt"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "Y Velocity Vector Component",
            "fieldnames": ["vy", "dy", "dy/dt"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "Z Velocity Vector Component",
            "fieldnames": ["vz", "dz", "dz/dt"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "X heliocentric position vector",
            "fieldnames": ["x_h", "X_h"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "Y heliocentric position vector",
            "fieldnames": ["y_h", "Y_h"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "Z heliocentric position vector",
            "fieldnames": ["z_h", "Z_h"],
            "provenance": ["orbit", "ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "Comet Total Absolute Magnitude",
            "fieldnames": ["m1", "M1"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.magnitude,
        },
        {
            "description": "Comet Nuclear Absolute Magnitude",
            "fieldnames": ["m2", "M2"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.magnitude,
        },
        {
            "description": "Total Magnitude Scaling Factor",
            "fieldnames": ["k1", "K1"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Nuclear Magnitude Scaling Factor",
            "fieldnames": ["k2", "K2"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Phase Coefficient",
            "fieldnames": ["phase_coeff", "Phase_coeff"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Information on Solar Presence",
            "fieldnames": ["solar_presence", "Solar_presence"],
            "provenance": ["ephem", "obs"],
            "dimension": None,
        },
        {
            "description": "Information on Moon and target status",
            "fieldnames": ["status_flag", "Status_flag"],
            "provenance": ["ephem", "obs"],
            "dimension": None,
        },
        {
            "description": "Apparent Right Ascension",
            "fieldnames": ["RA_app", "ra_app"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Apparent Declination",
            "fieldnames": ["DEC_app", "dec_app"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Azimuth Rate (dAZ*cosE)",
            "fieldnames": ["az_rate", "AZ_rate"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle_per_time,
        },
        {
            "description": "Elevation Rate (d(ELV)/dt)",
            "fieldnames": ["el_rate", "EL_rate"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle_per_time,
        },
        {
            "description": "Satellite Position Angle",
            "fieldnames": ["sat_pang", "Sat_pang"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Local Sidereal Time",
            "fieldnames": ["siderealtime", "Siderealtime"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.time,
        },
        {
            "description": "Target Optical Airmass",
            "fieldnames": ["airmass", "Airmass"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "V Magnitude Extinction",
            "fieldnames": ["vmagex", "Vmagex"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.magnitude,
        },
        {
            "description": "Surface Brightness",
            "fieldnames": ["Surfbright", "surfbright"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.magnitude_per_solid_angle,
        },
        {
            "description": "Fraction of Illumination",
            "fieldnames": ["frac_illum", "Frac_illum"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.percent,
        },
        {
            "description": "Illumination Defect",
            "fieldnames": ["defect_illum", "Defect_illum"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Target-primary angular separation",
            "fieldnames": ["targ_sep", "Targ_sep"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Target-primary visibility",
            "fieldnames": ["targ_vis", "Targ_vis"],
            "provenance": ["ephem", "obs"],
            "dimension": None,
        },
        {
            "description": "Angular width of target",
            "fieldnames": ["targ_width", "Targ_width"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Apparent planetodetic longitude",
            "fieldnames": ["pldetic_long", "Pldetic_long"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Apparent planetodetic latitude",
            "fieldnames": ["pldetic_lat", "Pldetic_lat"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Apparent planetodetic Solar longitude",
            "fieldnames": ["pltdeticSol_long", "PltdeticSol_long"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Apparent planetodetic Solar latitude",
            "fieldnames": ["pltdeticSol_lat", "PltdeticSol_lat"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Target sub-solar point position angle",
            "fieldnames": ["subsol_ang", "Subsol_ang"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Target sub-solar point angle distance",
            "fieldnames": ["subsol_dist", "Subsol_dist"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Target North pole position angle",
            "fieldnames": ["npole_angle", "Npole_angle"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Target North pole position distance",
            "fieldnames": ["npole_dist", "Npole_dist"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Observation centric ecliptic longitude",
            "fieldnames": ["obs_ecl_long", "Obs_ecl_long"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Observation centric ecliptic latitude",
            "fieldnames": ["obs_ecl_lat", "Obs_ecl_lat"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "One-way light time",
            "fieldnames": ["lighttime", "Lighttime"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.time,
        },
        {
            "description": "Target center velocity wrt Sun",
            "fieldnames": ["vel_sun", "Vel_sun"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "Target center velocity wrt Observer",
            "fieldnames": ["vel_obs", "Vel_obs"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "Lunar illumination",
            "fieldnames": ["lun_illum", "Lun_illum"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.percent,
        },
        {
            "description": "Apparent interfering body elongation wrt observer",
            "fieldnames": ["ib_elong", "IB_elong"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Interfering body illumination",
            "fieldnames": ["ib_illum", "IB_illum"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.percent,
        },
        {
            "description": "Observer primary target angle",
            "fieldnames": ["targ_angle_obs", "Targ_angle_obs"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Orbital plane angle",
            "fieldnames": ["orbangle_plane", "Orbangle_plane"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Constellation ID containing target",
            "fieldnames": ["constellation", "Constellation"],
            "provenance": ["ephem", "obs"],
            "dimension": None,
        },
        {
            "description": "Target North Pole RA",
            "fieldnames": ["targ_npole_ra", "targ_npole_RA"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Target North Pole DEC",
            "fieldnames": ["targ_npole_dec", "targ_npole_DEC"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Galactic Longitude",
            "fieldnames": ["glx_long", "Glx_long"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Galactic Latitude",
            "fieldnames": ["glx_lat", "Glx_lat"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Local apparent solar time",
            "fieldnames": ["solartime"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.time,
        },
        {
            "description": "Observer light time from Earth",
            "fieldnames": ["earthlighttime", "Earthlighttime"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.time,
        },
        {
            "description": "3 sigma positional uncertainty RA",
            "fieldnames": ["RA_3sigma", "ra_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "3 sigma positional uncertainty DEC",
            "fieldnames": ["DEC_3sigma", "dec_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "3 sigma positional uncertainty semi-major axis",
            "fieldnames": ["sma_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "3 sigma positional uncertainty semi-minor axis",
            "fieldnames": ["smi_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "3 sigma positional uncertainty position angle",
            "fieldnames": ["posangle_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "3 sigma positional uncertainty ellipse area",
            "fieldnames": ["area_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.solid_angle,
        },
        {
            "description": "3 sigma positional uncertainty root sum square",
            "fieldnames": ["rss_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "3 sigma range uncertainty",
            "fieldnames": ["r_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length,
        },
        {
            "description": "3 sigma range rate uncertainty",
            "fieldnames": ["r_rate_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.length_per_time,
        },
        {
            "description": "3 sigma doppler radar uncertainty at S-band",
            "fieldnames": ["sband_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.frequency,
        },
        {
            "description": "3 sigma doppler radar uncertainty at X-band",
            "fieldnames": ["xband_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.frequency,
        },
        {
            "description": "3 sigma doppler round-trip delay uncertainty",
            "fieldnames": ["dopdelay_3sigma"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.time,
        },
        {
            "description": "Local apparent hour angle",
            "fieldnames": ["locapp_hourangle"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.time,
        },
        {
            "description": "True phase angle",
            "fieldnames": ["true_phaseangle"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Phase angle bisector longitude",
            "fieldnames": ["pab_long"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Phase angle bisector latitude",
            "fieldnames": ["pab_lat"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Absolute V-band Magnitude",
            "fieldnames": ["abs_V", "abs_Vmag"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.magnitude,
        },
        {
            "description": "Satellite X-position",
            "fieldnames": ["sat_X", "sat_x"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Satellite Y-position",
            "fieldnames": ["sat_y", "sat_Y"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        {
            "description": "Atmospheric Refraction",
            "fieldnames": ["atm_refraction", "refraction"],
            "provenance": ["ephem", "obs"],
            "dimension": dimensions.angle,
        },
        # Physical properties (dependent on other properties)
        {
            "description": "Infrared Beaming Parameter",
            "fieldnames": ["eta", "Eta"],
            "provenance": ["ephem", "obs"],
            "dimension": None,
        },
        {
            "description": "Temperature",
            "fieldnames": ["temp", "Temp", "temperature", "Temperature"],
            "provenance": ["phys", "ephem", "obs"],
            "dimension": dimensions.temperature,
        },
        # Physical properties (static)
        {
            "description": "Effective Diameter",
            "fieldnames": ["d", "D", "diam", "diameter", "Diameter"],
            "provenance": ["phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Effective Radius",
            "fieldnames": ["R", "radius"],
            "provenance": ["phys"],
            "dimension": dimensions.length,
        },
        {
            "description": "Geometric Albedo",
            "fieldnames": ["pv", "pV", "p_v", "p_V", "geomalb"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Bond Albedo",
            "fieldnames": ["A", "bondalbedo"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Emissivity",
            "fieldnames": ["emissivity", "Emissivity"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Absolute Magnitude",
            "fieldnames": ["absmag", "H"],
            "provenance": ["phys", "ephem", "orbit"],
            "dimension": dimensions.magnitude,
        },
        {
            "description": "Photometric Phase Slope Parameter",
            "fieldnames": ["G", "slope"],
            "provenance": ["phys", "ephem", "orbit"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Molecule Identifier",
            "fieldnames": ["mol_tag", "mol_name"],
            "provenance": ["phys"],
            "dimension": None,
        },
        {
            "description": "Transition frequency",
            "fieldnames": ["t_freq"],
            "provenance": ["phys"],
            "dimension": dimensions.frequency,
            "equivalencies": u.spectral(),
        },
        {
            "description": "Integrated line intensity at 300 K",
            "fieldnames": ["lgint300"],
            "provenance": ["phys"],
            "dimension": None,
        },  # fix when intensity units are resolved
        {
            "description": "Integrated line intensity at designated Temperature",
            "fieldnames": ["intl", "lgint"],
            "provenance": ["phys"],
            "dimension": None,
        },  # fix when intensity units are resolved
        {
            "description": "Partition function at 300 K",
            "fieldnames": ["partfn300"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Partition function at designated temperature",
            "fieldnames": ["partfn"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Upper state degeneracy",
            "fieldnames": ["dgup"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Upper level energy in Joules",
            "fieldnames": ["eup_j", "eup_J"],
            "provenance": ["phys"],
            "dimension": dimensions.energy,
        },
        {
            "description": "Lower level energy in Joules",
            "fieldnames": ["elo_j", "elo_J"],
            "provenance": ["phys"],
            "dimension": dimensions.energy,
        },
        {
            "description": "Degrees of freedom",
            "fieldnames": ["degfr", "ndf", "degfreedom"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Einstein Coefficient",
            "fieldnames": ["au", "eincoeff"],
            "provenance": ["phys"],
            "dimension": dimensions.inverse_time,
        },
        {
            "description": "Timescale * r^2",
            "fieldnames": ["beta", "beta_factor"],
            "provenance": ["phys"],
            "dimension": dimensions.time_area,
        },
        {
            "description": "Total Number",
            "fieldnames": ["totnum", "total_number_nocd" "total_number"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        {
            "description": "Column Density from Bockelee Morvan et al. 2004",
            "fieldnames": ["cdensity", "col_density"],
            "provenance": ["phys"],
            "dimension": dimensions.inverse_area,
        },
        {
            "description": "Ratio of the force from radiation to the force from gravity",
            "fieldnames": ["beta_rad"],
            "provenance": ["phys"],
            "dimension": dimensions.dimensionless,
        },
        # {  # see module doc string
        #   'description': '',
        #   'fieldnames': [],
        #   'provenance': [],
        #   'dimension': None,
        #   'equivalencies': (astropy units equivalencies, e.g., u.spectral())
        # },
    ]

    # use this code snippet to identify duplicate field names:
    # from sbpy.data import Conf
    # import collections
    # a = sum(Conf.fieldnames, [])
    # print([item for item, count in collections.Counter(a).items()
    #        if count > 1])

    # list of fieldnames; each element a list of alternatives
    fieldnames = [prop["fieldnames"] for prop in fieldnames_info]

    fieldname_idx = {}
    for idx, field in enumerate(fieldnames):
        for alt in field:
            fieldname_idx[alt] = idx

    # field equivalencies defining conversions
    # key defines target quantity; dict with source quantity and function
    # for conversion
    # conversions considered as part of DataClass._translate_columns
    field_eq = {
        "R": {"d": lambda r: r / 2},
        # diameter to radius}
        "d": {"R": lambda d: d * 2},
    }

    # definitions for use of pyoorb in Orbits
    oorb_timeScales = {"UTC": 1, "UT1": 2, "TT": 3, "TAI": 4}
    oorb_elemType = {"CART": 1, "COM": 2, "KEP": 3, "DEL": 4, "EQX": 5}

    # field name, unit; in order as returned from oorb
    # However, in propagate, angular units are returned as deg.  This is
    # accounted for in Orbit.oo_propagate().
    oorb_orbit_fields = {
        "COM": (
            ("id", None),
            ("q", "au"),
            ("e", ""),
            ("incl", "rad"),
            ("Omega", "rad"),
            ("w", "rad"),
            ("Tp", "d"),
            ("orbtype", None),
            ("epoch", "d"),
            ("epoch_scale", None),
            ("H", "mag"),
            ("G", ""),
        ),
        "KEP": (
            ("id", None),
            ("a", "au"),
            ("e", ""),
            ("incl", "rad"),
            ("Omega", "rad"),
            ("w", "rad"),
            ("M", "rad"),
            ("orbtype", None),
            ("epoch", "d"),
            ("epoch_scale", None),
            ("H", "mag"),
            ("G", ""),
        ),
        "CART": (
            ("id", None),
            ("x", "au"),
            ("y", "au"),
            ("z", "au"),
            ("vx", "au/d"),
            ("vy", "au/d"),
            ("vz", "au/d"),
            ("orbtype", None),
            ("epoch", "d"),
            ("epoch_scale", None),
            ("H", "mag"),
            ("G", ""),
        ),
    }

    oorb_ephem_full_fields = [
        "MJD",
        "RA",
        "DEC",
        "RA*cos(Dec)_rate",
        "DEC_rate",
        "alpha",
        "elong",
        "r",
        "Delta",
        "V",
        "pa",
        "TopEclLon",
        "TopEclLat",
        "OppTopEclLon",
        "OppTopEclLat",
        "HelEclLon",
        "HelEclLat",
        "OppHelEclLon",
        "OppHelEclLat",
        "EL",
        "ELsun",
        "ELmoon",
        "lunarphase",
        "lunarelong",
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "obsx",
        "obsy",
        "obsz",
        "trueanom",
    ]

    oorb_ephem_full_units = [
        "d",
        "deg",
        "deg",
        "deg/d",
        "deg/d",
        "deg",
        "deg",
        "au",
        "au",
        "mag",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        "deg",
        None,
        "deg",
        "au",
        "au",
        "au",
        "au/d",
        "au/d",
        "au/d",
        "au",
        "au",
        "au",
        "deg",
    ]

    oorb_ephem_basic_fields = [
        "MJD",
        "RA",
        "DEC",
        "RA*cos(Dec)_rate",
        "DEC_rate",
        "alpha",
        "elong",
        "r",
        "Delta",
        "V",
        "trueanom",
    ]

    oorb_ephem_basic_units = [
        "d",
        "deg",
        "deg",
        "deg/d",
        "deg/d",
        "deg",
        "deg",
        "au",
        "au",
        "mag",
        "deg",
    ]

    # definitions for MPC orbits: MPC field name: [sbpy field name, unit]
    mpc_orbit_fields = {
        "absolute_magnitude": ["absmag", "mag"],
        "aphelion_distance": ["Q", "au"],
        "arc_length": ["arc", "day"],
        "argument_of_perihelion": ["w", "deg"],
        "ascending_node": ["Omega", "deg"],
        "delta_v": ["delta_v", "km/s"],
        "designation": ["desig", None],
        "earth_moid": ["moid_earth", "au"],
        "eccentricity": ["e", ""],
        "epoch_jd": ["epoch", "time_jd_utc"],
        "inclination": ["i", "deg"],
        "jupiter_moid": ["moid_jupiter", "au"],
        "mars_moid": ["moid_mars", "au"],
        "mean_anomaly": ["M", "deg"],
        "mercury_moid": ["moid_mercury", "au"],
        "name": ["name", None],
        "number": ["number", None],
        "orbit_type": ["mpc_orbit_type", None],
        "perihelion_date_jd": ["Tp", "time_jd_utc"],
        "perihelion_distance": ["q", "au"],
        "period": ["P", "year"],
        "phase_slope": ["G", ""],
        "saturn_moid": ["moid_saturn", "au"],
        "semimajor_axis": ["a", "au"],
        "tisserand_jupiter": ["Tj", ""],
        "uranus_moid": ["moid_uranus", "au"],
        "venus_moid": ["moid_venus", "au"],
    }


# clean namespace
del u, Time

from .core import (
    DataClass,
    DataClassError,
    QueryError,
    TimeScaleWarning,
)  # noqa: E402, E501
from .decorators import quantity_to_dataclass, dataclass_input  # noqa: E402
from .ephem import Ephem  # noqa: E402
from .orbit import Orbit  # noqa: E402
from .phys import Phys  # noqa: E402
from .obs import Obs  # noqa: E402
from .names import Names, natural_sort_key  # noqa: E402
