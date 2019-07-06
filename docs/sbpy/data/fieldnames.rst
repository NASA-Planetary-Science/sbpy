
.. _field name list:

sbpy Field Names
================

The following table lists field names that are recognized by `sbpy`
when accessing `~sbpy.data.DataClass` objects, i.e.,
`~sbpy.data.Ephem`, `~sbpy.data.Orbit`, or `~sbpy.data.Phys`
objects. Each row of the following table represents one property; for
each property it lists its description, acceptable field names,
provenance (which `~sbpy.data.DataClass` class should be used to store
this object so that `sbpy` uses it properly), and its physical
dimension (if any).

How do I use this Table?
------------------------

As an example, imagine you are interested in storing an object's right
ascension into a `~sbpy.data.DataClass` object. The field names table
tells you that you should name the field either ``ra`` or ``RA``, that
you should use either a `~sbpy.data.Ephem` or `~sbpy.data.Obs` object
to store the data in, and that the field data should be expressed as
angles. Based on this information, we can create a `~sbpy.data.Obs`
object (presuming that the data were derived from observations):

    >>> from sbpy.data import Obs
    >>> import astropy.units as u
    >>> obs = Obs.from_dict({'ra': [12.345, 12.346, 12.347]*u.deg})
    >>> obs['ra']  # doctest: +SKIP
    <Quantity [12.345, 12.346, 12.347] deg>

Since RA requires an angle as dimension, we use degrees, but we might
as well use radians - `sbpy` will convert the units where necessary.
RA has an alternative field name (``'RA'``), we can now use that name,
too, in order to retrieve the data:

    >>> obs['RA']  # doctest: +SKIP
    <Quantity [12.345, 12.346, 12.347] deg>


The field name list is always up to date, but it might not be
complete. If you think an important alternative name is missing,
please suggest it by opening an issue. However, keep in mind that each
alternative field name has to be **unique** and **unambiguous**. The
source list is located as ``sbpy.data.conf.fieldnames`` in
``sbpy/data/__init__.py``.


Field Name List
---------------

======================================================= =================================================================================================================== =========================================================================== ================
                                            Description                                                                                                         Field Names                                                                  Provenance        Dimension
======================================================= =================================================================================================================== =========================================================================== ================
                                  **Target Identifier**                                                                                  ``targetname``, ``id``, ``Object`` `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`, `~sbpy.data.Phys`             None
                                              **Epoch**                                          ``epoch``, ``datetime_jd``, ``JD``, ``Date``, ``date``, ``Time``, ``time``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`             time
                                    **Semi-Major Axis**                                                                                                      ``a``, ``sma``                                                          `~sbpy.data.Orbit`           length
                                       **Eccentricity**                                                                                                      ``e``, ``ecc``                                                          `~sbpy.data.Orbit`             None
                                        **Inclination**                                                                                            ``i``, ``inc``, ``incl``                                                          `~sbpy.data.Orbit`            angle
                    **Longitude of the Ascending Node**                                                                                   ``Omega``, ``longnode``, ``node``                                                          `~sbpy.data.Orbit`            angle
                          **Argument of the Periapsis**                                                                                                   ``w``, ``argper``                                                          `~sbpy.data.Orbit`            angle
                                       **Mean Anomaly**                                                                                                ``M``, ``mean anom``                                                          `~sbpy.data.Orbit`            angle
                                       **True Anomaly**                                                                                                ``v``, ``true_anom``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                              **Heliocentric Distance**                                                                               ``r``, ``rh``, ``r_hel``, ``heldist``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`           length
                           **Distance to the Observer**                                                                                   ``delta``, ``Delta``, ``obsdist``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`           length
                                    **Right Ascension**                                                                                                      ``ra``, ``RA``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                                        **Declination**                                                                                                    ``dec``, ``DEC``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                               **Right Ascension Rate**                                              ``ra_rate``, ``RA_rate``, ``ra_rates``, ``RA_rates``, ``dRA``, ``dra``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs` angular velocity
                                   **RA*cos(Dec) Rate**                        ``RA*cos(Dec)_rate``, ``dra cos(dec)``, ``dRA cos(Dec)``, ``dra*cos(dec)``, ``dRA*cos(Dec)``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs` angular velocity
                                   **Declination Rate** ``dec_rate``, ``DEC_rate``, ``Dec_rate``, ``dec_rates``, ``DEC_rates``, ``Dec_rates``, ``dDec``, ``dDEC``, ``ddec``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs` angular velocity
                                      **Proper Motion**                                                                                           ``mu``, ``Proper motion``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs` angular velocity
                            **Proper Motion Direction**                                                                                        ``Direction``, ``direction``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                                  **Solar Phase Angle**                                                                                ``alpha``, ``phaseangle``, ``Phase``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                             **Solar Elongation Angle**                                      ``elong``, ``solarelong``, ``solarelongation``, ``elongation``, ``Elongation``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                                   **V-band Magnitude**                                                                                                     ``V``, ``Vmag``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`        magnitude
                    **Heliocentric Ecliptic Longitude**                                                      ``hlon``, ``EclLon``, ``ecllon``, ``HelEclLon``, ``helecllon``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                     **Heliocentric Ecliptic Latitude**                                                      ``hlat``, ``EclLat``, ``ecllat``, ``HelEclLat``, ``helecllat``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                               **Horizontal Elevation**                                                  ``el``, ``EL``, ``elevation``, ``alt``, ``altitude``, ``Altitude``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                                 **Horizontal Azimuth**                                                                                         ``az``, ``AZ``, ``azimuth``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                                   **Lunar Elongation**                          ``lunar_elong``, ``elong_moon``, ``elongation_moon``, ``lunar_elongation``, ``lunarelong``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`            angle
                           **X State Vector Component**                                                                                             ``x``, ``X``, ``x_vec``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`           length
                           **Y State Vector Component**                                                                                             ``y``, ``Y``, ``y_vec``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`           length
                           **Z State Vector Component**                                                                                             ``z``, ``Z``, ``z_vec``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`           length
                        **X Velocity Vector Component**                                                                                           ``vx``, ``dx``, ``dx/dt``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`         velocity
                        **Y Velocity Vector Component**                                                                                           ``vy``, ``dy``, ``dy/dt``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`         velocity
                        **Z Velocity Vector Component**                                                                                           ``vz``, ``dz``, ``dz/dt``                    `~sbpy.data.Orbit`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`         velocity
                         **Infrared Beaming Parameter**                                                                                                    ``eta``, ``Eta``                                        `~sbpy.data.Ephem`, `~sbpy.data.Obs`             None
                                        **Temperature**                                                                ``temp``, ``Temp``, ``temperature``, ``Temperature``                     `~sbpy.data.Phys`, `~sbpy.data.Ephem`, `~sbpy.data.Obs`      temperature
                                 **Effective Diameter**                                                                  ``d``, ``D``, ``diam``, ``diameter``, ``Diameter``                                                           `~sbpy.data.Phys`           length
                                   **Effective Radius**                                                                                                   ``R``, ``radius``                                                           `~sbpy.data.Phys`           length
                                   **Geometric Albedo**                                                                       ``pv``, ``pV``, ``p_v``, ``p_V``, ``geomalb``                                                           `~sbpy.data.Phys`             None
                                        **Bond Albedo**                                                                                               ``A``, ``bondalbedo``                                                           `~sbpy.data.Phys`             None
                                         **Emissivity**                                                                                      ``emissivity``, ``Emissivity``                                                           `~sbpy.data.Phys`             None
                                **Molecule Identifier**                                                                                           ``mol_tag``, ``mol_name``                                                           `~sbpy.data.Phys`             None
                               **Transition frequency**                                                                                                          ``t_freq``                                                           `~sbpy.data.Phys`        frequency
                 **Integrated line intensity at 300 K**                                                                                                        ``lgint300``                                                           `~sbpy.data.Phys`        intensity
**Integrated line intensity at designated Temperature**                                                                                                 ``intl``, ``lgint``                                                           `~sbpy.data.Phys`        intensity
                        **Partition function at 300 K**                                                                                                       ``partfn300``                                                           `~sbpy.data.Phys`             None
       **Partition function at designated temperature**                                                                                                          ``partfn``                                                           `~sbpy.data.Phys`             None
                             **Upper state degeneracy**                                                                                                            ``dgup``                                                           `~sbpy.data.Phys`             None
                       **Upper level energy in Joules**                                                                                                ``eup_j``, ``eup_J``                                                           `~sbpy.data.Phys`           energy
                       **Lower level energy in Joules**                                                                                                ``elo_j``, ``elo_J``                                                           `~sbpy.data.Phys`           energy
                                 **Degrees of freedom**                                                                                  ``degfr``, ``ndf``, ``degfreedom``                                                           `~sbpy.data.Phys`             None
                               **Einstein Coefficient**                                                                                                ``au``, ``eincoeff``                                                           `~sbpy.data.Phys`             None
                                    **Timescale * r^2**                                                                                           ``beta``, ``beta_factor``                                                           `~sbpy.data.Phys`             time
                                       **Total Number**                                                                                   ``totnum``, ``total_number_nocd``                                                           `~sbpy.data.Phys`             None
======================================================= =================================================================================================================== =========================================================================== ================
