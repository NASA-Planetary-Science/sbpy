Activity Module (`sbpy.activity`)
=================================

Introduction
------------

`sbpy.activity` models cometary dust and gas activity.

Dust Comae: *Afρ* and *εfρ*
---------------------------

`sbpy` has two classes to assist with observations and modeling of coma continuum: `~sbpy.activity.dust.Afrho` and `~sbpy.activity.dust.Efrho`.

The *Afρ* parameter of A'Hearn et al (1984) is based on observations of idealized cometary dust comae.  It is proportional to the observed flux density within a circular aperture.  The quantity is the product of dust albedo, dust filling factor, and the radius of the aperture at the distance of the comet.  It carries the units of *ρ* (length), and under certain assumptions is proportional to the dust production rate of the comet.  See A'Hearn et al. (1984) and Fink & Rubin (2012) for more details.  The *εfρ* parameter is the thermal emission counterpart to *Afρ*, replacing albedo with IR emissivity (Kelley et al. 2013).

*Afρ* and *εfρ* are quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Afrho` and `Efrho` are subclasses of `astropy`'s `~astropy.units.Quantity` and carry units of length.

  >>> import numpy as np
  >>> import astropy.units as u
  >>> from sbpy.activity import Afrho, Efrho
  >>>
  >>> afrho = Afrho(100 * u.cm)
  >>> print(afrho)    # doctest: +FLOAT_CMP
  100.0 cm
  >>> efrho = Efrho(afrho * 3.5)
  >>> print(efrho)    # doctest: +FLOAT_CMP
  350.0 cm

They may be converted to other units of length just like any `Quantity`:

  >>> afrho.to('m')    # doctest: +FLOAT_CMP
  <Afrho 1. m>


Flux density
^^^^^^^^^^^^

The quantities may be initialized from flux densities.  Here, we reproduce one of the calculations from the original A'Hearn et al. (1984) work:

* :math:`r_h` = 4.785 au
* *Δ* = 3.822 au
* phase angle = 3.3°
* aperture radius = 9.8" (27200 km)
* :math:`\log{F_λ}` = -13.99 erg/(cm2 s Å) at *λ* = 5240 Å

The solar flux density at 1 au is also needed.  We use 1868 W/(m2 μm).

  >>> from sbpy.data import Ephem
  >>> from sbpy.calib import solar_fluxd
  >>>
  >>> solar_fluxd.set({
  ...     'λ5240': 1868 * u.W / u.m**2 / u.um,
  ...     'λ5240(lambda pivot)': 5240 * u.AA
  ... })              # doctest: +IGNORE_OUTPUT
  >>>
  >>> flam = 10**-13.99 * u.Unit('erg/(s cm2 AA)')
  >>> aper = 27200 * u.km
  >>>
  >>> eph = Ephem.from_dict({'rh': 4.785 * u.au, 'delta': 3.822 * u.au})
  >>>
  >>> afrho = Afrho.from_fluxd('λ5240', flam, aper, eph)
  >>> print(afrho)    # doctest: +FLOAT_CMP
  6029.90248952895 cm

Which is within a few percent of 6160 cm computed by A'Hearn et al.. The difference is likely due to the assumed solar flux density in the bandpass.

The `Afrho` class may be converted to a flux density, and the original value is recovered.

  >>> f = afrho.to_fluxd('λ5240', aper, eph).to('erg/(s cm2 AA)')
  >>> print(np.log10(f.value))    # doctest: +FLOAT_CMP
  -13.99

`Afrho` works seamlessly with `sbpy`'s spectral calibration framework (:ref:`sbpy_calib`) when the `astropy` affiliated package `synphot` is installed.  The solar flux density (via `~sbpy.calib.solar_fluxd`) is not required, but instead the spectral wavelengths or the system transmission of the instrument and filter:

.. doctest-requires:: synphot

   >>> wave = [0.4, 0.5, 0.6] * u.um
   >>> print(afrho.to_fluxd(wave, aper, eph))    # doctest: +FLOAT_CMP
   [7.76770018e-14 1.05542873e-13 9.57978939e-14] W / (m2 um)

.. doctest-requires:: synphot

   >>> from synphot import SpectralElement, Box1D
   >>> # bandpass width is a guess
   >>> bp = SpectralElement(Box1D, x_0=5240 * u.AA, width=50 * u.AA)
   >>> print(Afrho.from_fluxd(bp, flam, aper, eph))    # doctest: +FLOAT_CMP
   5994.110239075767 cm

*εfρ*
^^^^^
   
The `Efrho` class has the same functionality as the `Afrho` class.  The most important difference is that *εfρ* is calculated using a Planck function and temperature.  `sbpy` follows common practice and parameterizes the temperature as a constant scale factor of :math:`T_{BB} = 278\,r_h^{1/2}` K, the equilibrium temperature of a large blackbody sphere at a distance :math:`r_h` from the Sun.

Reproduce the *εfρ* of 246P/NEAT from Kelley et al. (2013).

.. doctest-requires:: synphot

  >>> wave = [15.8, 22.3] * u.um
  >>> fluxd = [25.75, 59.2] * u.mJy
  >>> aper = 11.1 * u.arcsec
  >>> eph = Ephem.from_dict({'rh': 4.28 * u.au, 'delta': 3.71 * u.au})
  >>> efrho = Efrho.from_fluxd(wave, fluxd, aper, eph)
  >>> for i in range(len(wave)):
  ...     print('{:5.1f} at {:.1f}'.format(efrho[i], wave[i]))    # doctest: +FLOAT_CMP
  406.2 cm at 15.8 um
  427.9 cm at 22.3 um

Compare to 397.0 cm and 424.6 cm listed in Kelley et al. (2013).


Magnitudes
^^^^^^^^^^

`Afrho` and `Efrho` work with `astropy`'s magnitude units.  If the conversion between Vega-based magnitudes is required, then `sbpy`'s calibration framework (:ref:`sbpy_spectral_standards`) will be used.

Estimate the *Afρ* of comet C/2012 S1 (ISON) based on Pan-STARRS 1 photometry in the *r* band (Meech et al. 2013)

.. doctest-requires:: synphot

  >>> w = 0.617 * u.um
  >>> m = 16.02 * u.ABmag
  >>> aper = 5 * u.arcsec
  >>> eph = {'rh': 5.234 * u.au, 'delta': 4.275 * u.au, 'phase': 2.6 * u.deg}
  >>> afrho = Afrho.from_fluxd(w, m, aper, eph)
  >>> print(afrho)    # doctest: +FLOAT_CMP
  1948.496075629444 cm
  >>> m2 = afrho.to_fluxd(w, aper, eph, unit=u.ABmag)    # doctest: +FLOAT_CMP
  >>> print(m2)
  16.02 mag(AB)


Phase angles and functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Phase angle was not used in the previous section.  The default behavior for `Afrho` is to compute *A(θ)fρ* (as opposed to *A(0°)fρ*).  Returning to the A'Hearn et al. data, we scale *Afρ* to 0° from 3.3° phase using the :func:`~sbpy.activity.Afrho.to_phase` method:

  >>> afrho = Afrho(6029.9 * u.cm)
  >>> print(afrho.to_phase(0 * u.deg, 3.3 * u.deg))    # doctest: +FLOAT_CMP
  6886.825981017757 cm

The default phase function is the Halley-Marcus composite phase function (:func:`~sbpy.activity.phase_HalleyMarcus`), but any callable that returns a scale factor from an angle may be used:

  >>> Phi = lambda phase: 10**(-0.016 / u.deg * phase.to('deg'))
  >>> print(afrho.to_phase(0 * u.deg, 3.3 * u.deg, Phi=Phi))    # doctest: +FLOAT_CMP
  6809.419810008357 cm

To make a phase function correction on an observed flux density, use the ``phasecor`` option of :func:`~sbpy.activity.Afrho.to_fluxd` and :func:`~sbpy.activity.Afrho.from_fluxd` methods:

  >>> flam = 10**-13.99 * u.Unit('erg/(s cm2 AA)')
  >>> aper = 27200 * u.km
  >>> eph = Ephem.from_dict({
  ...     'rh': 4.785 * u.au,
  ...     'delta': 3.822 * u.au,
  ...     'phase': 3.3 * u.deg
  ... })
  >>>
  >>> afrho = Afrho.from_fluxd('λ5240', flam, aper, eph, phasecor=True)
  >>> print(afrho)    # doctest: +FLOAT_CMP
  6886.828824340642 cm


Apertures
^^^^^^^^^

Other apertures may be used, as long as they can be converted into an equivalent radius, assuming a coma with a *1/ρ* surface brightness distribution.  `~sbpy.activity` has a collection of useful geometries.

  >>> from sbpy.activity import CircularAperture, AnnularAperture, RectangularAperture, GaussianAperture
  >>> apertures = (
  ...   ( '10" radius circle', CircularAperture(10 * u.arcsec)),
  ...   (    '5"–10" annulus', AnnularAperture([5, 10] * u.arcsec)),
  ...   (       '2"x10" slit', RectangularAperture([2, 10] * u.arcsec)),
  ...   ('σ=5" Gaussian beam', GaussianAperture(5 * u.arcsec))
  ... )
  >>> for name, aper in apertures:
  ...     afrho = Afrho.from_fluxd('λ5240', flam, aper, eph)
  ...     print('{:18s} = {:5.0f}'.format(name, afrho))    # doctest: +FLOAT_CMP
   10" radius circle =  5917 cm
      5"–10" annulus = 11834 cm
         2"x10" slit = 28114 cm
  σ=5" Gaussian beam =  9442 cm


Production Rate calculations
----------------------------

`~sbpy.activity.gas.productionrate` offers various functions that aid in the calculation
of production rates. `~sbpy.data.phys` has a function called `~sbpy.data.Phys.from_jplspec` which takes care of querying the JPL Molecular Spectral Catalog through the use of `~astroquery.jplspec` and calculates all the necessary constants needed for
production rate calculations in this module. Yet, the option for the user to
provide their own molecular data is possible through the use of an `~sbpy.data.phys` object, as long as it has the required information. It is imperative to read
the documentation of the functions in this section to understand what is needed
for each. If the user does not have the necessary data, they can build an object
using JPLSpec:

.. doctest-skip::

    >>> from sbpy.data.phys import Phys
    >>> import astropy.units as u
    >>> temp_estimate = 47. * u.K
    >>> transition_freq = (230.53799 * u.GHz).to('MHz')
    >>> mol_tag = '^CO$'
    >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)

Having this information, we can move forward towards the calculation of
production rate. The functions that sbpy currently provides to calculate
production rates are listed below.

Integrated Line Intensity Conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The JPL Molecular Spectroscopy Catalog offers the integrated line intensity
at 300 K for a molecule. Yet, in order to calculate production rate, we need
to know the integrated line intensity at a given temperature. This function
takes care of converting the integrated line intensity at 300 K to its equivalent
in the desired temperature using equations provided by the JPLSpec documentation.
For more information on the needed parameters for this function follow the link
for `~sbpy.activity.intensity_conversion` under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import intensity_conversion
    >>> intl = intensity_conversion(mol_data)
    >>> mol_data.add_column([intl.value] * intl.unit, name='intl')
     11
    >>> intl
     <Quantity 0.00280051 MHz nm2>


Einstein Coefficient Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike catalogs like LAMDA, JPLSpec does not offer the Eistein coefficient
and it must be calculated using equations provided by the JPL Molecular
Spectroscopy Catalog. These equations have been compared to established LAMDA
values of the Einstein Coefficient for HCN and CO, and no more than
a 24% difference has been found between the calculation from JPLSpec and the
LAMDA catalog value. Since JPLSpec and LAMDA are two very different
catalogs with different data, the difference is expected, and the user is
allowed to provide their own Einstein Coefficient if they want. If the user
does want to provide their own Einstein Coefficient, they may do so simply
by appending their value with the unit 1/s to the `~sbpy.data.Phys` object, called
`mol_data` in these examples. For more information on the needed parameters
for this function follow the link for `~sbpy.activity.einstein_coeff`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import einstein_coeff
    >>> au = einstein_coeff(mol_data)
    >>> mol_data.add_column([au.value] * au.unit, name = 'Einstein Coefficient')
     12
    >>> au
      <Quantity 7.03946054e-07 1 / s>


Beta Factor Calculation
^^^^^^^^^^^^^^^^^^^^^^^

Returns beta factor based on timescales from `~sbpy.activity.gas` and distance
from the Sun using an `~sbpy.data.ephem` object. The calculation is
parent photodissociation timescale * (distance from comet to Sun)**2
and it accounts for certain photodissociation and geometric factors needed
in the calculation of total number of molecules `~sbpy.activity.total_number_nocd`
If you wish to provide your own beta factor, you can calculate the equation
expressed in units of AU**2 * s , all that is needed is the timescale
of the molecule and the distance of the comet from the Sun. Once you
have the beta factor you can append it to your mol_data phys object
with the name 'beta' or any of its alternative names. For more information on
the needed parameters for this function follow the link for
`~sbpy.activity.beta_factor` under Reference/API section.

.. doctest-skip::

    >>> from astropy.time import Time
    >>> from sbpy.data import Ephem
    >>> from sbpy.activity import beta_factor
    >>> target = 'C/2016 R2'
    >>> time = Time('2017-12-22 05:24:20', format = 'iso')
    >>> ephemobj = Ephem.from_horizons(target, epochs=time.jd)
    >>> beta = beta_factor(mol_data, ephemobj)
    >>> mol_data.add_column([beta.value] * beta.unit, name='beta')
     13
    >>> beta
     <Quantity [13333365.25745597] AU2 s>


Total Number
^^^^^^^^^^^^

In order to obtain our total number of molecules from flux data,
with no column density information, we can use equation 10 from `Bockelee-Morvan
et al. 2004 <https://ui.adsabs.harvard.edu/#abs/2004come.book..391B>`_ and the
millimeter/submillimeter spectroscopy beam factors explained and detailed
in equation 1.3 from:

    | Drahus, M. (2010). Microwave observations and modeling of the molecular
    | coma in comets. PhD Thesis, Georg-August-Universität Göttingen.

If the user prefers to give the total number, they may do so by appending
to the mol_data `~sbpy.data.phys` object with the name `total_number_nocd` or
any of its alternative names. For more information on the needed parameters
for this function follow the link for `~sbpy.activity.total_number_nocd`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import total_number_nocd
    >>> integrated_flux = 0.26 * u.K * u.km / u.s
    >>> b = 0.74
    >>> aper = 10 * u.m
    >>> tnum = total_number_nocd(integrated_flux, mol_data, aper, b)
    >>> mol_data.add_column([tnum], name='total_number_nocd')
     14
    >>> tnum
     <Quantity [2.93988826e+26]>


Simplified Model for Production Rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`~sbpy.activity` provides several models to calculate production rates in comets.
One of the models followed by this module is based on the following paper:

| Drahus et al. September 2012. The Sources of HCN and CH3OH and the
| Rotational Temperature in Comet 103P/Hartley 2 from Time-resolved
| Millimeter Spectroscopy. The Astrophysical Journal, Volume 756,
| Issue 1.

The following example shows the usage of the function. This LTE model does not
include photodissociation, but it does serve as way to obtain educated
first guesses for other models within sbpy. For more information on the
parameters that are needed for the function follow
the link for the function `from_Drahus` in `~sbpy.activity.LTE`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import LTE
    >>> vgas = 0.5 * u.km / u.s
    >>> lte = LTE()
    >>> q = lte.from_Drahus(integrated_flux, mol_data, ephemobj, vgas, aper, b=b)
    >>> q
     <Quantity 3.59397119e+28 1 / s>


Haser Model for Production Rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another model included in the module is based off of the model in the following
literature:

| Haser 1957, Bulletin de la Societe Royale des Sciences de Liege 43, 740.
| Newburn and Johnson 1978, Icarus 35, 360-368.

This model is well-known as the Haser model. In the case of our implementation
the function takes in an initial guess for the production rate, and uses the
module found in `~sbpy.activity.gas` to find a ratio between the model
total number of molecules and the number of molecules calculated from the data
to scale the model Q and output the new production rate from the result. This
LTE model does account for the effects of photolysis. For more information
on the parameters that are needed for the function follow
the link for the function `from_Haser` in `~sbpy.activity.LTE`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import Haser, photo_timescale
    >>> Q_estimate = 3.5939*10**(28) / u.s
    >>> parent = photo_timescale('CO') * vgas
    >>> coma = Haser(Q_estimate, vgas, parent)
    >>> lte = LTE()
    >>> Q = lte.from_Haser(coma, mol_data, aper=aper)
    >>> Q
     <Quantity [[9.35795579e+27]] 1 / s>


Reference/API
-------------
.. automodapi:: sbpy.activity
    :no-heading:
