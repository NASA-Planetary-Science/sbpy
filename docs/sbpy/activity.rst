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
  >>> print(afrho)
  100.0 cm
  >>> efrho = Efrho(afrho * 3.5)
  >>> print(efrho)
  350.0 cm

They may be converted to other units of length just like any `Quantity`:

  >>> afrho.to('m')
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

  >>> wave = 5240 * u.AA
  >>> flam = 10**-13.99 * u.Unit('erg/(s cm2 AA)')
  >>> aper = 27200 * u.km
  >>> eph = dict(rh=4.785 * u.au, delta=3.822 * u.au)
  >>> Slam = 1868 * u.W / u.m**2 / u.um
  >>> 
  >>> afrho = Afrho.from_fluxd(wave, flam, aper, eph, S=Slam)
  >>> print(afrho)    # doctest: +FLOAT_CMP
  6029.90248952895 cm

Which is within a few percent of 6160 cm computed by A'Hearn et al.. The difference is likely due to the assumed solar flux density in the bandpass.

The `Afrho` class may be converted to a flux density, and the original value is recovered.

  >>> f = afrho.to_fluxd(wave, aper, eph, S=Slam).to('erg/(s cm2 AA)')
  >>> print('''     F_λ = {}
  ... log(F_λ) = {}'''.format(f, np.log10(f.value)))    # doctest: +FLOAT_CMP
  F_λ = 1.0232929922807537e-14 erg / (Angstrom cm2 s)
  log(F_λ) = -13.99

The :func:`~sbpy.activity.Afrho.to_fluxd` and :func:`~sbpy.activity.Afrho.from_fluxd` methods work with units of flux density per wavelength or frequency.

  >>> fnu = flam.to('Jy', u.spectral_density(wave))
  >>> print(fnu)    # doctest: +FLOAT_CMP
  0.009372206976883997 Jy
  >>> Snu = 1.711e14 * u.Jy
  >>> print(Afrho.from_fluxd(wave, fnu, aper, eph, S=Snu))    # doctest: +FLOAT_CMP
  6029.468388208903 cm

`Afrho` works seamlessly with `sbpy`'s calibration framework (:ref:`sbpy_spectral_standards`) when the `astropy` affiliated package `synphot` is installed.  To convert to flux density using the default solar spectrum omit the `S` parameter:

.. doctest-requires:: synphot

   >>> wave = [0.4, 0.5, 0.6] * u.um
   >>> print(afrho.to_fluxd(wave, aper, eph))    # doctest: +FLOAT_CMP
   [7.76770018e-14 1.05542873e-13 9.57978939e-14] W / (m2 um)

To use the Kurucz (1993) model:

.. doctest-requires:: synphot

   >>> from sbpy.spectroscopy.sun import default_sun
   >>> with default_sun.set('Kurucz1993'):            # doctest: +REMOTE_DATA
   ...     print(afrho.to_fluxd(wave, aper, eph))    # doctest: +FLOAT_CMP
   [7.62582935e-14 1.06322888e-13 9.55650074e-14] W / (m2 um)
   
The `Efrho` class has the same functionality as the `Afrho` class.  The most important difference is that *εfρ* is calculated using a Planck function and temperature.  `sbpy` follows common practice and parameterizes the temperature as a constant scale factor of :math:`T_{BB} = 278\,r_h^{1/2}` K, the equilibrium temperature of a large blackbody sphere at a distance :math:`r_h` from the Sun.

Reproduce the *εfρ* of 246P/NEAT from Kelley et al. (2013).

  >>> wave = [15.8, 22.3] * u.um
  >>> fluxd = [25.75, 59.2] * u.mJy
  >>> aper = 11.1 * u.arcsec
  >>> eph = {'rh': 4.28 * u.au, 'delta': 3.71 * u.au}
  >>> efrho = Efrho.from_fluxd(wave, fluxd, aper, eph)
  >>> for i in range(len(wave)):
  ...     print('{:5.1f} at {:.1f}'.format(efrho[i], wave[i]))  # doctest: +FLOAT_CMP
  406.2 cm at 15.8 um
  427.9 cm at 22.3 um

Compare to 397.0 cm and 424.6 cm listed in Kelley et al. (2013).


Magnitudes
^^^^^^^^^^

`Afrho` and `Efrho` work seamlessly with `astropy`'s magnitude units.  If the conversion between Vega-based magnitudes is required, then `sbpy`'s calibration framework (:ref:`sbpy_spectral_standards`) will be used.

Estimate the *Afρ* of comet C/2012 S1 (ISON) based on Pan-STARRS 1 photometry in the *r* band (Meech et al. 2013)

.. doctest-requires:: synphot

  >>> w = 0.617 * u.um
  >>> m = 16.02 * u.ABmag
  >>> aper = 5 * u.arcsec
  >>> eph = {'rh': 5.234 * u.au, 'delta': 4.275 * u.au, 'phase': 2.6 * u.deg}
  >>> afrho = Afrho.from_fluxd(w, m, aper, eph)
  >>> print(afrho)    # doctest: +FLOAT_CMP
  1948.496075629444 cm
  >>> m2 = afrho.to_fluxd(w, aper, eph, unit=u.ABmag)   # doctest: +FLOAT_CMP
  >>> print(m2)
  16.02 mag(AB)


Phase angles and functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Phase angle was not used in the previous section.  The default behavior for `Afrho` is to compute *A(θ)fρ* (as opposed to *A(0°)fρ*).  Returning to the A'Hearn et al. data, we scale *Afρ* to 0° from to 3.3° phase using the :func:`~sbpy.activity.Afrho.to_phase` method:

  >>> afrho = Afrho(6029.9 * u.cm)
  >>> print(afrho.to_phase(0 * u.deg, 3.3 * u.deg))  # doctest: +FLOAT_CMP
  6886.825981017757 cm

The default phase function is the Halley-Marcus composite phase function (:func:`~sbpy.activity.phase_HalleyMarcus`), but any callable that returns a scale factor from an angle:

  >>> Phi = lambda phase: 10**(-0.016 / u.deg * phase.to('deg'))
  >>> print(afrho.to_phase(0 * u.deg, 3.3 * u.deg, Phi=Phi))    # doctest: +FLOAT_CMP
  6809.419810008357 cm

*Afρ* can also be scaled with the ``phasecor`` option in the :func:`~sbpy.activity.Afrho.to_fluxd` and :func:`~sbpy.activity.Afrho.from_fluxd` methods:

  >>> wave = 5240 * u.AA
  >>> flam = 10**-13.99 * u.Unit('erg/(s cm2 AA)')
  >>> aper = 27200 * u.km
  >>> eph = dict(rh=4.785 * u.au, delta=3.822 * u.au, phase=3.3 * u.deg)
  >>> Slam = 1868 * u.W / u.m**2 / u.um
  >>> afrho = Afrho.from_fluxd(wave, flam, aper, eph, S=Slam, phasecor=True)
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
  ...     afrho = Afrho.from_fluxd(wave, flam, aper, eph, S=Slam)
  ...     print('{:18s} = {:5.0f}'.format(name, afrho))    # doctest: +FLOAT_CMP
   10" radius circle =  5917 cm
      5"–10" annulus = 11834 cm
         2"x10" slit = 28114 cm
  σ=5" Gaussian beam =  9442 cm



Reference/API
-------------
.. automodapi:: sbpy.activity
    :no-heading:
