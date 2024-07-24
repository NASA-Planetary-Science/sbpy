Dust comae and tails (`sbpy.activity.dust`)
*******************************************

Cometary dust is the refractory material released by comets.  This sub-module provides simple photometric models of cometary dust comae.


*Afρ* and *εfρ* models
======================

`sbpy` has two classes to support observations and models of a coma continuum: `~sbpy.activity.dust.Afrho` and `~sbpy.activity.dust.Efrho`.

The *Afρ* parameter of A'Hearn et al (1984) is based on observations of idealized cometary dust comae.  It is proportional to the observed flux density within a circular aperture.  The quantity is the product of dust albedo, dust filling factor, and the radius of the aperture at the distance of the comet.  It carries the units of *ρ* (length), and under certain assumptions is proportional to the dust production rate of the comet:

.. math::

   Afρ = \frac{4 Δ^2 r_h^2}{ρ}\frac{F_c}{F_⊙}

where *Δ* and *ρ* have the same (linear) units, but :math:`r_h` is in units of au.  :math:`F_c` * is the flux density of the comet in the aperture, and :math:`F_⊙` is that of the Sun at 1 au in the same units.  See A'Hearn et al. (1984) and Fink & Rubin (2012) for more details.

The *εfρ* parameter is the thermal emission counterpart to *Afρ*, replacing albedo with IR emissivity, *ε*, and the solar spectrum with the Planck function, *B*:

.. math::

   εfρ = \frac{F_c Δ^2}{π ρ B(T_c)}

where :math:`T_c` is the spectral temperature of the continuum (Kelley et al. 2013).

*Afρ* and *εfρ* are quantities
------------------------------

``Afrho`` and ``Efrho`` are subclasses of `astropy`'s `~astropy.units.Quantity` and carry units of length.

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

They may be converted to other units of length just like any `~astropy.units.Quantity`:

   >>> afrho.to('m')    # doctest: +FLOAT_CMP
   <Afrho 1. m>

.. _afrho-to-from-flux-density:

Convert to/from flux density
----------------------------

The quantities may be initialized from flux densities.  Here, we reproduce one of the calculations from the original A'Hearn et al. (1984) work:

* :math:`r_h` = 4.785 au
* *Δ* = 3.822 au
* phase angle = 3.3°
* aperture radius = 9.8" (27200 km)
* :math:`\log{F_λ}` = -13.99 erg/(cm\ :sup:`2` s Å) at *λ* = 5240 Å

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

The ``Afrho`` class may be converted to a flux density, and the original value is recovered.

   >>> f = afrho.to_fluxd('λ5240', aper, eph).to('erg/(s cm2 AA)')
   >>> print(np.log10(f.value))    # doctest: +FLOAT_CMP
   -13.99

``Afrho`` works seamlessly with `sbpy`'s spectral calibration framework (:ref:`sbpy-calib`) when the `astropy` affiliated package `synphot` is installed.  The solar flux density (via `~sbpy.calib.solar_fluxd`) is not required, but instead the spectral wavelengths or the system transmission of the instrument and filter:

.. doctest-requires:: synphot

   >>> wave = [0.4, 0.5, 0.6] * u.um
   >>> print(afrho.to_fluxd(wave, aper, eph))    # doctest: +FLOAT_CMP
   [7.76770018e-14 1.05542873e-13 9.57978939e-14] W / (um m2)

.. doctest-requires:: synphot

   >>> from synphot import SpectralElement, Box1D
   >>> # bandpass width is a guess
   >>> bp = SpectralElement(Box1D, x_0=5240 * u.AA, width=50 * u.AA)
   >>> print(Afrho.from_fluxd(bp, flam, aper, eph))    # doctest: +FLOAT_CMP
   5994.110239075767 cm

Thermal emission with *εfρ*
---------------------------

The ``Efrho`` class has the same functionality as the ``Afrho`` class.  The most important difference is that *εfρ* is calculated using a Planck function and temperature.  `sbpy` follows common practice and parameterizes the temperature as a constant scale factor of :math:`T_{BB} = 278\,r_h^{1/2}`\  K, the equilibrium temperature of a large blackbody sphere at a distance :math:`r_h` from the Sun.

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


To/from magnitudes
------------------

``Afrho`` and ``Efrho`` also work with `astropy`'s magnitude units.  If the conversion between Vega-based magnitudes is required, then `sbpy`'s calibration framework (:ref:`sbpy-calib`) will be used.

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
--------------------------

Phase angle was not used in the previous section.  In the *Afρ* formalism, "albedo" includes the scattering phase function, and is more precisely written *A(θ)*, where *θ* is the phase angle.  The default behavior for ``Afrho`` is to compute *A(θ)fρ* as opposed to *A(0°)fρ*.  Returning to the A'Hearn et al. data, we scale *Afρ* to 0° from 3.3° phase using the :func:`~sbpy.activity.dust.Afrho.to_phase` method:

.. doctest-requires:: scipy

   >>> afrho = Afrho(6029.9 * u.cm)
   >>> print(afrho.to_phase(0 * u.deg, 3.3 * u.deg))    # doctest: +FLOAT_CMP
   6886.825981017757 cm

The default phase function is the Halley-Marcus composite phase function (:func:`~sbpy.activity.phase_HalleyMarcus`).  Any function or callable object that accepts an angle as a `~astropy.units.Quantity` and returns a scalar value may be used:

.. doctest-requires:: scipy

   >>> Phi = lambda phase: 10**(-0.016 / u.deg * phase.to('deg'))
   >>> print(afrho.to_phase(0 * u.deg, 3.3 * u.deg, Phi=Phi))    # doctest: +FLOAT_CMP
   6809.419810008357 cm

To correct an observed flux density for the phase function, use the ``phasecor`` option of :func:`~sbpy.activity.dust.Afrho.to_fluxd` and :func:`~sbpy.activity.dust.Afrho.from_fluxd` methods:

.. doctest-requires:: scipy

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


Using apertures
---------------

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




Reference/API
=============

.. automodapi:: sbpy.activity.dust
   :no-main-docstr:
   :inherited-members:
