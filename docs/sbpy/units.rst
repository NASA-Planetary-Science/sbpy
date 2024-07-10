Units Module (`sbpy.units`)
===========================

Introduction
------------

`~sbpy.units` provides common planetary astronomy units and unit conversions, including Vega-based magnitudes.  sbpy units may be added to the top-level `astropy.units` namespace via:

  >>> import sbpy.units
  >>> sbpy.units.enable()    # doctest: +IGNORE_OUTPUT

This will make them searchable via strings, such as:

  >> print(u.Unit("mag(VEGA)"))
  VEGAmag


Spectral Gradient Units
-----------------------

Spectral gradients are commonly expressed as % per 100 nm.  This is too subtle a concept for astropy's `~astropy.units`:

  >>> import astropy.units as u
  >>> print(u.percent / (100 * u.nm))
  0.01 % / nm
  >>> print(u.Unit('% / (100 nm)'))    # doctest: +SKIP
  ValueError: '% / (100 nm)' did not parse as unit: Syntax error parsing unit '% / 100 nm'

As a convenience, sbpy defines the `~sbpy.units.hundred_nm` unit, which has an appropriate string representation:

  >>> from sbpy.units import hundred_nm
  >>> print(u.percent / hundred_nm)
  % / (100 nm)

.. _vega-magnitudes:

Vega Magnitudes
---------------

With the synphot package, sbpy has the capability to convert between flux densities and Vega-based magnitude systems.  The conversions require a spectrum of Vega, which is provided by sbpy's :ref:`calibration system <sbpy-calib>`.

sbpy defines two new spectral flux density units: ``VEGA`` and ``JM``.  ``VEGA`` represents the flux density of Vega.  ``JM`` represents the flux density zeropoint of the Johnson-Morgan system, assuming Vega has a magnitude of 0.03 at all wavelengths (`Johnson et al. 1966 <https://ui.adsabs.harvard.edu/abs/1966CoLPL...4...99J>`_, `Bessell & Murphy 2012 <https://ui.adsabs.harvard.edu/abs/2012PASP..124..140B>`_).  Two magnitude units are also defined: ``VEGAmag`` and ``JMmag``.

Unit conversions between flux density and Vega-based magnitudes use the `astropy.units equivalency system <https://docs.astropy.org/en/stable/units/equivalencies.html#unit-equivalencies>`_.  sbpy's :func:`~sbpy.units.spectral_density_vega` provides the equivalencies, which astropy would use to convert the units.  The function requires a wavelength, frequency, or bandpass:

.. doctest-requires:: synphot

  >>> import astropy.units as u
  >>> from sbpy.units import VEGAmag, spectral_density_vega
  >>> wave = 5500 * u.AA
  >>> m = 0 * VEGAmag
  >>> fluxd = m.to('erg/(cm2 s AA)', spectral_density_vega(wave))
  >>> fluxd.value                                  # doctest: +FLOAT_CMP
  3.5469235179497687e-09
  >>> m = fluxd.to(VEGAmag, spectral_density_vega(wave))
  >>> m.value                                      # doctest: +FLOAT_CMP
  0.0

To use a bandpass, define and pass a `synphot.spectrum.SpectralElement`.  A limited set of bandpasses are distributed with sbpy (see :ref:`filter-bandpasses`):

.. doctest-requires:: synphot

  >>> from sbpy.units import VEGAmag, spectral_density_vega
  >>> from sbpy.photometry import bandpass
  >>> V = bandpass('Johnson V')
  >>> m = 0.0 * VEGAmag
  >>> fluxd = m.to('erg/(cm2 s AA)', spectral_density_vega(V))
  >>> fluxd.value                   # doctest: +FLOAT_CMP
  3.5469235114856157e-09

.. _reflectance-equivalencies:

Reflectance Units
-----------------

``sbpy.units`` supports the conversions between reflected flux (often
expressed in magnitude), and bidirectional reflectance (in unit of 1/sr) and
scattering cross-section through function `~sbpy.units.reflectance`.

For example, the absolute magnitude of Ceres in V-band is 3.4 in ``VEGAmag``
system, the radius is 460 km, its disk-averaged bidirectional reflectance at 0
phase angle can be calculated:

  >>> import numpy as np
  >>> from astropy import units as u
  >>> from sbpy.calib import solar_fluxd, vega_fluxd
  >>> from sbpy.units import reflectance, VEGAmag, spectral_density_vega
  >>>
  >>> solar_fluxd.set({"V": -26.775 * VEGAmag, "V(lambda pivot)": 0.55 * u.um})
  <ScienceState solar_fluxd: {'V': <Magnitude -26.775 mag(VEGA)>, 'V(lambda pivot)': <Quantity 0.55 um>}>
  >>> vega_fluxd.set({"V": 3.5885e-08 * u.Unit("W / (m2 um)"), "V(lambda pivot)": 0.55 * u.um})
  <ScienceState vega_fluxd: {'V': <Quantity 3.5885e-08 W / (um m2)>, 'V(lambda pivot)': <Quantity 0.55 um>}>
  >>> mag = 3.4 * VEGAmag
  >>> radius = 460 * u.km
  >>> cross_sec = np.pi * (radius)**2
  >>> ref = mag.to('1/sr', reflectance('V', cross_section=cross_sec))
  >>> print('{0:.4f}'.format(ref))
  0.0287 1 / sr

`~sbpy.units.reflectance` works with `sbpy`'s spectral calibration system (see :ref:`sbpy-calib`):

.. doctest-requires:: synphot

  >>> from sbpy.photometry import bandpass
  >>> V = bandpass('Johnson V')
  >>> ref = 0.0287 / u.sr
  >>> cross_sec = mag.to('km2', reflectance(V, reflectance=ref))
  >>> radius = np.sqrt(cross_sec / np.pi)
  >>> print('{0:.2f}'.format(radius))
  459.69 km

`~sbpy.units.reflectance` also supports conversion between a flux spectrum and a reflectance spectrum:

.. doctest-requires:: synphot

  >>> wave = [1.046, 1.179, 1.384, 1.739, 2.416] * u.um
  >>> flux = [1.636e-18, 1.157e-18, 8.523e-19, 5.262e-19, 1.9645e-19] \
  ...     * u.Unit('W/(m2 um)')
  >>> xsec = 0.0251 * u.km**2
  >>> ref = flux.to('1/sr', reflectance(wave, cross_section=xsec))
  >>> print(ref)  # doctest: +FLOAT_CMP
  [0.0021763  0.00201223 0.0022041  0.00269637 0.00292785] 1 / sr

`~sbpy.units.reflectance` also supports dimensionless logarithmic unit for
solar flux, which can be specified with `~sbpy.calib.solar_fluxd.set`:

  >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
  ...     mag = 3.4 * u.mag
  ...     ref = mag.to('1/sr', reflectance('V', cross_section =
  ...         np.pi*(460*u.km)**2))
  >>> print(ref)  # doctest: +FLOAT_CMP
  0.028786262941247264 1 / sr

Projected Sizes
---------------

With the `~sbpy.units.projected_size` equivalencies, one may convert between angles and lengths, for a given distance:

  >>> import astropy.units as u
  >>> import sbpy.units as sbu
  >>> (1 * u.arcsec).to('km', sbu.projected_size(1 * u.au))
  ... # doctest: +FLOAT_CMP
  <Quantity [725.27094381] km>
  >>> (725.27 * u.km).to('arcsec', sbu.projected_size(1 * u.au))
  ... # doctest: +FLOAT_CMP
  <Quantity [1.00] arcsec>


Reference/API
-------------
.. automodapi:: sbpy.units
    :no-heading:
    :include-all-objects:
