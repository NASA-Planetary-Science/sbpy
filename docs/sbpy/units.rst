Photometry Module (`sbpy.units`)
=====================================

Introduction
------------

`~sbpy.units` provides common planetary astronomy units and Vega-based
magnitude conversions.  ``sbpy`` units may be added to the top-level
``astropy.units`` namespace via:

  >>> import sbpy.units
  >>> sbpy.units.enable()    # doctest: +SKIP

This will make them searchable via strings, such as
``u.Unit('mag(VEGA)')``.


Spectral Gradient Units
-----------------------

Spectral gradients are commonly expressed as % per 100 nm.  This is a
subtle concept for ``astropy``'s `~astropy.units`:

  >>> import astropy.units as u
  >>> print(u.percent / (100 * u.nm))
  0.01 % / nm
  >>> print(u.Unit('% / (100 nm)'))    # doctest: +SKIP
  ValueError: '% / (100 nm)' did not parse as unit: Syntax error parsing unit '% / 100 nm'

As a convenience, ``sbpy`` defines a ``hundred_nm`` unit that has an
appropriate string representation:

  >>> from sbpy.units import hundred_nm
  >>> print(u.percent / hundred_nm)
  % / 100 nm


Vega Magnitudes
---------------

With the `~synphot` package, ``sbpy`` has the capability to convert
between flux densities and Vega-based magnitude systems.
The conversions require a spectrum of Vega.
``sbpy`` includes the dust-free spectrum from `Bohlin 2014 <https://ui.adsabs.harvard.edu/#abs/2014AJ....147..127B/abstract>`_, but any user-provided spectrum may be used (`:ref:_sbpy_spectral_standards`).

``sbpy`` defines two new units: ``VEGA`` and ``JM``.  ``VEGA`` represents the flux density of Vega.  ``JM`` represents the flux density zeropoint of the Johnson-Morgan system, assuming Vega has a magnitude of 0.03 at all wavelengths (`Johnson et al. 1966 <https://ui.adsabs.harvard.edu/#abs/1966CoLPL...4...99J/abstract>`_, `Bessell & Murphy 2012 <https://ui.adsabs.harvard.edu/#abs/2012PASP..124..140B/abstract>`_).  Two magnitude units are also defined: ``VEGAmag`` and ``JMmag``.

Unit conversions between flux density and Vega-based magnitudes use the
`~astropy.units` `equivalency system
<http://docs.astropy.org/en/stable/units/equivalencies.html#unit-equivalencies>`_
and ``sbpy``'s :func:`~sbpy.units.spectral_density_vega`.
``spectral_density_vega`` requires a wavelength, frequency, or bandpass
for the transformation:

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

To use a bandpass, pass the name as a string (see
:func:`~synphot.spectrum.SpectralElement.from_filter`) or an instance
of `~synphot.spectrum.SpectralElement`:

  >>> from sbpy.units import VEGAmag, spectral_density_vega
  >>> m = 0.0 * VEGAmag
  >>> fluxd = m.to('erg/(cm2 s AA)',
  ...     spectral_density_vega('johnson_v'))  # doctest: +REMOTE_DATA
  >>> fluxd.value                   # doctest: +FLOAT_CMP +REMOTE_DATA
  3.588524658721229e-09


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
  >>> from sbpy.units import reflectance, VEGAmag, spectral_density_vega
  >>> mag = 3.4 * VEGAmag
  >>> radius = 460 * u.km
  >>> cross_sec = np.pi * (radius)**2
  >>> ref = mag.to('1/sr', reflectance(cross_section=cross_sec))
  >>> print('{0:.4f}'.format(ref))
  0.0287 1 / sr

Users can supply solar flux via keyword ``f_sun``, or solar magnitude in any
magnitude system via keyword ``M_sun`` as long as the quantity to be converted
is in the same unit as the supplied solar flux or magnitude, respectively:

  >>> ref = 0.0287 / u.sr
  >>> flux_ceres = mag.to('W/(m2 um)', spectral_density_vega('johnson_v')) # doctest: +IGNORE_OUTPUT
  >>> flux_sun = 1839.93273227 * u.Unit('W/(m2 um)')
  >>> cross_sec = flux_ceres.to('km2', reflectance(reflectance=ref,
  ...     f_sun=flux_sun))
  >>> radius = np.sqrt(cross_sec/np.pi)
  >>> print('{0:.2f}'.format(radius))
  459.69 km

`~sbpy.units.reflectance` can also support conversion between flux spectrum and
reflectance spectrum:

  >>> wave = [1.046, 1.179, 1.384, 1.739, 2.416] * u.um
  >>> flux = [1.636e-18, 1.157e-18, 8.523e-19, 5.262e-19, 1.9645e-19] \
  ...     * u.Unit('W/(m2 um)')
  >>> xsec = 0.0251 * u.km**2
  >>> ref = flux.to('1/sr', reflectance(cross_section=xsec, wfb=wave))
  >>> print(ref)  # doctest: +FLOAT_CMP
  [0.0021763  0.00201223 0.0022041  0.00269637 0.00292785] 1 / sr


Reference/API
-------------
.. automodapi:: sbpy.units
    :no-heading:
    :include-all-objects:
