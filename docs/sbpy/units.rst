Photometry Module (`sbpy.units`)
=====================================

Introduction
------------

`~sbpy.units` provides common planetary astronomy units and Vega-based
magnitude conversions.


Vega Magnitudes
---------------

With the `~synphot` package, `sbpy` has the capability to convert
between flux densities and the Vega-based magnitude system (`Bessell &
Murphy 2012
<https://ui.adsabs.harvard.edu/#abs/2012PASP..124..140B/abstract>`_).
The Vega-based system requires a spectrum of Vega for a standard.
`sbpy` includes the dust-free spectrum from `Bohlin 2014
<https://ui.adsabs.harvard.edu/#abs/2014AJ....147..127B/abstract>`_
and assumes it has a magnitude of 0.03 at all wavelengths (`Bessell &
Murphy 2012
<https://ui.adsabs.harvard.edu/#abs/2012PASP..124..140B/abstract>`_).
`synphot` defines a special unit, `VEGAMAG`, to represent this system
(imported into `sbpy.units` for convenience).

Unit conversion between flux density and `VEGAMAG` uses the
`~astropy.units` `equivalency system
<http://docs.astropy.org/en/stable/units/equivalencies.html#unit-equivalencies>`_
and `sbpy`'s :func:`~sbpy.units.spectral_density_vega`.
`spectral_density_vega` requires a wavelength, frequency, or bandpass
for the transformation:

  >>> import astropy.units as u
  >>> from sbpy.units import VEGAMAG, spectral_density_vega
  >>> m = 0.03 * VEGAMAG
  >>> fluxd = m.to('erg/(cm2 s AA)', spectral_density_vega(5500 * u.AA))
  >>> fluxd.value                                  # doctest: +FLOAT_CMP
  3.5469235179497687e-09
  >>> m = fluxd.to(VEGAMAG, spectral_density_vega(5500 * u.AA))
  >>> m.value                                      # doctest: +FLOAT_CMP
  0.03

To use a bandpass, pass the name as a string (see
:func:`~synphot.spectrum.SpectralElement.from_filter`) or an instance
of `~synphot.spectrum.SpectralElement`:

  >>> from sbpy.units import VEGAMAG, spectral_density_vega
  >>> m = 0.03 * VEGAMAG
  >>> fluxd = m.to('erg/(cm2 s AA)',
  ...     spectral_density_vega('johnson_v'))  # doctest: +REMOTE_DATA
  >>> fluxd.value                   # doctest: +FLOAT_CMP +REMOTE_DATA
  3.588524658721229e-09


Reference/API
-------------
.. automodapi:: sbpy.units
    :no-heading:
