Photometry Module (`sbpy.units`)
=====================================

Introduction
------------

`~sbpy.units` provides common planetary astronomy units and Vega-based
magnitude conversions.  ``sbpy`` units may be added to the top-level
``astropy.units`` namespace via:

  >>> import sbpy.units
  >>> sbpy.units.enable()

This will make them searchable via strings, such as
``u.Unit('mag(VEGA)')``.


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
  >>> m = 0 * vega_mag
  >>> fluxd = m.to('erg/(cm2 s AA)', spectral_density_vega(wave))
  >>> fluxd.value                                  # doctest: +FLOAT_CMP
  3.5469235179497687e-09
  >>> m = fluxd.to(VEGAmag, spectral_density_vega(wave))
  >>> m.value                                      # doctest: +FLOAT_CMP
  0.0

To use a bandpass, pass the name as a string (see
:func:`~synphot.spectrum.SpectralElement.from_filter`) or an instance
of `~synphot.spectrum.SpectralElement`:

  >>> from sbpy.units import vega_mag, spectral_density_vega
  >>> m = 0.0 * vega_mag
  >>> fluxd = m.to('erg/(cm2 s AA)',
  ...     spectral_density_vega('johnson_v'))  # doctest: +REMOTE_DATA
  >>> fluxd.value                   # doctest: +FLOAT_CMP +REMOTE_DATA
  3.588524658721229e-09


Reference/API
-------------
.. automodapi:: sbpy.units
    :no-heading:
