
Spectroscopy Module (`sbpy.spectroscopy`)
=========================================

Introduction
------------

`~sbpy.spectroscopy` provides routines for the modeling and analysis of emission spectra of gas in comets and reflection spectra of asteroid surfaces and cometary comae.  Sub-module `sources` control `sbpy`'s photometric calibration using spectra of the Sun and Vega.


Spectral Gradients
------------------

Spectral gradient or slope is commonly expressed as a percent change
per wavelength interval, usually % per 100 nm or % per 0.1 μm.  The
class `~sbpy.photometry.SpectralGradient` enables easy conversion
between spectral gradient and color index (magnitudes), and
re-normalization to other wavelengths.

Initialize a spectral gradient object using ``astropy``'s `~astropy.units`:

  >>> import astropy.units as u
  >>> from sbpy.spectroscopy import SpectralGradient
  >>> from sbpy.units import hundred_nm
  >>> S = SpectralGradient(10 * u.percent / hundred_nm)
  >>> print(S)
  10.0 % / 100 nm

Initialize a spectral gradient from a color index:

  >>> w = (550, 650) * u.nm
  >>> SpectralGradient.from_color(w, 0.1 * u.mag)  # doctest: +FLOAT_CMP
  <SpectralGradient 9.20383492 % / 100 nm>

Convert spectral gradient (normalized to 550 nm) to a color index:

  >>> S = SpectralGradient(10 * u.percent / hundred_nm, wave0=550 * u.nm)
  >>> S.to_color((500, 600) * u.nm)  # doctest: +FLOAT_CMP
  <Quantity 0.10866423 mag>

Renormalize to 1.0 μm:

  >>> S = SpectralGradient(10 * u.percent / hundred_nm, wave0=550 * u.nm)
  >>> S.renormalize(1 * u.um)  # doctest: +FLOAT_CMP
  <SpectralGradient 6.89655172 % / 100 nm>



Reference/API
-------------
.. automodapi:: sbpy.spectroscopy
    :no-heading:
