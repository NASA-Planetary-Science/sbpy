
Spectroscopy Module (`sbpy.spectroscopy`)
=========================================

Introduction
------------

`~sbpy.spectroscopy` provides routines for working with emission and reflected/scattered light from asteroids and comets.


Spectral Gradients
------------------

Spectral gradient or slope is a measurement of the shape of a
spectrum.  It is commonly expressed as a percent change per wavelength
interval, e.g., % per 100 nm or % per 0.1 μm.  Equation 1 of A'Hearn et
al. (1984):

.. math::
   S = \frac{R(λ_1) - R(λ_0)}{R(λ_1) + R(λ_0)} \frac{2}{λ_1 - λ_0}

where *R(λ)* is the reflectivity at wavelength *λ*.

The class `~sbpy.spectroscopy.SpectralGradient` enables easy
conversion between spectral gradient and color index (magnitudes), and
re-normalization to other wavelengths.

`SpectralGradient` is a `~astropy.units.Quantity` with units of
inverse length.  For convenience, `sbpy` includes a
`~sbpy.units.hundred_nm` unit, which is equal to 100 nm:

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

Note we use the dimensionless magnitude unit from `astropy`, i.e., not
one that carries flux density units such as `astropy.units.ABmag`.

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
