
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
  10.0 % / (100 nm)

Initialize a spectral gradient from a color index:

.. doctest-requires:: synphot

  >>> w = (550, 650) * u.nm
  >>> SpectralGradient.from_color(w, 0.1 * u.mag)  # doctest: +FLOAT_CMP
  <SpectralGradient 9.20383492 % / (100 nm)>

Note we use the dimensionless magnitude unit from `astropy`, i.e., not
one that carries flux density units such as `astropy.units.ABmag`.

Convert spectral gradient (normalized to 550 nm) to a color index:

.. doctest-requires:: synphot

  >>> S = SpectralGradient(10 * u.percent / hundred_nm, wave0=550 * u.nm)
  >>> S.to_color((500, 600) * u.nm)  # doctest: +FLOAT_CMP
  <Quantity 0.10866423 mag>

Renormalize to 1.0 μm:

  >>> S = SpectralGradient(10 * u.percent / hundred_nm, wave0=550 * u.nm)
  >>> S.renormalize(1 * u.um)  # doctest: +FLOAT_CMP
  <SpectralGradient 6.89655172 % / (100 nm)>


Spectral Reddening
------------------

Linear spectral reddening is enabled by the class `~sbpy.spectroscopy.Reddening`,
which is based on `~synphot.spectrum.BaseUnitlessSpectrum`.

Initialize a `~sbpy.spectroscopy.Reddening` class from a spectral gradient:


.. plot::
  :include-source: True

  import matplotlib.pyplot as plt
  import astropy.units as u
  from sbpy.spectroscopy import SpectralGradient, Reddening
  from sbpy.units import hundred_nm

  S = SpectralGradient(10 * u.percent / hundred_nm, wave0=0.55 * u.um)
  linear_reddening = Reddening(S)
  linear_reddening.plot()
  plt.gca().grid()


This class can then be used to linearly redden a spectrum as a
`~synphot.SourceSpectrum` class instance:


.. plot::
  :include-source: True

  import numpy as np
  import matplotlib.pyplot as plt
  import astropy.units as u
  from synphot import SourceSpectrum
  from synphot.models import BlackBodyNorm1D
  from sbpy.spectroscopy import SpectralGradient, Reddening
  from sbpy.units import hundred_nm

  S = SpectralGradient(10 * u.percent / hundred_nm, wave0=0.55 * u.um)
  linear_reddening = Reddening(S)
  spec = SourceSpectrum(BlackBodyNorm1D, temperature=5500 * u.K)
  reddened = spec * linear_reddening
  wv = np.linspace(0.3, 1, 100) * u.um

  plt.plot(wv, spec(wv))
  plt.plot(wv, reddened(wv))
  plt.legend(['Original', 'Reddened'])
  plt.setp(plt.gca(), xlabel='Wavelength (um)', ylabel='Flux (PHOTLAM)')
  plt.tight_layout()

``sbpy`` ``SpectralSource`` objects have a ``redden()`` method for reddening.
The following example reddens a solar spectrum:


.. plot::
  :include-source: True

  import numpy as np
  import matplotlib.pyplot as plt
  import astropy.units as u
  from sbpy.calib import Sun
  from sbpy.spectroscopy import SpectralGradient
  from sbpy.units import hundred_nm

  sun = Sun.from_builtin('E490_2014LR') # low-resolution solar spectrum
  S = SpectralGradient(10 * u.percent / hundred_nm, wave0=0.55 * u.um)
  red_sun = sun.redden(S)

  wv = np.linspace(0.3, 1, 100) * u.um
  fluxd_unit = u.Unit('W/(m2 um)')

  plt.plot(wv, sun.observe(wv, unit=fluxd_unit))
  plt.plot(wv, red_sun.observe(wv, unit=fluxd_unit))
  plt.legend(['Original', 'Reddened'])
  plt.setp(plt.gca(), xlabel='Wavelength (um)',
      ylabel='Flux density ({})'.format(fluxd_unit))
  plt.tight_layout()


Reference/API
-------------
.. automodapi:: sbpy.spectroscopy
    :no-heading:
