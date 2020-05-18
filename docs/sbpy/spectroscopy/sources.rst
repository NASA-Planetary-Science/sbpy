Spectroscopy Sources Module (`sbpy.spectroscopy.sources`)
=========================================================

`sbpy` uses the `~sbpy.spectroscopy.sources` module for evaluation and
binning of spectra and spectral models.

Blackbody spectra
-----------------

The `~sbpy.spectroscopy.sources.BlackbodySource` encapsulates a
blackbody spectral model.  It represents the spectral flux density
from a sphere, i.e., *πB(T)*, where *B* is the Planck function.

Initialize from temperature:

.. doctest-requires:: synphot

  >>> import astropy.units as u
  >>> from sbpy.spectroscopy.sources import BlackbodySource
  >>>
  >>> B = BlackbodySource(T=278 * u.K)
  >>> print(B)
  <BlackbodySource: T=278.0 K>

Observe the source through a fictitious 10-μm filter:

.. doctest-requires:: synphot

  >>> from synphot import SpectralElement, Box1D
  >>> N = SpectralElement(Box1D, x_0=10.5 * u.um, width=1 * u.um)
  >>> print(B.observe(N, unit='MJy'))  # doctest: +FLOAT_CMP
  783577717.2783749 MJy

Observe the source through a low-resolution spectrometer:

.. doctest-requires:: synphot

  >>> import numpy as np
  >>> import matplotlib.pyplot as plt
  >>>
  >>> wave = np.logspace(0.5, 1.5, 100) * u.um
  >>> fluxd = B.observe(wave, unit='MJy')
  >>>
  >>> plt.plot(wave, fluxd, drawstyle='steps-mid', label=str(B.T))
  ...     # doctest: +IGNORE_OUTPUT
  >>> plt.setp(plt.gca(), xlabel='Wavelength (μm)', ylabel='$F_ν$ (MJy)')
  ...     # doctest: +IGNORE_OUTPUT

.. plot::

   import numpy as np
   import matplotlib.pyplot as plt
   import astropy.units as u
   from sbpy.spectroscopy.sources import BlackbodySource

   B = BlackbodySource(T=278 * u.K)
   wave = np.logspace(0.5, 1.5, 100) * u.um
   fluxd = B.observe(wave, unit='MJy')
  
   plt.plot(wave, fluxd, drawstyle='steps-mid', label=str(B.T))
   plt.setp(plt.gca(), xlabel='Wavelength (μm)', ylabel='$F_ν$ (MJy)')


Developers' Notes
-----------------

Developers wanting to use the spectrophotometric mechanics should
create classes that inherit from
`sbpy.spectroscopy.sources.SpectralSource` using
`~sbpy.spectroscopy.sources.BlackbodySource` as an example.


Reference/API
-------------
.. automodapi:: sbpy.spectroscopy.sources
    :no-heading:
