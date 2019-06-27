.. _sbpy_spectral_standards:

Spectroscopy Sources Module (`sbpy.spectroscopy.sources`)
=========================================================

Spectral standards and photometric calibration
----------------------------------------------
`sbpy`'s photometric calibration is partially based on spectra of the Sun and Vega.  There are built-in spectra for each, and users may provide their own spectra.

The spectrum of `Bohlin (2014) <https://dx.doi.org/10.1088/0004-6256/147/6/127>`_ is the default and only built-in spectrum for Vega.  It is distributed with `sbpy`.  Four solar spectra are built-in:

  * E490_2014 - E490 (2014) standard.
  * E490_2014LR - A low resolution version of the E490 standard.
  * Kurucz1993 - Kurucz (1993) model.
  * Castelli1996 - Castelli model from Colina et al. (1996).

The E490 spectra are included with `sbpy`, and the Kurucz and Castell spectra are downloaded as needed from `STScI's reference data system <http://www.stsci.edu/hst/observatory/crds/astronomical_catalogs.html>`_.

Each star has a class for use within `sbpy`.  The classes can be initialized with the default spectrum using :func:`~sbpy.spectroscopy.sun.Sun.from_default`:

.. doctest-requires:: synphot

  >>> from sbpy.spectroscopy import Sun
  >>> sun = Sun.from_default()
  >>> print(sun)
  <Sun: E490-00a (2014) reference solar spectrum (Table 3)>

The names of the built-in sources are stored as an internal array.  They can be discovered with :func:`~sbpy.spectroscopy.sun.Sun.show_builtin`, and used to initialize an object with :func:`~sbpy.spectroscopy.sun.Sun.from_builtin`:

.. doctest-requires:: synphot

  >>> from sbpy.spectroscopy import Sun
  >>> Sun.show_builtin()
      name                                description
  ------------ -----------------------------------------------------------------
     E490_2014                E490-00a (2014) reference solar spectrum (Table 3)
   E490_2014LR E490-00a (2014) low resolution reference solar spectrum (Table 4)
    Kurucz1993               Kurucz (1993) model, scaled by Colina et al. (1996)
  Castelli1996      Castelli model, scaled and presented by Colina et al. (1996)
  >>> sun = Sun.from_builtin('E490_2014LR')
  >>> print(sun)
  <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>

The solar spectrum in current use is controlled with an `astropy`  `~astropy.utils.state.ScienceState` named `~sbpy.spectroscopy.sun.default_sun`:

.. doctest-requires:: synphot

  >>> from sbpy.spectroscopy import Sun, default_sun
  >>> default_sun.set('E490_2014LR')
  <ScienceState default_sun: <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>>
  >>> # E490 low-resolution spectrum in effect for all of sbpy
  >>> sun = Sun.from_default()
  >>> print(sun)
  <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>

`default_sun` can also be used as a context manager to temporarily change the default spectrum:

.. doctest-requires:: synphot

  >>> from sbpy.spectroscopy import Sun, default_sun
  >>> default_sun.set('E490_2014')  # E490 in effect
  <ScienceState default_sun: <Sun: E490-00a (2014) reference solar spectrum (Table 3)>>
  >>> with default_sun.set('E490_2014LR'):
  ...   # E490 low-res in effect
  ...   print(Sun.from_default())
  <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>
  >>> # Back to module default.
  >>> print(Sun.from_default())
  <Sun: E490-00a (2014) reference solar spectrum (Table 3)>

Provide your own solar spectrum with the `Sun` class:

.. doctest-requires:: synphot

  >>> from sbpy.spectroscopy import Sun, default_sun
  >>> with default_sun.set(Sun.from_file('sun.txt')):  # doctest: +SKIP
  ...   # sun.txt in effect

See `~sbpy.spectroscopy.sun.Sun` for more information on ways to create solar spectra.

In a similar manner, Vega spectra are accessed and controlled via `~sbpy.spectroscopy.vega.Vega` and `~sbpy.spectroscopy.vega.default_vega`:

.. doctest-requires:: synphot

  >>> from sbpy.spectroscopy import Vega, default_vega
  >>> print(Vega.from_default())     # doctest: +SKIP
  <Vega: Dust-free template spectrum of Bohlin 2014>
  >>> with default_vega.set(Vega.from_file('vega.txt')):  # doctest: +SKIP
  ...   # vega.txt in effect


Observe the Sun through a filter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`sbpy` can simulate observations of comets and asteroids through filter bandpasses.  To support this functionality, the `Sun` and `Vega` classes have :func:`~sbpy.spectroscopy.SpectralStandard.filt` and :func:`sbpy.spectroscopy.SpectralStandard.observe` methods that return simulated flux densities.  `filt` is specific for filter bandpasses and returns effective wavelengths, `observe` is for general spectral rebinning (see next section).

Get the default solar spectrum, observe it through the Johnson V-band filter (distributed with `sbpy`), returning the result as a Vega-based magnitude in the Johnson-Morgan system:

.. doctest-requires:: synphot

  >>> from sbpy.spectroscopy import Sun
  >>> from sbpy.utils import get_bandpass
  >>> from sbpy.units import JMmag
  >>> sun = Sun.from_default()
  >>> bp = get_bandpass('Johnson V')
  >>> wave, mag = sun.filt(bp, unit=JMmag)
  >>> print('{}\n- Johnson V = {:.2f} at {:.0f}'.format(sun.description, mag, wave))   # doctest: +FLOAT_CMP
  E490-00a (2014) reference solar spectrum (Table 3)
   - Johnson V = -26.74 mag(JM) at 5502 Angstrom


Plot solar spectra
^^^^^^^^^^^^^^^^^^

Solar spectra in Sun objects can be plotted at the native resolution of the data, or rebinned. Plot the solar spectrum at the native resolution, and at a resolution of ~25:

.. doctest-requires:: synphot

  >>> import astropy.units as u
  >>> import numpy as np
  >>> import matplotlib.pyplot as plt
  >>> from sbpy.spectroscopy import Sun
  >>> # Create an array of wavelengths at R~25
  >>> wrange = 0.3, 0.8  # wavelength range
  >>> d = 1 + 1 / 25
  >>> n = int(np.ceil(np.log(wrange[1] / wrange[0]) / np.log(d)))
  >>> wave_binned = wrange[0] * d**np.arange(n) * u.um
  >>> # Get the default solar spectrum, and rebin it
  >>> sun = Sun.from_default()
  >>> fluxd_binned = sun.observe(wave_binned, unit='W / (m2 um)')
  >>> # Plot
  >>> plt.plot(sun.wave.to('um'), sun.fluxd.to('W/(m2 um)'),
  ...          ls='steps-mid', color='#1f77b4', label='Native resolution')
  >>> plt.plot(wave_binned, fluxd_binned, ls='steps-mid',
  ...          color='#ff7f0e', label='R~25')
  >>> plt.setp(plt.gca(), xlim=wrange, xlabel='Wavelength (μm)',
  ...          ylabel='Flux density (W/(m2 μm)')
  >>> plt.legend()
  >>> plt.tight_layout()

.. plot::

  import astropy.units as u
  import numpy as np
  import matplotlib.pyplot as plt
  from sbpy.spectroscopy.sun import Sun
  wrange = 0.3, 0.8  # wavelength range
  d = 1 + 1 / 25
  n = int(np.ceil(np.log(wrange[1] / wrange[0]) / np.log(d)))
  wave_binned = wrange[0] * d**np.arange(n) * u.um
  sun = Sun.from_default()
  fluxd_binned = sun.observe(wave_binned, unit='W / (m2 um)')
  plt.plot(sun.wave.to('um'), sun.fluxd.to('W/(m2 um)'), ls='steps-mid', color='#1f77b4', label='Native resolution')
  plt.plot(wave_binned, fluxd_binned, ls='steps-mid', color='#ff7f0e', label='R~25')
  plt.setp(plt.gca(), xlim=wrange, xlabel='Wavelength (μm)', ylabel='Flux density (W/(m2 μm)')
  plt.legend()
  plt.tight_layout()


Reference/API
-------------
.. automodapi:: sbpy.spectroscopy.sources
    :no-heading:
