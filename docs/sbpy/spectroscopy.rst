Spectroscopy Module (`sbpy.spectroscopy`)
=========================================

Introduction
------------

`sbpy.spectroscopy` provides routines for the modeling and analysis of emission spectra of gas in comets and reflection spectra of asteroid surfaces.  Sub-modules `sun` and `vega` control `sbpy`'s photometric calibration.

Spectral standards and photometric calibration
----------------------------------------------
`sbpy` includes a set of solar spectra and a spectrum of Vega used for the photometric calibration of the module.  It also includes support for users to provide their own spectra.  The solar spectra included are::

  * E490_2014 - E490 (2014) standard.
  * E490_2014LR - A low resolution version of the E490 standard.
  * Kurucz1993 - Kurucz (1993) model.
  * Castelli1996 - Castelli model from Colina et al. (1996).

The E490 spectra are included with `sbpy`, and the Kurucz and Castell spectra are downloaded as needed from `STScI's reference data <http://www.stsci.edu/hst/observatory/crds/astronomical_catalogs.html>`_. The spectrum in use is controlled and revealed through `default_sun`::
  
  >>> from sbpy.spectroscopy.sun import default_sun
  >>> default_sun.set('E490_2014')
  >>> # E490 in effect for all of sbpy
  >>> sun = default_sun.get() # Get the spectrum as a `Sun` class

`default_sun` can be used as a context manager controls to temporarily change the default spectrum::

  >>> from sbpy.spectroscopy.sun import default_sun
  >>> with default_sun.set('E490_2014LR'):
  ...   # E490 low-res in effect
  ...   sun = default_sun.get()
  >>> sun = default_sun.get()  # Back to module default.

Provide your own solar spectrum with the `Sun` class::

  >>> from sbpy.spectroscopy.sun import Sun, default_sun
  >>> with default_sun.set(Sun.from_file('sun.txt')):  # doctest: +SKIP
  ...   # vega.txt in effect

See `Sun` for more information on ways to create solar spectra.

In a similar manner, Vega spectra are controlled via `default_vega` and `Vega`::

  >>> from sbpy.spectroscopy.vega import Vega, default_vega
  >>> with default_vega.set(Vega.from_file('vega.txt')):  # doctest: +SKIP
  ...   # vega.txt in effect
  
Reference/API
-------------
.. automodapi:: sbpy.spectroscopy
    :no-heading:
