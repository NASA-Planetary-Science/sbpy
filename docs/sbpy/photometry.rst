Photometry Module (`sbpy.photometry`)
=====================================

Introduction
------------

`~sbpy.photometry` provides routines for photometric modeling of small
bodies, and conversions between flux density and magnitude units.

Getting started
---------------

Magnitude conversions
^^^^^^^^^^^^^^^^^^^^^

Converting between flux density and magnitudes for any bandpass or
single wavelength is done with `convert_mag`:

  >>> import astropy.units as u
  >>> from sbpy.photometry import convert_unit
  >>> m0 = 0 * u.ABmag
  >>> fluxd = convert_mag(m0, 'erg/(s cm2 Hz)')
  >>> fluxd.value  # doctest: +FLOAT_CMP
  3.6307805477010035e-20
  >>> m = convert_mag(fluxd, u.ABmag)
  >>> m.value  # doctest: +FLOAT_CMP
  0
  

`sbpy` supports Vega-based (`Bessell & Murphy 2012
<https://ui.adsabs.harvard.edu/#abs/2012PASP..124..140B/abstract>`_),
AB (or AB:sub:`ν`; `Oke & Gunn 1983
<https://ui.adsabs.harvard.edu/#abs/1983ApJ...266..713O/abstract>`_),
and ST (or ST:sub:`λ`, AB:sub:`λ`) magnitude systems.  For a
comparison of these systems, see Bessell & Murphy (2012) or the
`synphot documentation
<https://synphot.readthedocs.io/en/latest/synphot/units.html#counts-and-magnitudes>`_.

The AB and ST systems are flux density based, and conversion to/from
magnitudes is straightforward.  The conversion constants are built
into `astropy`'s units module.

The Vega-based system requires a spectrum of Vega for a standard.
`sbpy` includes the dust-free spectrum from `Bohlin 2014
<https://ui.adsabs.harvard.edu/#abs/2014AJ....147..127B/abstract>`_
and assumes Vega has a magnitude of 0.03 at all wavelengths (`Bessell
& Murphy 2012
<https://ui.adsabs.harvard.edu/#abs/2012PASP..124..140B/abstract>`_).
A special unit, `VegaMag` is defined by `sbpy`.  Conversions to/from
the Vega system require the `synphot` package:

  >>> from sbpy.units import VegaMag
  >>> from sbpy.photometry import convert_mag
  >>> fluxd = convert_mag(0 * VegaMag, 'Jy', wave=5500 * u.AA)
  >>> fluxd.value  # doctest: +FLOAT_CMP
  3630.7805477010033

`synphot`'s `VEGAMAG` may also be used, but `convert_mag` will always
return `sbpy`'s `VegaMag`.

Reference/API
-------------
.. automodapi:: sbpy.photometry
    :no-heading:

.. automodapi:: sbpy.photometry.calibration
    :no-heading:
