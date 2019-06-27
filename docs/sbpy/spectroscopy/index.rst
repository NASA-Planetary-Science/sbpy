
Spectroscopy Module (`sbpy.spectroscopy`)
=========================================

Introduction
------------

`~sbpy.spectroscopy` provides routines for the modeling and analysis of emission spectra of gas in comets and reflection spectra of asteroid surfaces and cometary comae.  Sub-module `sources` control `sbpy`'s photometric calibration using spectra of the Sun and Vega.

JPLSpec Constants and Line Intensity Integral Conversion
--------------------------------------------------------

`sbpy.spectroscopy` has a function called ``molecular_data`` which takes care of
querying the JPL Molecular Spectral Catalog through the use of
`astroquery.jplspec` and calculates all the necessary constants needed for both
production rate and Einstein coefficient calculations. ``molecular_data``
returns an `sbpy.data.phys` instance with quantities in the order of: Transition
Frequency, Temperature, Integrated line intensity at 300 K, and Partition
function at 300 K. The function ``intensity_conversion`` takes care of
converting the intensity line integral at 300 K found in JPL Spec catalog and
convert it closer to the temperature given by the user.

.. code-block:: python
.. doctest-skip::

   >>> from sbpy.spectroscopy import molecular_data, intensity_conversion
   >>> import astropy.units as u
   >>> temp_estimate = 33. * u.K
   >>> vgas = 0.8 * u.km / u.s
   >>> mol_tag = 27001
   >>> transition_freq = 265.886434 * u.MHz
   >>> mol_data = molecular_data(temp_estimate, transition_freq, mol_tag, vgas)
   >>> intl = intensity_conversion(temp_estimate, transition_freq, mol_tag, vgas)
   >>> co
     [<Quantity 265886.18 MHz>,
      <Quantity 37.5 K>,
      <Quantity 0.08412014 MHz nm2>,
      424.32634842758375,
      53.91380704697852,
      1.7317,
      <Quantity 6.62607004e-34 J s>,
      <Quantity 1.38064852e-23 J / K>,
      <Quantity 2.99792458e+08 m / s>,
      <Quantity 800. m / s>,
      <Quantity 3.52359898e-22 J>,
      <Quantity 1.76181853e-22 J>]

.. doctest-skip::
   >>> intl
      <Quantity 3.99731047636546 MHz nm2>


Einstein coefficient
--------------------

`sbpy.spectroscopy` offers a function to calculate the Einstein coefficient
for a specific molecule and transition frequency.

.. code-block:: python
.. doctest-skip::

   >>> from sbpy.spectroscopy import einstein_coeff
   >>> e = einstein_coeff(temp_estimate, transition_freq, mol_tag, vgas)
   >>> e
      <Quantity 0.0008601294364222543 1 / s>


Production Rate Calculations
----------------------------

`sbpy.spectroscopy` provides several models to calculate production rates in comets.
One of the models followed by this module is based on the following paper:

| Drahus et al. September 2012. The Sources of HCN and CH3OH and the
| Rotational Temperature in Comet 103P/Hartley 2 from Time-resolved
| Millimeter Spectroscopy. The Astrophysical Journal, Volume 756,
| Issue 1.

This model does not include photodissociation. The model uses the functions for
constants, conversion, and Einstein coefficients as well as a JPLHorizons
lookup to obtain the distance to the comet at observation time. The following
example shows the usage of the module. For more information on the parameters
that are optional or needed for the module follow the link under Reference/API
section.

.. code-block:: python
.. doctest-skip::

  >>> from sbpy.spectroscopy import prodrate_np
  >>> from astropy.time import Time
  >>> from sbpy.data import Ephem
  >>> temp_estimate = 33. * u.K
  >>> target = '103P'
  >>> vgas = 0.8 * u.km / u.s
  >>> aper = 30 * u.m
  >>> b = 1.13
  >>> mol_tag = 27001
  >>> transition_freq = 265.886434 * u.MHz
  >>> spectra = 1.22 * u.K * u.km / u.s
  >>> time = Time('2010-11-3 00:48:06', format='iso')
  >>> ephemobj = Ephem(target, epochs=time.jd, id_type='id')
  >>> q = prodrate_np(spectra, temp_estimate, transition_freq,
  ...                 mol_tag, ephemobj, vgas, aper, b=b)
  >>> q
  <Quantity 1.0432591198553935e+25 1 / s>

Another model included in the module is based off of the model in the following
literature:

| Haser 1957, Bulletin de la Societe Royale des Sciences de Liege 43, 740.
| Newburn and Johnson 1978, Icarus 35, 360-368.

This model takes in an initial guess for the production rate, and uses the
module found in ``sbpy.activity.gas`` to find a ratio between the model model
total number of molecules and the number of molecules calculated from the data
to scale the model Q and output the new production rate from the result. This
model does account for the effects of photolysis.

.. code-block:: python
.. doctest-skip::
   
  >>> from sbpy.activity.gas import Haser
  >>> from astropy.time import Time
  >>> from sbpy.data import Ephem
  >>> Q = spec.production_rate(coma, molecule='H2O')
  >>> coma = Haser(Q, v, parent)
  >>>
  >>> Q_estimate = 2.8*10**(28) / u.s
  >>> transition_freq = (230.53799 * u.GHz).to('MHz')
  >>> aper = 10 * u.m
  >>> mol_tag = 28001
  >>> temp_estimate = 25. * u.K
  >>> target = 'C/2016 R2'
  >>> b = 0.74
  >>> vgas = 0.5 * u.km / u.s
  >>>
  >>> time = '2017-12-22 05:24:20'
  >>> ephemobj = Ephem(target, epochs=time.jd)
  >>> spectra = 0.26 * u.K * u.km / u.s
  >>>
  >>> parent = photo_timescale('CO') * vgas
  >>>
  >>> coma = Haser(Q_estimate, vgas, parent)
  >>>
  >>> Q = spec.production_rate(coma, spectra, temp_estimate,
                               transition_freq, mol_tag, ephemobj,
                               aper=aper, b=b)
  >>>
  >>> print(Q)
      <Quantity [1.64403219e+28] 1 / s>


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
