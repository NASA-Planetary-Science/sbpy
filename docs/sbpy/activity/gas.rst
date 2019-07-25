Gas Comae (`sbpy.activity.gas`)
===============================

Photolysis
----------

Two functions provide reference data for the photolysis of gas molecules in optically thin comae: :func:`~sbpy.activity.gas.photo_lengthscale` and :func:`~sbpy.activity.gas.photo_timescale`.  The source data is stored in `sbpy.activity.gas.data`.

:func:`~sbpy.activity.gas.photo_lengthscale` provides empirical comae lengthscales (defaults to Cochran and Schleicher 1993)):

  >>> from sbpy.activity import gas
  >>> gas.photo_lengthscale(None)
  Traceback (most recent call last):
      ...
  ValueError: Invalid species None.  Choose from:
  H2O [CS93]
  OH [CS93]
  >>> gas.photo_lengthscale('H2O')  # doctest: +FLOAT_CMP
  <Quantity 24000. km>

Use :func:`~sbpy.activity.gas.photo_timescale` to retrieve photolysis timescales:

  >>> gas.photo_timescale(None)
  Traceback (most recent call last):
      ...
  ValueError: Invalid species None.  Choose from:
  CH3OH [C94]
  CN [H92]
  CO [CE83]
  CO2 [CE83]
  H2CO [C94]
  H2O [CS93]
  HCN [C94]
  OH [CS93]
  >>> gas.photo_timescale('H2O')
  <Quantity 52000. s>

Some sources provide values for the quiet and active Sun (Huebner et al. 1992):

  >>> gas.photo_timescale('CN', source='H92')
  <Quantity [315000., 135000.] s>


With the :doc:`../bib`, the citation may be discovered:

  >>> from sbpy import bib
  >>> with bib.Tracking():
  ...    tau = gas.photo_timescale('H2O')
  >>> print(bib.to_text())    # doctest: +REMOTE_DATA
  sbpy.activity.gas.core.photo_timescale:
    H2O photodissociation timescale:
        Cochran & Schleicher 1993, Icarus, Vol 105, 1, 235


Fluorescence
------------
Reference data for fluorescence band emission is available via :func:`~sbpy.activity.gas.fluorescence_band_strength`.  Compute the fluorescence band strength (luminosity per molecule) of the OH 0-0 band at 1 au from the Sun, moving towards the Sun at 1 km/s (defaults to Schleicher and A'Hearn 1988):

  >>> import astropy.units as u
  >>> LN = gas.fluorescence_band_strength('OH 0-0', -1 * u.km / u.s)
  >>> print(LN)    # doctest: +FLOAT_CMP
  [1.54e-15] erg / s
  

Gas coma models
---------------

The Haser (1957) model for parent and daughter species is included, with some calculation enhancements based on Newburn and Johnson (1978).  With `~sbpy.activity.gas.Haser`, we may compute the column density and total number of molecules within aperture:

  >>> Q = 1e28 / u.s        # production rate
  >>> v = 0.8 * u.km / u.s  # expansion speed
  >>> parent = gas.photo_lengthscale('H2O')
  >>> daughter = gas.photo_lengthscale('OH')
  >>> coma = gas.Haser(Q, v, parent, daughter)
  >>> print(coma.column_density(10 * u.km))    # doctest: +FLOAT_CMP
  7.099280153851781e+17 1 / m2
  >>> print(coma.total_number(1000 * u.km))    # doctest: +FLOAT_CMP
  1.161357452192558e+30

The gas coma models work with sbpy's apertures:

  >>> from sbpy.activity import AnnularAperture
  >>> ap = AnnularAperture((5000, 10000) * u.km)
  >>> print(coma.total_number(ap))    # doctest: +FLOAT_CMP
  3.8133654170856037e+31


Production Rate calculations
----------------------------

`~sbpy.activity.gas.productionrate` offers various functions that aid in the calculation
of production rates. `~sbpy.data.phys` has a function called `~sbpy.data.Phys.from_jplspec` which takes care of querying the JPL Molecular Spectral Catalog through the use of `~astroquery.jplspec` and calculates all the necessary constants needed for
production rate calculations in this module. Yet, the option for the user to
provide their own molecular data is possible through the use of an `~sbpy.data.phys` object, as long as it has the required information. It is imperative to read
the documentation of the functions in this section to understand what is needed
for each. If the user does not have the necessary data, they can build an object
using JPLSpec:

.. doctest-skip::

    >>> from sbpy.data.phys import Phys
    >>> import astropy.units as u
    >>> temp_estimate = 47. * u.K
    >>> transition_freq = (230.53799 * u.GHz).to('MHz')
    >>> mol_tag = '^CO$'
    >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)

Having this information, we can move forward towards the calculation of
production rate. The functions that sbpy currently provides to calculate
production rates are listed below.

Integrated Line Intensity Conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The JPL Molecular Spectroscopy Catalog offers the integrated line intensity
at 300 K for a molecule. Yet, in order to calculate production rate, we need
to know the integrated line intensity at a given temperature. This function
takes care of converting the integrated line intensity at 300 K to its equivalent
in the desired temperature using equations provided by the JPLSpec documentation.
For more information on the needed parameters for this function follow the link
for `~sbpy.activity.intensity_conversion` under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import intensity_conversion
    >>> intl = intensity_conversion(mol_data)
    >>> mol_data.add_column([intl.value] * intl.unit, name='intl')
     11
    >>> intl
     <Quantity 0.00280051 MHz nm2>


Einstein Coefficient Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike catalogs like LAMDA, JPLSpec does not offer the Eistein coefficient
and it must be calculated using equations provided by the JPL Molecular
Spectroscopy Catalog. These equations have been compared to established LAMDA
values of the Einstein Coefficient for HCN and CO, and no more than
a 24% difference has been found between the calculation from JPLSpec and the
LAMDA catalog value. Since JPLSpec and LAMDA are two very different
catalogs with different data, the difference is expected, and the user is
allowed to provide their own Einstein Coefficient if they want. If the user
does want to provide their own Einstein Coefficient, they may do so simply
by appending their value with the unit 1/s to the `~sbpy.data.Phys` object, called
`mol_data` in these examples. For more information on the needed parameters
for this function follow the link for `~sbpy.activity.einstein_coeff`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import einstein_coeff
    >>> au = einstein_coeff(mol_data)
    >>> mol_data.add_column([au.value] * au.unit, name = 'Einstein Coefficient')
     12
    >>> au
      <Quantity 7.03946054e-07 1 / s>


Beta Factor Calculation
^^^^^^^^^^^^^^^^^^^^^^^

Returns beta factor based on timescales from `~sbpy.activity.gas` and distance
from the Sun using an `~sbpy.data.ephem` object. The calculation is
parent photodissociation timescale * (distance from comet to Sun)**2
and it accounts for certain photodissociation and geometric factors needed
in the calculation of total number of molecules `~sbpy.activity.total_number_nocd`
If you wish to provide your own beta factor, you can calculate the equation
expressed in units of AU**2 * s , all that is needed is the timescale
of the molecule and the distance of the comet from the Sun. Once you
have the beta factor you can append it to your mol_data phys object
with the name 'beta' or any of its alternative names. For more information on
the needed parameters for this function follow the link for
`~sbpy.activity.beta_factor` under Reference/API section.

.. doctest-skip::

    >>> from astropy.time import Time
    >>> from sbpy.data import Ephem
    >>> from sbpy.activity import beta_factor
    >>> target = 'C/2016 R2'
    >>> time = Time('2017-12-22 05:24:20', format = 'iso')
    >>> ephemobj = Ephem.from_horizons(target, epochs=time.jd)
    >>> beta = beta_factor(mol_data, ephemobj)
    >>> mol_data.add_column([beta.value] * beta.unit, name='beta')
     13
    >>> beta
     <Quantity [13333365.25745597] AU2 s>


Total Number
^^^^^^^^^^^^

In order to obtain our total number of molecules from flux data,
with no column density information, we can use equation 10 from `Bockelee-Morvan
et al. 2004 <https://ui.adsabs.harvard.edu/#abs/2004come.book..391B>`_ and the
millimeter/submillimeter spectroscopy beam factors explained and detailed
in equation 1.3 from:

    | Drahus, M. (2010). Microwave observations and modeling of the molecular
    | coma in comets. PhD Thesis, Georg-August-Universität Göttingen.

If the user prefers to give the total number, they may do so by appending
to the mol_data `~sbpy.data.phys` object with the name `total_number_nocd` or
any of its alternative names. For more information on the needed parameters
for this function follow the link for `~sbpy.activity.total_number_nocd`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import total_number_nocd
    >>> integrated_flux = 0.26 * u.K * u.km / u.s
    >>> b = 0.74
    >>> aper = 10 * u.m
    >>> tnum = total_number_nocd(integrated_flux, mol_data, aper, b)
    >>> mol_data.add_column([tnum], name='total_number_nocd')
     14
    >>> tnum
     <Quantity [2.93988826e+26]>


Simplified Model for Production Rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`~sbpy.activity` provides several models to calculate production rates in comets.
One of the models followed by this module is based on the following paper:

| Drahus et al. September 2012. The Sources of HCN and CH3OH and the
| Rotational Temperature in Comet 103P/Hartley 2 from Time-resolved
| Millimeter Spectroscopy. The Astrophysical Journal, Volume 756,
| Issue 1.

The following example shows the usage of the function. This LTE model does not
include photodissociation, but it does serve as way to obtain educated
first guesses for other models within sbpy. For more information on the
parameters that are needed for the function follow
the link for the function `from_Drahus` in `~sbpy.activity.LTE`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import LTE
    >>> vgas = 0.5 * u.km / u.s
    >>> lte = LTE()
    >>> q = lte.from_Drahus(integrated_flux, mol_data, ephemobj, vgas, aper, b=b)
    >>> q
     <Quantity 3.59397119e+28 1 / s>


Haser Model for Production Rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another model included in the module is based off of the model in the following
literature:

| Haser 1957, Bulletin de la Societe Royale des Sciences de Liege 43, 740.
| Newburn and Johnson 1978, Icarus 35, 360-368.

This model is well-known as the Haser model. In the case of our implementation
the function takes in an initial guess for the production rate, and uses the
module found in `~sbpy.activity.gas` to find a ratio between the model
total number of molecules and the number of molecules calculated from the data
to scale the model Q and output the new production rate from the result. This
LTE model does account for the effects of photolysis. For more information
on the parameters that are needed for the function follow
the link for the function `from_Haser` in `~sbpy.activity.LTE`
under Reference/API section.

.. doctest-skip::

    >>> from sbpy.activity import Haser, photo_timescale
    >>> Q_estimate = 3.5939*10**(28) / u.s
    >>> parent = photo_timescale('CO') * vgas
    >>> coma = Haser(Q_estimate, vgas, parent)
    >>> lte = LTE()
    >>> Q = lte.from_Haser(coma, mol_data, aper=aper)
    >>> Q
     <Quantity [[9.35795579e+27]] 1 / s>


Reference/API
-------------
.. automodapi:: sbpy.activity.gas
    :no-heading:
