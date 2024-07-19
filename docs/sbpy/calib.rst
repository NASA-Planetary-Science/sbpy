.. _sbpy-calib:

Spectral Standards and Photometric Calibration (`sbpy.calib`)
=============================================================

sbpy's photometric calibration is based on spectra of the Sun and Vega.  For example, they are used to convert between :ref:`reflectance, cross-section, and magnitude <reflectance-equivalencies>`, between :ref:`Afρ and spectral flux density <afrho-to-from-flux-density>`, and between :ref:`Vega-based and other magnitude systems <vega-magnitudes>`.   sbpy has built-in spectra for each, and users may provide their own.

The spectrum of `Bohlin (2014) <https://dx.doi.org/10.1088/0004-6256/147/6/127>`_ is the default and only built-in spectrum for Vega.  It is distributed with sbpy.  Five solar spectra are built-in:

  * Castelli1996 - Castelli model from Colina et al. (1996).
  * E490_2014 - E490 (2014) standard.
  * E490_2014LR - A low resolution version of the E490 standard.
  * Kurucz1993 - Kurucz (1993) model.
  * calspec - R=5000, created by R. Bohlin from Kurucz Special Model

The E490 spectra are included with sbpy, and the others are downloaded as needed from MAST's `Spectral Atlas Files for Synphot Software (REFERENCE-ATLASES) <https://archive.stsci.edu/hlsp/reference-atlases>`_ or STScI's `CALSPEC Database <https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec>`_.

Each star has a class for use within sbpy.  The classes can be initialized with the default spectrum using :func:`~sbpy.calib.SpectralStandard.from_default`:

.. doctest-requires:: synphot

    >>> from sbpy.calib import Sun
    >>> sun = Sun.from_default()
    >>> print(sun)
    <Sun: E490-00a (2014) reference solar spectrum (Table 3)>

The names of the built-in sources are stored as an internal array.  They can be discovered with :func:`~sbpy.calib.SpectralStandard.show_builtin`, and used to initialize an object with :func:`~sbpy.calib.SpectralStandard.from_builtin`:

.. doctest-requires:: synphot

    >>> from sbpy.calib import Sun
    >>> Sun.show_builtin()
        name                                description
    ------------ -----------------------------------------------------------------
    Castelli1996      Castelli model, scaled and presented by Colina et al. (1996)
       E490_2014                E490-00a (2014) reference solar spectrum (Table 3)
     E490_2014LR E490-00a (2014) low resolution reference solar spectrum (Table 4)
      Kurucz1993               Kurucz (1993) model, scaled by Colina et al. (1996)
         calspec            R=5000, created by R. Bohlin from Kurucz Special Model
    >>> sun = Sun.from_builtin('E490_2014LR')
    >>> print(sun)
    <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>

Controlling the default spectra
-------------------------------

The Vega and solar spectra in current use are respectively controlled with `~astropy.utils.state.ScienceState` objects named `~sbpy.calib.vega_spectrum` and `~sbpy.calib.solar_spectrum`:

.. doctest-requires:: synphot

    >>> from sbpy.calib import Sun, solar_spectrum
    >>> solar_spectrum.set('E490_2014LR')
    <ScienceState solar_spectrum: <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>>
    >>> # E490 low-resolution spectrum in effect for all of sbpy
    >>> sun = Sun.from_default()
    >>> print(sun)
    <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>

`~sbpy.calib.vega_spectrum` and `~sbpy.calib.solar_spectrum` can also be used as a context manager to temporarily change the default spectrum:

.. doctest-requires:: synphot

    >>> from sbpy.calib import Sun, solar_spectrum
    >>> solar_spectrum.set('E490_2014')  # E490 in effect
    <ScienceState solar_spectrum: <Sun: E490-00a (2014) reference solar spectrum (Table 3)>>
    >>> with solar_spectrum.set('E490_2014LR'):
    ...   # E490 low-res in effect
    ...   print(Sun.from_default())
    <Sun: E490-00a (2014) low resolution reference solar spectrum (Table 4)>
    >>> # Back to module default.
    >>> print(Sun.from_default())
    <Sun: E490-00a (2014) reference solar spectrum (Table 3)>

Provide your own solar spectrum with the `~sbpy.calib.Sun` class:

.. doctest-requires:: synphot

    >>> from sbpy.calib import Sun, solar_spectrum
    >>> with solar_spectrum.set(Sun.from_file('sun.txt')):  # doctest: +SKIP
    ...   # sun.txt in effect

See `~sbpy.calib.Sun` for more information on ways to create solar spectra.

An example showing how to change the default Vega spectrum:

.. doctest-requires:: synphot

    >>> from sbpy.calib import Vega, vega_spectrum
    >>> print(Vega.from_default())     # doctest: +SKIP
    <Vega: Dust-free template spectrum of Bohlin 2014>
    >>> with vega_spectrum.set(Vega.from_file('vega.txt')):  # doctest: +SKIP
    ...   # vega.txt in effect


Photometric calibration (without spectra or `synphot`)
------------------------------------------------------

The `~astropy.utils.state.ScienceState` objects `~sbpy.calib.solar_fluxd` and `~sbpy.calib.vega_fluxd` control photometric calibration by filter name.  These are completely independent of the spectroscopic calibration and can be used without the optional `synphot` package.  The spectral flux densities (per unit wavelength) of the Sun and Vega are provided and enabled by default.  Values and filters are from Willmer (2018):

    >>> from sbpy.calib import Sun, solar_fluxd, vega_fluxd
    >>> import sbpy.units as sbu
    >>>
    >>> solar_fluxd.set('Willmer2018')   # doctest: +IGNORE_OUTPUT
    >>> sun = Sun(None)
    >>> print(sun.observe('PS1 r'))    # doctest: +FLOAT_CMP
    167.49428760264365 erg / (Angstrom s cm2)
    >>> vega_fluxd.set('Willmer2018')   # doctest: +IGNORE_OUTPUT
    >>> print(sun.observe('PS1 r', unit=sbu.VEGAmag))    # doctest: +FLOAT_CMP
    -27.05 mag(VEGA)

Use ``solar_fluxd.get('Willmer2018')`` to discover all built-in values.

Users wanting to calibrate data with their own flux densities may do so.  For example, set the *V*-band apparent magnitude of the Sun to that in Colina et al. (1996).  Observations through the ``'V'`` filter will use the specified value:

    >>> solar_fluxd.set({'V': -26.75 * sbu.VEGAmag})  # doctest: +IGNORE_OUTPUT
    >>> sun = Sun.from_default()
    >>> print(sun.observe('V'))
    -26.75 mag(VEGA)

Some sbpy calculations will require the effective wavelength or the pivot wavelength.  These are optional parameters that may be specified with `~sbpy.calib.solar_fluxd` and `~sbpy.calib.vega_fluxd`:

    >>> import astropy.units as u
    >>> from sbpy.calib import vega_fluxd, Vega
    >>>
    >>> # values from Willmer (2018)
    >>> vega_fluxd.set({
    ...     'V': 3674.73 * u.Jy,
    ...     'V(lambda eff)': 5476 * u.AA
    ... })    # doctest: +IGNORE_OUTPUT
    <ScienceState vega_fluxd: {'V': <Quantity 3674.73 Jy>, 'V(lambda eff)': <Quantity 5476. Angstrom>}>
    >>> vega = Vega.from_default()
    >>> print(vega.observe('V'))
    3674.73 Jy
    >>> print(vega.observe('V', unit='erg/(s cm2 AA)'))
    ...     # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    UnitConversionError: 'Jy' (spectral flux density) and 'erg / (Angstrom s cm2)' (spectral flux density wav) are not convertible  Is "V(lambda pivot)" required and was it provided?
    >>> vega_fluxd.set({
    ...     'V': 3674.73 * u.Jy,
    ...     'V(lambda eff)': 5476 * u.AA,
    ...     'V(lambda pivot)': 5511 * u.AA
    ... })    # doctest: +IGNORE_OUTPUT
    <ScienceState vega_fluxd: {'V': <Quantity 3674.73 Jy>, 'V(lambda eff)': <Quantity 5476. Angstrom>, 'V(lambda pivot)': <Quantity 5511. Angstrom>}>
    >>> print(vega.observe('V', unit='erg/(s cm2 AA)'))   # doctest: +FLOAT_CMP
    3.62701e-9 erg / (Angstrom s cm2)

Observe the Sun
---------------

sbpy can simulate observations of comets and asteroids through spectrometers and filter bandpasses.  To support this functionality, the `~sbpy.calib.Sun` and `~sbpy.calib.Vega` classes have the :func:`~sbpy.calib.SpectralStandard.observe` method that returns simulated flux densities.  Users may request observations through filter bandpasses, or at a set of wavelengths (the default is to rebin the source spectrum).

Get the default solar spectrum, observe it through the Johnson V-band filter (distributed with sbpy), returning the result as a Vega-based magnitude in the Johnson-Morgan system:

.. doctest-requires:: synphot

    >>> from sbpy.calib import Sun
    >>> from sbpy.photometry import bandpass
    >>> from sbpy.units import JMmag
    >>>
    >>> sun = Sun.from_default()
    >>> bp = bandpass('Johnson V')
    >>> fluxd = sun.observe(bp, unit=JMmag)
    >>> print(fluxd)    # doctest: +FLOAT_CMP
    -26.744715028702647 mag(JM)


Binning versus interpolation with ``observe()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a user requests a series of wavelengths or frequencies with a `~astropy.units.Quantity` object, the default for :func:`~sbpy.calib.SpectralStandard.observe` is to rebin the source spectrum using the requested values as bin centers.  This behavior is appropriate when the source spectrum is at a higher spectral resolution than the requested wavelengths.  This is because `~synphot` assumes source spectra are continuous functions, rather than observations taken through a spectrometer (binned data).

When the requested spectral resolution is comparable to the spectral resolution of the source, rebinning may result in errors at the percent-level or more.  Instead, use the ``interpolate=True`` parameter for ``observe``.

Compare interpolation and rebinning for the E490 low-resolution solar spectrum, using the stored wavelengths of the spectrum.  Initialize a `~sbpy.calib.Sun` object with the low-resolution spectrum.

.. doctest-requires:: synphot

    >>> import numpy as np
    >>> from sbpy.calib import Sun
    >>> sun = Sun.from_builtin('E490_2014LR')

Inspect a sub-set of the data for this example.

.. doctest-requires:: synphot

    >>> wave = sun.wave[430:435]
    >>> S = sun.fluxd[430:435]
    >>> print(wave)    # doctest: +FLOAT_CMP
    [5495. 5505. 5515. 5525. 5535.] Angstrom
    >>> print(S)       # doctest: +FLOAT_CMP
    [1895. 1862. 1871. 1846. 1882.] W / (um m2)

Interpolate with observe() and compare to the original values.

Re-bin with observe using the same wavelengths as band centers.

.. doctest-requires:: synphot

    >>> S_rebin = sun.observe(wave)
    >>> np.allclose(S.value, S_rebin.value)
    False

Inspect the differences.

.. doctest-requires:: synphot

    >>> print((S_rebin - S) / (S_rebin + S) * 2)    # doctest: +FLOAT_CMP
    [-0.00429693  0.00281266 -0.00227604  0.00412338 -0.00132301]


Plot solar spectra
^^^^^^^^^^^^^^^^^^

Solar spectra in Sun objects can be plotted at the native resolution of the data, or rebinned. Plot the solar spectrum at the native resolution, and at a resolution of ~25:

.. doctest-skip::

    >>> import astropy.units as u
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from sbpy.calib import Sun
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
    ...          drawstyle='steps-mid', color='#1f77b4',
    ...          label='Native resolution')
    >>> plt.plot(wave_binned, fluxd_binned, drawstyle='steps-mid',
    ...          color='#ff7f0e', label='R~25')
    >>> plt.setp(plt.gca(), xlim=wrange, xlabel='Wavelength (μm)',
    ...          ylabel='Flux density (W/(m2 μm)')
    >>> plt.legend()
    >>> plt.tight_layout()

.. plot::

    import astropy.units as u
    import numpy as np
    import matplotlib.pyplot as plt
    from sbpy.calib import Sun
    wrange = 0.3, 0.8  # wavelength range
    d = 1 + 1 / 25
    n = int(np.ceil(np.log(wrange[1] / wrange[0]) / np.log(d)))
    wave_binned = wrange[0] * d**np.arange(n) * u.um
    sun = Sun.from_default()
    fluxd_binned = sun.observe(wave_binned, unit='W / (m2 um)')
    plt.plot(sun.wave.to('um'), sun.fluxd.to('W/(m2 um)'), drawstyle='steps-mid', color='#1f77b4', label='Native resolution')
    plt.plot(wave_binned, fluxd_binned, drawstyle='steps-mid', color='#ff7f0e', label='R~25')
    plt.setp(plt.gca(), xlim=wrange, xlabel='Wavelength (μm)', ylabel='Flux density (W/(m2 μm)')
    plt.legend()
    plt.tight_layout()


Reference/API
-------------
.. automodapi:: sbpy.calib
   :no-heading:
   :inherited-members:
