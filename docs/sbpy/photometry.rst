Photometry Module (`sbpy.photometry`)
=====================================

Introduction
------------

`~sbpy.photometry` provides routines for photometric phase curve
modeling of small bodies, and for retrieving `sbpy`'s built-in filter
bandpasses.


Disk-integrated Phase Function Models
-------------------------------------

Disk-integrated phase function refers to the brightness or reflectance of
an asteroid, often measured in magnitude, with respect to phase angle.  The IAU
adopted a number of models to describe the phase function of asteroids,
including the two-parameter H, G system (Bowell et al. 1989), the
three-parameter H, G1, G2, and the two-parameter H, G12 system derived from
the three-parameter system (Muinonen et al. 2010).  The H, G12 system was
further revised by Penttilä et al. (2016).  `~sbpy.photometry` provides
classes for all four models `~sbpy.photometry.HG`, `~sbpy.photometry.HG1G2`,
`~sbpy.photometry.HG12`, and `~sbpy.photometry.HG12_Pen16`, respectively.

The photometric model class can be initialized with the default parameters,
by supplying the model parameters as either dimensionless numbers or
`~astropy.units.Quantity`:

  >>> import astropy.units as u
  >>> from sbpy.photometry import HG, HG1G2, HG12, HG12_Pen16
  >>> m = HG()
  >>> print(m)
  Model: HG
  Inputs: ('x',)
  Outputs: ('y',)
  Model set size: 1
  Parameters:
       H   G
      --- ---
      8.0 0.4

  >>> m = HG(H = 3.34 * u.mag, G = 0.12, radius = 460 * u.km, wfb = 'V')
  >>> print(m)
  Model: HG
  Inputs: ('x',)
  Outputs: ('y',)
  Model set size: 1
  Parameters:
       H    G
      mag
      ---- ----
      3.34 0.12

The calculations that involve conversion between magnitude and reflectance
requires valid object size and `wfb` (wavelength/frequency/band) parameter to
be set for the photometric model.  The corresponding solar flux at the `wfb`
of the photometric model object has to be available through `~sbpy.calib`
system, or set by `~sbpy.calib.solar_fluxd.set`, which works with context
management `with` syntax.

Calculate geometric albedo, Bond albedo, and phase integral:

.. doctest-requires:: scipy

  >>> import astropy.units as u
  >>> from sbpy.calib import solar_fluxd
  >>> solar_fluxd.set({'V': -26.77 * u.mag})
  <ScienceState solar_fluxd: {'V': <Quantity -26.77 mag>}>
  >>> print(m.geomalb)  # doctest: +FLOAT_CMP
  0.09557298727307795
  >>> print(m.bondalb)  # doctest: +FLOAT_CMP
  0.03482207291799989
  >>> print(m.phaseint)  # doctest: +FLOAT_CMP
  0.3643505755292945

Note that the current version of `astropy.modeling.Model` doesn't support
`astropy.units.MagUnit` instance as model parameters.  For now one has to use
the dimensionless magnitude `~astropy.units.mag` in the phase function
parameter, and manually set solar flux in order for the conversion between
magnitude and reflectance to work.

The model class can also be initialized by a subclass of ``sbpy``'s
`~sbpy.data.DataClass`, such as `~sbpy.data.Phys`, as long as it contains the
model parameters:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

  >>> from sbpy.data import Phys
  >>> phys = Phys.from_sbdb('Ceres')
  >>> m = HG.from_phys(phys)
  INFO: Model initialized for 1 Ceres (A801 AA). [sbpy.photometry.core]
  >>> print(m)                          # doctest: +SKIP
  Model: HG
  Inputs: ('x',)
  Outputs: ('y',)
  Model set size: 1
  Parameters:
       H    G 
      mag
      ---- ----
      3.31 0.12

Note that in this case, model set is not supported.  Only one model can be
initialized with the first set of valid parameters in the input
`~sbpy.data.DataClass`.

To fit a photometric model, one may follow the standard procedure defined in
`astropy.modeling` submodule, first initializing a model, then using one of
the fitter classes defined in `astropy.modeling.fitting`
submodule, such as `~astropy.modeling.fitting.SLSQPLSQFitter`.  Note that
the `~sbpy.photometry.HG1G2` model requires that the G1 and G2 parameters
are bounded within [0, 1], as well as an inequality constraint,
0 <= 1 - G1 - G2 <= 1.  These constraints are implemented in sbpy via the
`bounds` parameter of `~astropy.modeling.Parameter` and the `ineqcons`
parameter of `~astropy.modeling.Model`.  Some fitters, such as
`astropy.modeling.LevMarLSQFitter`, do not support constrained fit via
the `ineqcons` parameter, though.

.. doctest-requires:: scipy

  >>> import numpy as np
  >>> import astropy.units as u
  >>> from astropy.modeling.fitting import SLSQPLSQFitter
  >>> # generate data to be fitted
  >>> model1 = HG(3.34 * u .mag, 0.12)
  >>> alpha = np.linspace(0, 40, 20) * u.deg
  >>> mag = model1(alpha) + (np.random.rand(20)*0.2 - 0.1) * u.mag
  >>> # fit new model
  >>> fitter = SLSQPLSQFitter()
  >>> model2 = HG()
  >>> model2 = fitter(model2, alpha, mag, verblevel=0)

Alternatively, one may use the class method
`~sbpy.photometry.DiskIntegratedPhaseFunc.from_obs` to
initialize a model directly from an `~sbpy.data.Obs` object by fitting the
data contained therein.

.. doctest-requires:: scipy

  >>> # use class method .from_obs
  >>> from astropy.modeling.fitting import SLSQPLSQFitter
  >>> fitter = SLSQPLSQFitter()
  >>> from sbpy.data import Obs
  >>> obs = Obs.from_dict({'alpha': alpha, 'mag': mag})
  >>> model3 = HG12.from_obs(obs, fitter, 'mag', verblevel=0)

One can also initialize a model set from multiple columns in the input
`~sbpy.data.Obs` object if it contains more than one columns of brightness
measurements.  The columns to be fitted are specified by a keyward argument
``fields``.  By default, the column ``'mag'`` will be fitted.

.. doctest-requires:: scipy

  >>> # Initialize model set
  >>> model4 = HG(5.2 * u.mag, 0.18)
  >>> mag4 = model4(alpha) + (np.random.rand(20)*0.2 - 0.1) * u.mag
  >>> fitter = SLSQPLSQFitter()
  >>> obs = Obs.from_dict({'alpha': alpha, 'mag': mag, 'mag1': mag4})
  >>> model5 = HG.from_obs(obs, fitter, fields=['mag', 'mag1'], verblevel=0)

.. _filter-bandpasses:

Filter Bandpasses
-----------------
A few filter bandpasses are included with `sbpy` for internal tests and your convenience.  The function `~sbpy.photometry.bandpass` will return the filter transmission as a `~synphot.spectrum.SpectralElement` (requires `synphot`):

.. doctest-requires:: synphot

  >>> from sbpy.photometry import bandpass
  >>> bp = bandpass('Cousins R')
  >>> print(bp.avgwave())    # doctest: +FLOAT_CMP
  6499.914781904409 Angstrom

For other bandpasses, obtain the photon-counting relative spectral response curves as a two-column file.  If the first column is wavelength in Angstroms, and the second is the response, then read the file with:

.. doctest-requires:: synphot

  >>> from synphot import SpectralElement             # doctest: +SKIP
  >>> bp = SpectralElement.from_file('filename.txt')  # doctest: +SKIP

See `synphot.spectrum.SpectralElement` for other options and file format details.

Reference/API
-------------
.. automodapi:: sbpy.photometry
    :no-heading:
