Photometry Module (`sbpy.photometry`)
=====================================

Introduction
------------

`~sbpy.photometry` provides routines for photometric phase curve
modeling of small bodies.


Disk-integrated Phase Function Models
-------------------------------------

Disk-integrated phase function refers to the brightness or reflectance of
an asteroid, often measured in magnitude, with respect to phase angle.  The IAU
adopted a number of models to describe the phase function of asteroids,
including the two-parameter H, G system (Bowell et al. 1989), the
three-parameter H, G1, G2, and the two-parameter H, G12 system derived from
the three-parameter system (Muinonen et al. 2010).  The H, G12 system was
further revised by Penttilä et al. (2016).  `~sbpy.photometry` provides all
four models classes `~sbpy.photometry.HG`, `~sbpy.photometry.HG1G2`,
`~sbpy.photometry.HG12`, and `~sbpy.photometry.HG12_Pen16`, respectively.

The photometric model class can be initialized with the default parameters,
by supplying the model parameters as either dimensionless numbers or the
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
  >>> m = HG(H=3.34*u.mag, G=0.12*u.dimensionless_unscaled, radius=460*u.km)
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

Calculate geometric albedo, Bond albedo, and phase integral:

  >>> print(m.geoalb)  # doctest: +FLOAT_CMP
  0.09825058857735823
  >>> print(m.bondalb)  # doctest: +FLOAT_CMP
  0.035797658494252343
  >>> print(m.phaseint)  # doctest: +FLOAT_CMP
  0.36435057552929395

Users can supply a solar magnitude corresponding to the magnitude system of
the H-parameter.  The default is the apparent V-magnitude of the Sun,
M_sun = -26.74 mag.

The model class can also be initialized by a subclass of ``sbpy``'s
`~sbpy.data.DataClass`, such as `~sbpy.data.Phys`, as long as it contains the
model parameters:

  >>> from sbpy.data import Phys
  >>> phys = Phys.from_sbdb('Ceres')    # doctest: +REMOTE_DATA
  >>> m = HG.from_phys(phys)            # doctest: +REMOTE_DATA
  INFO: Model initialized for 1 Ceres. [sbpy.photometry.core]
  >>> print(m)                          # doctest: +REMOTE_DATA
  Model: HG
  Inputs: ('x',)
  Outputs: ('y',)
  Model set size: 1
  Parameters:
       H    G
      ---- ----
      3.34 0.12

Note that  model set is not supported.  Only one model can be initialized with
the first set of valid parameters in the input `~sbpy.data.DataClass`.

To fit a photometric model, one may follow the standard procedure defined in
``astropy.modeling`` submodule, first initializing a model, then using one of
the fitter classes defined in ``astropy``'s `~astropy.modeling.fitting`
submodule, such as `~astropy.modeling.fitting.LevMarLSQFitter`.

  >>> import numpy as np
  >>> from astropy.modeling.fitting import LevMarLSQFitter
  >>> # generate data to be fitted
  >>> model1 = HG(3.34, 0.12)
  >>> alpha = np.linspace(0, 40, 20) * u.deg
  >>> mag = model1(alpha) + np.random.rand(20)*0.2-0.1
  >>> # fit new model
  >>> fitter = LevMarLSQFitter()
  >>> model2 = HG()
  >>> model2 = fitter(model2, alpha, mag)
  >>> print(model2)  # doctest: +SKIP
  Model: HG
  Inputs: ('x',)
  Outputs: ('y',)
  Model set size: 1
  Parameters:
              H                  G

      ----------------- -------------------
      3.305001580933622 0.08532754955207918

Alternatively, one may use the wrapper method
`~sbpy.photometry.DiskIntegratedPhaseFunc.fit` pre-defined in the photometric
model classes to fit magnitude data and return a model class object, or use
class method `~sbpy.photometry.DiskIntegratedPhaseFunc.from_data` to
initialize a model directly from data by fitting.

  >>> # use .fit method
  >>> model3 = HG1G2().fit(alpha, mag)
  >>> print(model3)  # doctest: +SKIP
  Model: HG1G2
  Inputs: ('x',)
  Outputs: ('y',)
  Model set size: 1
  Parameters:
              H                  G1                  G2

      ------------------ ------------------ -------------------
      3.3391210178858013 0.6551118013498317 0.09728839079940735

  >>> # use class method .from_data
  >>> model4 = HG12.from_data(alpha, mag)
  >>> print(model4)  # doctest: +SKIP
  Model: HG12
  Inputs: ('x',)
  Outputs: ('y',)
  Model set size: 1
  Parameters:
              H                G12

      ----------------- ------------------
      3.424576008941485 0.8052670564595159


Reference/API
-------------
.. automodapi:: sbpy.photometry
    :no-heading:
