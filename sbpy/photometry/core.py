# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Photometry Module

created on June 23, 2017

"""

__all__ = ['DiskIntegratedPhaseFunc', 'LinearPhaseFunc', 'HG', 'HG12BaseClass',
           'HG12', 'HG1G2', 'HG12_Pen16', 'NonmonotonicPhaseFunctionWarning']

from collections import OrderedDict
import warnings
import numpy as np
from numbers import Number
from scipy.integrate import quad
from astropy.modeling import (Fittable1DModel, Parameter)
from astropy.table import Column
import astropy.units as u
from astropy import log
from ..data import (DataClass, Phys, Obs, Ephem, dataclass_input,
                    quantity_to_dataclass)
from ..bib import cite
from ..units import reflectance
from ..exceptions import SbpyWarning


class _spline(object):

    """Cubic spline class

    Spline function is defined by function values at nodes and the first
    derivatives at both ends.  Outside the range of nodes, the extrapolations
    are linear based on the first derivatives at the corresponding ends.
    """

    def __init__(self, x, y, dy):
        """
        Spline initialization

        Parameters
        ----------
        x, y : array_like float
            The (x, y) values at nodes that defines the spline
        dy : array_like float with two elements
            The first derivatives of the left and right ends of the nodes
        """
        from numpy.linalg import solve
        from numpy.polynomial.polynomial import Polynomial
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.dy = np.asarray(dy)
        n = len(self.y)
        h = self.x[1:]-self.x[:-1]
        r = (self.y[1:]-self.y[:-1])/(self.x[1:]-self.x[:-1])
        B = np.zeros((n-2, n))
        for i in range(n-2):
            k = i+1
            B[i, i:i+3] = [h[k], 2*(h[k-1]+h[k]), h[k-1]]
        C = np.empty((n-2, 1))
        for i in range(n-2):
            k = i+1
            C[i] = 3*(r[k-1]*h[k]+r[k]*h[k-1])
        C[0] = C[0]-self.dy[0]*B[0, 0]
        C[-1] = C[-1]-self.dy[1]*B[-1, -1]
        B = B[:, 1:n-1]
        dys = solve(B, C)
        dys = np.array(
            [self.dy[0]] + [tmp for tmp in dys.flatten()] + [self.dy[1]])
        A0 = self.y[:-1]
        A1 = dys[:-1]
        A2 = (3*r-2*dys[:-1]-dys[1:])/h
        A3 = (-2*r+dys[:-1]+dys[1:])/h**2
        self.coef = np.array([A0, A1, A2, A3]).T
        self.polys = [Polynomial(c) for c in self.coef]
        self.polys.insert(0, Polynomial(
            [self.y[0]-self.x[0]*self.dy[0], self.dy[0]]))
        self.polys.append(Polynomial(
            [self.y[-1]-self.x[-1]*self.dy[-1], self.dy[-1]]))

    def __call__(self, x):
        x = np.asarray(x)
        out = np.zeros_like(x)
        idx = x < self.x[0]
        if idx.any():
            out[idx] = self.polys[0](x[idx])
        for i in range(len(self.x)-1):
            idx = (self.x[i] <= x) & (x < self.x[i+1])
            if idx.any():
                out[idx] = self.polys[i+1](x[idx]-self.x[i])
        idx = (x >= self.x[-1])
        if idx.any():
            out[idx] = self.polys[-1](x[idx])
        return out


class NonmonotonicPhaseFunctionWarning(SbpyWarning):
    pass


class DiskIntegratedPhaseFunc(Fittable1DModel):
    """Base class for disk-integrated phase function model

    Examples
    --------
    Define a linear phase function with phase slope 0.04 mag/deg, and
    study its properties:

    >>> # Define a disk-integrated phase function model
    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.modeling import Parameter
    >>> from sbpy.calib import solar_fluxd
    >>> from sbpy.photometry import DiskIntegratedPhaseFunc
    >>>
    >>> class LinearPhaseFunc(DiskIntegratedPhaseFunc):
    ...
    ...     _unit = 'mag'
    ...     H = Parameter()
    ...     S = Parameter()
    ...
    ...     @staticmethod
    ...     def evaluate(a, H, S):
    ...         return H + S * a
    ...
    >>> linear_phasefunc = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg,
    ...     radius = 300 * u.km, wfb = 'V')
    >>> pha = np.linspace(0, 180, 200) * u.deg
    >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
    ...     mag = linear_phasefunc.to_mag(pha)
    ...     ref = linear_phasefunc.to_ref(pha)
    ...     geomalb = linear_phasefunc.geomalb
    ...     phaseint = linear_phasefunc.phaseint
    ...     bondalb = linear_phasefunc.bondalb
    >>> print('Geometric albedo is {0:.3}'.format(geomalb))
    Geometric albedo is 0.0487
    >>> print('Bond albedo is {0:.3}'.format(bondalb))
    Bond albedo is 0.0179
    >>> print('Phase integral is {0:.3}'.format(phaseint))
    Phase integral is 0.367

    Initialization with subclass of `~sbpy.data.DataClass`:

    The subclassed models can either be initialized by model parameters, or by
    subclass of `~sbpy.data.DataClass`.  Below example uses the `HG` model
    class.

    >>> from sbpy.photometry import HG
    >>> from sbpy.data import Phys, Orbit, Ephem
    >>>
    >>> # Initialize from physical parameters pulled from JPL SBDB
    >>> phys = Phys.from_sbdb('Ceres')       # doctest: +REMOTE_DATA
    >>> print(phys['targetname','H','G'])    # doctest: +REMOTE_DATA
    <QTable length=1>
    targetname    H       G
       str7    float64 float64
    ---------- ------- -------
       1 Ceres    3.34    0.12
    >>> m = HG.from_phys(phys)                  # doctest: +REMOTE_DATA
    INFO: Model initialized for 1 Ceres. [sbpy.photometry.core]
    >>> print(m)                             # doctest: +REMOTE_DATA
    Model: HG
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Parameters:
         H    G
        ---- ----
        3.34 0.12
    >>> print(m.meta['targetname'])          # doctest: +REMOTE_DATA
    1 Ceres
    >>> print(m.radius)                      # doctest: +REMOTE_DATA
    469.7 km
    >>>
    >>> # Initialize from orbital elements pulled from JPL Horizons that also
    >>> # contain the H and G parameters
    >>> elem = Orbit.from_horizons('Ceres')  # doctest: +REMOTE_DATA
    >>> print(elem['targetname','H','G'])    # doctest: +REMOTE_DATA
    <QTable masked=True length=1>
    targetname    H       G
                 mag
       str7    float64 float64
    ---------- ------- -------
       1 Ceres    3.34    0.12
    >>> m = HG.from_phys(elem)                    # doctest: +REMOTE_DATA
    INFO: Model initialized for 1 Ceres. [sbpy.photometry.core]
    >>>
    >>> # Failed initialization due to the lack of field 'G'
    >>> phys = Phys.from_sbdb('12893')       # doctest: +REMOTE_DATA
    >>> print('G' in phys.field_names)      # doctest: +REMOTE_DATA
    False
    >>> m = HG(data=phys)                    # doctest: +SKIP
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    KeyError: 'field G not available.'
    """

    # Some phase function models are defined in magnitude space, such as the
    # IAU H, G system.  Some phase function models are defined in reflectance
    # space, such as the disk-integrated phase function of the Hapke model.
    # _unit defines which unit the model is defined in.
    _unit = None

    # The default unit for model input when the model is dimensional
    input_units = {'x': u.rad}

    # Whether or not the model input is allowed to be dimensionless
    input_units_allow_dimensionless = {'x': True}

    @u.quantity_input(radius=u.km)
    def __init__(self, *args, radius=None, wfb=None, **kwargs):
        """Initialize DiskIntegratedPhaseFunc

        Parameters
        ----------
        radius : astropy.units.Quantity, optional
            Radius of object.  Required if conversion between magnitude and
            reflectance is involved.
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, string
            Wavelengths, frequencies, or bandpasses.  Bandpasses may
            be a filter name (string).  Required if conversion between
            magnitude and reflectance is involved.
        **kwargs : optional parameters accepted by
            `astropy.modeling.Model.__init__()`
        """
        super().__init__(*args, **kwargs)
        self.radius = radius
        self.wfb = wfb

    def _check_unit(self):
        if self._unit is None:
            raise ValueError('the unit of phase function is unknown')

    @property
    def geomalb(self):
        """Geometric albedo"""
        alb = np.pi*self.to_ref(0.*u.rad)
        if hasattr(alb, 'unit') and (alb.unit == 1/u.sr):
            alb = alb*u.sr
        return alb

    @property
    def bondalb(self):
        """Bond albedo"""
        return self.geomalb*self.phaseint

    @property
    def phaseint(self):
        """Phase integral"""
        return self._phase_integral()

    @classmethod
    def from_phys(cls, phys, **kwargs):
        """Initialize an object from `~sbpy.data.Phys` object

        Parameters
        ----------
        phys : `~sbpy.data.Phys`
            Contains the parameters needed to initialize the model class
            object.  If the required field is not found, then an `KeyError`
            exception will be thrown.
        **kwargs : optional parameters accepted by
            `astropy.modeling.Model.__init__()`

        Returns
        ------
        Object of `DiskIntegratedPhaseFunc` subclass
            The phase function model object

        Examples
        --------
        Initialization with `~sbpy.data.Phys`.  This example uses the `HG`
        model class.

        >>> from sbpy.photometry import HG
        >>> from sbpy.data import Phys
        >>>
        >>> # Initialize from physical parameters pulled from JPL SBDB
        >>> phys = Phys.from_sbdb('Ceres')      # doctest: +REMOTE_DATA
        >>> print(phys['targetname','H','G'])   # doctest: +REMOTE_DATA
        <QTable length=1>
        targetname    H       G
           str7    float64 float64
        ---------- ------- -------
           1 Ceres    3.34    0.12
        >>> m = HG.from_phys(phys)              # doctest: +REMOTE_DATA
        INFO: Model initialized for 1 Ceres. [sbpy.photometry.core]
        >>> print(m)                            # doctest: +REMOTE_DATA
        Model: HG
        Inputs: ('x',)
        Outputs: ('y',)
        Model set size: 1
        Parameters:
             H    G
            ---- ----
            3.34 0.12
        >>> print(m.meta['targetname'])         # doctest: +REMOTE_DATA
        1 Ceres
        >>> print(m.radius)                     # doctest: +REMOTE_DATA
        469.7 km
        >>>
        >>> # Failed initialization due to the lack of field 'G'
        >>> phys = Phys.from_sbdb('12893')      # doctest: +REMOTE_DATA
        >>> print('G' in phys.field_names)      # doctest: +REMOTE_DATA
        False
        >>> m = HG.from_phys(phys)              # doctest: +REMOTE_DATA +SKIP
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        KeyError: 'field G not available.'
        """

        par = {}
        valid = np.ones(len(phys), dtype=bool)
        for p in cls.param_names:
            par[p] = phys[p]
            valid = valid & np.isfinite(par[p])
        if valid.any():
            valid = list(valid).index(True)
            for p in cls.param_names:
                par[p] = par[p][valid]
            meta = kwargs.pop('meta', OrderedDict())
            if 'targetname' in phys.field_names:
                meta.update({'targetname': phys['targetname'][valid]})
            kwargs['meta'] = meta
            for p in cls.param_names:
                val = kwargs.pop(p, None)
            try:
                par['radius'] = phys['diameter'][valid]/2
            except KeyError:
                pass
            if 'targetname' in meta.keys():
                log.info("Model initialized for {}.".format(
                    meta['targetname']))
            else:
                log.info("Model initialized.")
            kwargs.update(par)
        else:
            raise ValueError(
                'no valid model parameters found in `data` keyword')
        out = cls(**kwargs)
        return out

    def to_phys(self):
        """Wrap the model into a `sbpy.data.Phys` object

        Returns
        -------
        `~sbpy.data.Phys` object

        Examples
        --------
        >>> import astropy.units as u
        >>> from sbpy.calib import solar_fluxd
        >>> from sbpy.photometry import HG
        >>> from sbpy.data import Phys
        >>>
        >>> # Initialize from physical parameters pulled from JPL SBDB
        >>> phys = Phys.from_sbdb('Ceres')             # doctest: +REMOTE_DATA
        >>> print(phys['targetname','radius','H','G']) # doctest: +REMOTE_DATA
        <QTable length=1>
        targetname  radius    H       G
                      km
           str7    float64 float64 float64
        ---------- ------- ------- -------
           1 Ceres   469.7    3.34    0.12
        >>> m = HG.from_phys(phys)   # doctest: +REMOTE_DATA
        INFO: Model initialized for 1 Ceres. [sbpy.photometry.core]
        >>> m.wfb = 'V'              # doctest: +REMOTE_DATA
        >>> m.H = m.H * u.mag        # doctest: +REMOTE_DATA
        >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
        ...     p = m.to_phys()      # doctest: +REMOTE_DATA
        >>> print(type(p))           # doctest: +REMOTE_DATA
        <class 'sbpy.data.phys.Phys'>
        >>> print(p)                 # doctest: +REMOTE_DATA
        <QTable length=1>
        targetname diameter    H       G             pv                  A
                      km      mag
           str7    float64  float64 float64       float64             float64
        ---------- -------- ------- ------- ------------------- -------------------
           1 Ceres    939.4    3.34    0.12 0.09166630037900923 0.03339866929973315
        """
        cols = {}
        if (self.meta is not None) and ('targetname' in self.meta.keys()):
            val = self.meta['targetname']
            if isinstance(val, str):
                val = [val]
            cols['targetname'] = val
        if self.radius is not None:
            cols['diameter'] = self.radius * 2
        for p in self.param_names:
            val = getattr(self, p)
            if val.quantity is None:
                cols[p] = val.value
            else:
                cols[p] = val.quantity
        try:
            cols['pv'] = self.geomalb
            cols['A'] = self.bondalb
        except ValueError:
            pass
        return Phys.from_dict(cols)

    @classmethod
    @dataclass_input(obs=Obs)
    def from_obs(cls, obs, fitter, fields='mag', init=None, **kwargs):
        """Instantiate a photometric model class object from data

        Parameters
        ----------
        obs : `~sbpy.data.DataClass`, dict_like
            If `~sbpy.data.DataClass` or dict_like, must contain
            ``'phaseangle'`` or the equivalent names (see
            `~sbpy.data.DataClass`).  If any distance (heliocentric and
            geocentric) is provided, then they will be used to correct
            magnitude to 1 au before fitting.
        fitter : `~astropy.modeling.fitting.Fitter`
            The fitter to be used for fitting.
        fields : str or array_like of str
            The field name or names in ``obs`` to be fitted.  If an array_like
            str, then multiple fields will be fitted one by one and a model
            set will be returned.  In this case, ``.meta['fields']`` of the
            returned object contains the names of fields fitted.
        init : numpy array, `~astropy.units.Quantity`, optional
            The initial parameters for model fitting.  Its first dimension has
            the length of the model parameters, and its second dimension has
            the length of ``n_model`` if multiple models are fitted.
        **kwargs : optional parameters accepted by `fitter()`.
            Note that the magnitude uncertainty can also be supplied to the fit
            via `weights` keyword for all fitters provided by
            `~astropy.modeling.fitting`.

        Returns
        -------
        Object of `DiskIntegratedPhaseFunc` subclass
            The best-fit model class object.

        Examples
        --------
        >>> from sbpy.photometry import HG # doctest: +SKIP
        >>> from sbpy.data import Misc # doctest: +SKIP
        >>> from astropy.modeling.fitting import LevMarLSQFitter
        >>> fitter = LevMarLSQFitter()
        >>> obs = Misc.mpc_observations('Bennu') # doctest: +SKIP
        >>> hg = HG() # doctest: +SKIP
        >>> best_hg = hg.from_obs(obs, eph['mag'], fitter) # doctest: +SKIP
        """
        pha = obs['alpha']

        if isinstance(fields, (str, bytes)):
            n_models = 1
        else:
            n_models = len(fields)
        if init is not None:
            init = np.asanyarray(init)

        dist_corr = cls()._distance_module(obs)
        if n_models == 1:
            mag = obs[fields]
            if isinstance(mag, u.Quantity):
                dist_corr = u.Quantity(dist_corr).to(u.mag, u.logarithmic())
            else:
                dist_corr = -2.5 * alog10(dist_corr)
            mag0 = mag + dist_corr
            if init is None:
                m0 = cls()
            else:
                m0 = cls(*init)
            return fitter(m0, pha, mag0, **kwargs)
        else:
            if init is not None:
                sz1 = init.shape
                sz2 = len(cls.param_names), n_models
                if sz1 != sz2:
                    raise ValueError('`init` must have a shape of ({}, {}),'
                                     ' shape {} is given.'.format(sz2[0],
                                     sz2[1], sz1))
            par = np.zeros((len(cls.param_names), n_models))
            for i in range(n_models):
                mag = obs[fields[i]]
                if isinstance(mag, u.Quantity):
                    dist_corr1 = u.Quantity(dist_corr).to(u.mag,
                            u.logarithmic())
                else:
                    dist_corr1 = -2.5 * alog10(dist_corr)
                mag0 = mag + dist_corr1
                if init is None:
                    m0 = cls()
                else:
                    m0 = cls(*init[:, i])
                m = fitter(m0, pha, mag0, **kwargs)
                par[:, i] = m.parameters
            pars_list = []
            for i, p_name in enumerate(cls.param_names):
                p = getattr(m, p_name)
                if p.unit is None:
                    pars_list.append(par[i])
                else:
                    pars_list.append(par[i]*p.unit)
            model = cls(*pars_list, n_models=n_models)
            if not isinstance(model.meta, dict):
                model.meta = OrderedDict()
            model.meta['fields'] = fields
            return model

    @dataclass_input(eph=Ephem)
    def _distance_module(self, eph):
        """Return the correction magnitude or factor for heliocentric distance
        and observer distance

        Parameters
        ----------
        eph : any type
            If `~sbpy.data.Ephem` or dict_like, then the relevant fields, such
            as 'rh' and 'delta' or the equivalent will be searched and, if
            exist, used to calculate distance correction.  If non-exist, then
            no correction will be included for the corresponding field.  If no
            unit is provided via type `~astropy.units.Quantity`, then the
            distance is assumed to be in unit of au.  For any other data type,
            a factor 1 or magnitude of 0 will be returned (implying no
            correction).

        Returns
        -------
        float or numpy array
            Linear factors to be applied to flux to correct to heliocentric
            distance and observer distance of both 1 au.
        """
        module = 1.
        try:
            rh = eph['r']
            if isinstance(rh, u.Quantity):
                rh = rh.to('au').value
            module = module * rh * rh
        except (KeyError, TypeError):
            pass
        try:
            delta = eph['delta']
            if isinstance(delta, u.Quantity):
                delta = delta.to('au').value
            module = module * delta * delta
        except (KeyError, TypeError):
            pass
        return np.asarray(module)

    @quantity_to_dataclass(eph=(Ephem, 'alpha'))
    def to_mag(self, eph, unit=None, append_results=False, **kwargs):
        """Calculate phase function in magnitude

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`, numbers, iterables of numbers, or
            `~astropy.units.Quantity`
            If `~sbpy.data.Ephem` or dict_like, ephemerides of the object that
            can include phase angle, heliocentric and geocentric distances via
            keywords `phase`, `r` and `delta`.  If float or array_like, then
            the phase angle of object.  If any distance (heliocentric and
            geocentric) is not provided, then it will be assumed to be 1 au.
            If no unit is provided via type `~astropy.units.Quantity`, then
            radians is assumed for phase angle, and au is assumed for
            distances.
        unit : `astropy.units.mag`, `astropy.units.MagUnit`, optional
            The unit of output magnitude.  The corresponding solar magnitude
            must be available either through `~sbpy.calib.sun` module or set
            by `~sbpy.calib.solar_fluxd.set`.
        append_results : bool, optional
            Controls the return of this method.
        **kwargs : optional parameters accepted by
            `astropy.modeling.Model.__call__`

        Returns
        -------
        `~astropy.units.Quantity`, array if ``append_results == False``
        `~sbpy.data.Ephem` if ``append_results == True``

        When ``apend_results == False``: The calculated magnitude will be
        returned.

        When ``append_results == True``:  If ``eph`` is a `~sbpy.data.Ephem`
        object, then the calculated magnitude will be appended to ``eph`` as
        a new column.  Otherwise a new `~sbpy.data.Ephem` object is created to
        contain the input ``eph`` and the calculated magnitude in two columns.

        Examples
        --------
        >>> import numpy as np
        >>> from astropy import units as u
        >>> from sbpy.photometry import HG
        >>> from sbpy.data import Ephem
        >>> ceres_hg = HG(3.34 * u.mag, 0.12)
        >>> # parameter `eph` as `~sbpy.data.Ephem` type
        >>> eph = Ephem.from_dict({'alpha': np.linspace(0,np.pi*0.9,200)*u.rad,
        ...              'r': np.repeat(2.7*u.au, 200),
        ...              'delta': np.repeat(1.8*u.au, 200)})
        >>> mag1 = ceres_hg.to_mag(eph)
        >>> # parameter `eph` as numpy array
        >>> pha = np.linspace(0, 170, 200) * u.deg
        >>> mag2 = ceres_hg.to_mag(pha)
        """
        self._check_unit()
        pha = eph['alpha']
        if len(pha) == 1:
            pha = pha[0]
        out = self(pha, **kwargs)
        if self._unit == 'ref':
            if unit is None:
                raise ValueError('Magnitude unit is not specified.')
            if self.radius is None:
                raise ValueError(
                    'Cannot calculate phase funciton in magnitude because the'
                    ' size of object is unknown.')
            if self.wfb is None:
                raise ValueError('Wavelength/Frequency/Band is unknown.')
            out = out.to(unit, reflectance(self.wfb,
                    cross_section=np.pi * self.radius**2))
        dist_corr = self._distance_module(eph)
        dist_corr = u.Quantity(dist_corr).to(u.mag, u.logarithmic())
        out = out - dist_corr
        if append_results:
            name = 'mag'
            i = 1
            while name in eph.field_names:
                name = 'mag'+str(i)
                i += 1
            eph.table.add_column(Column(out, name=name))
            return eph
        else:
            return out

    @quantity_to_dataclass(eph=(Ephem, 'alpha'))
    def to_ref(self, eph, normalized=None, append_results=False, **kwargs):
        """Calculate phase function in average bidirectional reflectance

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`, numbers, iterables of numbers, or
            `~astropy.units.Quantity`
            If `~sbpy.data.Ephem` or dict_like, ephemerides of the object that
            can include phase angle, heliocentric and geocentric distances via
            keywords `phase`, `r` and `delta`.  If float or array_like, then
            the phase angle of object.  If any distance (heliocentric and
            geocentric) is not provided, then it will be assumed to be 1 au.
            If no unit is provided via type `~astropy.units.Quantity`, then
            radians is assumed for phase angle, and au is assumed for
            distances.
        normalized : number, `~astropy.units.Quantity`
            The angle to which the reflectance is normalized.
        append_results : bool
            Controls the return of this method.
        **kwargs : optional parameters accepted by
            `astropy.modeling.Model.__call__`

        Returns
        -------
        `~astropy.units.Quantity`, array if ``append_results == False``
        `~sbpy.data.Ephem` if ``append_results == True``

        When ``apend_results == False``: The calculated reflectance will be
        returned.

        When ``append_results == True``:  If ``eph`` is a `~sbpy.data.Ephem`
        object, then the calculated reflectance will be appended to ``eph`` as
        a new column.  Otherwise a new `~sbpy.data.Ephem` object is created to
        contain the input ``eph`` and the calculated reflectance in two
        columns.

        Examples
        --------
        >>> import numpy as np
        >>> from astropy import units as u
        >>> from sbpy.calib import solar_fluxd
        >>> from sbpy.photometry import HG
        >>> from sbpy.data import Ephem
        >>> ceres_hg = HG(3.34 * u.mag, 0.12, radius = 480 * u.km, wfb= 'V')
        >>> # parameter `eph` as `~sbpy.data.Ephem` type
        >>> eph = Ephem.from_dict({'alpha': np.linspace(0,np.pi*0.9,200)*u.rad,
        ...              'r': np.repeat(2.7*u.au, 200),
        ...              'delta': np.repeat(1.8*u.au, 200)})
        >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
        ...     ref1 = ceres_hg.to_ref(eph)
        ...     # parameter `eph` as numpy array
        ...     pha = np.linspace(0, 170, 200) * u.deg
        ...     ref2 = ceres_hg.to_ref(pha)
        """
        self._check_unit()
        pha = eph['alpha']
        if len(pha) == 1:
            pha = pha[0]
        out = self(pha, **kwargs)
        if normalized is not None:
            norm = self(normalized, **kwargs)
        if self._unit == 'ref':
            if normalized is not None:
                out /= norm
        else:
            if normalized is None:
                if self.radius is None:
                    raise ValueError(
                        'Cannot calculate phase function in reflectance unit'
                        ' because the size of object is unknown.  Normalized'
                        ' phase function can be calculated.')
                if self.wfb is None:
                    raise ValueError('Wavelength/Frequency/Band is unknown.')
                out = out.to('1/sr', reflectance(self.wfb,
                        cross_section=np.pi*self.radius**2))
            else:
                out = out - norm
                out = out.to('', u.logarithmic())
        if append_results:
            name = 'ref'
            i = 1
            while name in eph.field_names:
                name = 'ref'+str(i)
                i += 1
            eph.table.add_column(Column(out, name=name))
            return eph
        else:
            return out

    def _phase_integral(self, integrator=quad):
        """Calculate phase integral with numerical integration

        Parameters
        ----------
        integrator : function, optinonal
            Numerical integrator, default is `~scipy.integrate.quad`.
            If caller supplies a numerical integrator, it must has the same
            return signature as `~scipy.integrator.quad`, i.e., a tuple of
            (y, ...), where `y` is the result of numerical integration

        Returns
        -------
        Float, phase integral

        Examples
        --------
        >>> import astropy.units as u
        >>> from sbpy.calib import solar_fluxd
        >>> from sbpy.photometry import HG
        >>> ceres_hg = HG(3.34 * u.mag, 0.12, radius = 480 * u.km, wfb = 'V')
        >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
        ...     print('{0:.3}'.format(ceres_hg._phase_integral()))
        0.364
        """
        def integrand(x):
            return 2*self.to_ref(x * u.rad, normalized=0. * u.rad) * \
                np.sin(x * u.rad)
        return integrator(integrand, 0, np.pi)[0]


class LinearPhaseFunc(DiskIntegratedPhaseFunc):
    """Linear phase function model

    Examples
    --------
    >>> # Define a linear phase function model with absolute magnitude
    >>> # H = 5 and slope = 0.04 mag/deg = 2.29 mag/rad
    >>> import astropy.units as u
    >>> from sbpy.calib import solar_fluxd
    >>> from sbpy.photometry import LinearPhaseFunc
    >>>
    >>> linear_phasefunc = LinearPhaseFunc(5 * u.mag, 0.04 * u.mag/u.deg,
    ...     radius = 300 * u.km, wfb = 'V')
    >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
    ...     pha = np.linspace(0, 180, 200) * u.deg
    ...     mag = linear_phasefunc.to_mag(pha)
    ...     ref = linear_phasefunc.to_ref(pha)
    ...     geomalb = linear_phasefunc.geomalb
    ...     phaseint = linear_phasefunc.phaseint
    ...     bondalb = linear_phasefunc.bondalb
    >>> print('Geometric albedo is {0:.3}'.format(geomalb))
    Geometric albedo is 0.0487
    >>> print('Bond albedo is {0:.3}'.format(bondalb))
    Bond albedo is 0.0179
    >>> print('Phase integral is {0:.3}'.format(phaseint))
    Phase integral is 0.367

    """

    _unit = 'mag'
    H = Parameter(description='Absolute magnitude')
    S = Parameter(description='Linear slope (mag/deg)')
    input_units = {'x': u.deg}

    @staticmethod
    def evaluate(a, H, S):
        return H + S * a

    @staticmethod
    def fit_deriv(a, H, S):
        if hasattr(a, '__iter__'):
            ddh = np.ones_like(a)
        else:
            ddh = 1.
        dds = a
        return [ddh, dds]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict([('H', outputs_unit['y']),
                            ('S', outputs_unit['y']/inputs_unit['x'])])


class HG(DiskIntegratedPhaseFunc):
    """HG photometric phase model (Bowell et al. 1989)

    Examples
    --------

    >>> # Define the phase function for Ceres with H = 3.34, G = 0.12
    >>> import astropy.units as u
    >>> from sbpy.calib import solar_fluxd
    >>> from sbpy.photometry import HG
    >>> ceres = HG(3.34 * u.mag, 0.12, radius = 480 * u.km, wfb = 'V')
    >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
    ...     print('geometric albedo = {0:.4f}'.format(ceres.geomalb))
    ...     print('phase integral = {0:.4f}'.format(ceres.phaseint))
    geometric albedo = 0.0878
    phase integral = 0.3644

    """

    _unit = 'mag'
    H = Parameter(description='H parameter', default=8)
    G = Parameter(description='G parameter', default=0.4)

    @cite({'definition': '1989aste.conf..524B'})
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @G.validator
    def G(self, value):
        """Validate parameter G to avoid non-monotonic phase function

        If G > 1.194, the phase function could potentially be non-monotoic,
        and a warning will be issued.
        """
        if np.any(value > 1.194):
            warnings.warn(
                'G parameter could result in a non-monotonic phase function',
                NonmonotonicPhaseFunctionWarning)

    @staticmethod
    def _hgphi(pha, i):
        """Core function in IAU HG phase function model

        Parameters
        ----------
        pha : float or array_like of float
            Phase angle
        i   : int in [1, 2]
            Choose the form of function

        Returns
        -------
        numpy array of float

        Note
        ----
        See Bowell et al. (1989), Eq. A4.
        """

        if i not in [1, 2]:
            raise ValueError('i needs to be 1 or 2, {0} received'.format(i))

        a, b, c = [3.332, 1.862], [0.631, 1.218], [0.986, 0.238]
        pha_half = pha*0.5
        sin_pha = np.sin(pha)
        tan_pha_half = np.tan(pha_half)
        w = np.exp(-90.56 * tan_pha_half * tan_pha_half)
        phiis = 1 - c[i-1]*sin_pha/(0.119+1.341*sin_pha -
                                    0.754*sin_pha*sin_pha)
        phiil = np.exp(-a[i-1] * tan_pha_half**b[i-1])
        return w*phiis + (1-w)*phiil

    @staticmethod
    def evaluate(pha, hh, gg):
        func = (1-gg)*HG._hgphi(pha, 1)+gg*HG._hgphi(pha, 2)
        if isinstance(func, u.Quantity):
            func = func.value
        func = -2.5 * np.log10(func)
        if isinstance(hh, u.Quantity):
            func = func * hh.unit
        return hh + func

    @staticmethod
    def fit_deriv(pha, hh, gg):
        if hasattr(pha, '__iter__'):
            ddh = np.ones_like(pha)
        else:
            ddh = 1.
        phi1 = HG._hgphi(pha, 1)
        phi2 = HG._hgphi(pha, 2)
        ddg = 1.085736205*(phi1-phi2)/((1-gg)*phi1+gg*phi2)
        return [ddh, ddg]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict([('H', outputs_unit['y']),
                            ('G', u.dimensionless_unscaled)])


class HG12BaseClass(DiskIntegratedPhaseFunc):
    """Base class for IAU HG1G2 model and HG12 model"""

    _unit = 'mag'

    @cite({'definition': '2010Icar..209..542M'})
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def _G1(self):
        return None

    @property
    def _G2(self):
        return None

    @property
    def phaseint(self):
        """Phase integral, q
        Based on Muinonen et al. (2010) Eq. 22
        """
        return 0.009082+0.4061*self._G1+0.8092*self._G2

    @property
    def phasecoeff(self):
        """Phase coefficient, k
        Based on Muinonen et al. (2010) Eq. 23
        """
        return -(30*self._G1+9*self._G2)/(5*np.pi*float(self._G1+self._G2))

    @property
    def oe_amp(self):
        """Opposition effect amplitude, :math:`\zeta-1`
        Based on Muinonen et al. (2010) Eq. 24)
        """
        tmp = float(self._G1+self._G2)
        return (1-tmp)/tmp

    class _spline_positive(_spline):
        """
        Define a spline class that clips negative function values
        """

        def __call__(self, x):
            y = super().__call__(x)
            if hasattr(y, '__iter__'):
                y[y < 0] = 0
            else:
                if y < 0:
                    y = 0
            return y

    _phi1v = (np.deg2rad([7.5, 30., 60, 90, 120, 150]),
              [7.5e-1, 3.3486016e-1, 1.3410560e-1,
               5.1104756e-2, 2.1465687e-2, 3.6396989e-3],
              [-1.9098593, -9.1328612e-2])
    _phi1 = _spline_positive(*_phi1v)
    _phi2v = (np.deg2rad([7.5, 30., 60, 90, 120, 150]),
              [9.25e-1, 6.2884169e-1, 3.1755495e-1,
               1.2716367e-1, 2.2373903e-2, 1.6505689e-4],
              [-5.7295780e-1, -8.6573138e-8])
    _phi2 = _spline_positive(*_phi2v)
    _phi3v = (np.deg2rad([0.0, 0.3, 1., 2., 4., 8., 12., 20., 30.]),
              [1., 8.3381185e-1, 5.7735424e-1, 4.2144772e-1, 2.3174230e-1,
               1.0348178e-1, 6.1733473e-2, 1.6107006e-2, 0.],
              [-1.0630097, 0])
    _phi3 = _spline_positive(*_phi3v)


class HG1G2(HG12BaseClass):
    """HG1G2 photometric phase model (Muinonen et al. 2010)

    Examples
    --------

    >>> # Define the phase function for Themis with
    >>> # H = 7.063, G1 = 0.62, G2 = 0.14
    >>>
    >>> import astropy.units as u
    >>> from sbpy.calib import solar_fluxd
    >>> from sbpy.photometry import HG1G2
    >>> themis = HG1G2(7.063 * u.mag, 0.62, 0.14, radius = 100 * u.km,
    ...     wfb = 'V')
    >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
    ...     print('geometric albedo = {0:.4f}'.format(themis.geomalb))
    ...     print('phase integral = {0:.4f}'.format(themis.phaseint))
    geometric albedo = 0.0656
    phase integral = 0.3742
    """

    H = Parameter(description='H parameter', default=8)
    G1 = Parameter(description='G1 parameter', default=0.2)
    G2 = Parameter(description='G2 parameter', default=0.2)

    @G1.validator
    def G1(self, value):
        """Validate parameter G1 to avoid non-monotonic phase function

        If G1 < 0 or G2 < 0 or G1 + G2 > 1, the phase function could
        potentially be non-monotoic, and a warning will be issued.
        """
        if np.any(value < 0) or np.any(value + self.G2 > 1):
            warnings.warn(
                'G1, G2 parameter combination might result in a non-monotonic'
                ' phase function', NonmonotonicPhaseFunctionWarning)

    @G2.validator
    def G2(self, value):
        """Validate parameter G1 to avoid non-monotonic phase function

        If G1 < 0 or G2 < 0 or G1 + G2 > 1, the phase function could
        potentially be non-monotoic, and a warning will be issued.
        """
        if np.any(value < 0) or np.any(value + self.G1 > 1):
            warnings.warn(
                'G1, G2 parameter combination might result in a non-monotonic'
                ' phase function', NonmonotonicPhaseFunctionWarning)

    @property
    def _G1(self):
        return self.G1.value

    @property
    def _G2(self):
        return self.G2.value

    @staticmethod
    def evaluate(ph, h, g1, g2):
        func = g1*HG1G2._phi1(ph)+g2*HG1G2._phi2(ph)+(1-g1-g2)*HG1G2._phi3(ph)
        if isinstance(func, u.Quantity):
            func = func.value
        func = -2.5 * np.log10(func)
        if isinstance(h, u.Quantity):
            func = func * h.unit
        return h + func

    @staticmethod
    def fit_deriv(ph, h, g1, g2):
        if hasattr(ph, '__iter__'):
            ddh = np.ones_like(ph)
        else:
            ddh = 1.
        phi1 = HG1G2._phi1(ph)
        phi2 = HG1G2._phi2(ph)
        phi3 = HG1G2._phi3(ph)
        dom = (g1*phi1+g2*phi2+(1-g1-g2)*phi3)
        ddg1 = 1.085736205*(phi3-phi1)/dom
        ddg2 = 1.085736205*(phi3-phi2)/dom
        return [ddh, ddg1, ddg2]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict([('H', outputs_unit['y']),
                            ('G1', u.dimensionless_unscaled),
                            ('G2', u.dimensionless_unscaled)])


class HG12(HG12BaseClass):
    """HG12 photometric phase model (Muinonen et al. 2010)

    This system is adopted by IAU as the "standard" model for disk-integrated
    phase functions of planetary objects.  Note that there is a discontinuity
    in the derivative for parameter G12, sometimes making the model fitting
    difficult.  Penttil\"a et al. (2016, Planet. Space Sci. 123, 117-125)
    revised the H, G12 system such that the G12 parameter has a continuous
    derivative.  The revised model is implemented in class `G12_Pen16`.

    Examples
    --------

    >>> # Define the phase function for Themis with
    >>> # H = 7.121, G12 = 0.68
    >>>
    >>> import astropy.units as u
    >>> from sbpy.calib import solar_fluxd
    >>> from sbpy.photometry import HG12
    >>> themis = HG12(7.121 * u.mag, 0.68, radius = 100 * u.km, wfb = 'V')
    >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
    ...     print('geometric albedo = {0:.4f}'.format(themis.geomalb))
    ...     print('phase integral = {0:.4f}'.format(themis.phaseint))
    geometric albedo = 0.0622
    phase integral = 0.3949

    """

    H = Parameter(description='H parameter', default=8)
    G12 = Parameter(description='G12 parameter', default=0.3)

    @G12.validator
    def G12(self, value):
        """Validate parameter G12 to avoid non-monotonic phase function

        If G12 < -0.70 or G12 > 1.30, the phase function could potentially be
        non-monotoic, and a warning will be issued.
        """
        if np.any(value < -0.70) or np.any(value > 1.30):
            warnings.warn(
                'G12 parameter could result in a non-monotonic phase function',
                NonmonotonicPhaseFunctionWarning)

    @property
    def _G1(self):
        return self._G12_to_G1(self.G12.value)

    @property
    def _G2(self):
        return self._G12_to_G2(self.G12.value)

    @staticmethod
    def _G12_to_G1(g12):
        """Calculate G1 from G12"""
        if g12 < 0.2:
            return 0.7527*g12+0.06164
        else:
            return 0.9529*g12+0.02162

    @staticmethod
    def _G12_to_G2(g12):
        """Calculate G2 from G12"""
        if g12 < 0.2:
            return -0.9612*g12+0.6270
        else:
            return -0.6125*g12+0.5572

    @staticmethod
    def evaluate(ph, h, g):
        g1 = HG12._G12_to_G1(g)
        g2 = HG12._G12_to_G2(g)
        return HG1G2.evaluate(ph, h, g1, g2)

    @staticmethod
    def fit_deriv(ph, h, g):
        if hasattr(ph, '__iter__'):
            ddh = np.ones_like(ph)
        else:
            ddh = 1.
        g1 = HG12._G12_to_G1(g)
        g2 = HG12._G12_to_G2(g)
        phi1 = HG1G2._phi1(ph)
        phi2 = HG1G2._phi2(ph)
        phi3 = HG1G2._phi3(ph)
        dom = (g1*phi1+g2*phi2+(1-g1-g2)*phi3)
        if g < 0.2:
            p1 = 0.7527
            p2 = -0.9612
        else:
            p1 = 0.9529
            p2 = -0.6125
        ddg = 1.085736205*((phi3-phi1)*p1+(phi3-phi2)*p2)/dom
        return [ddh, ddg]

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict([('H', outputs_unit['y']),
                            ('G12', u.dimensionless_unscaled)])


class HG12_Pen16(HG12):
    """Revised H, G12 model by Penttil\"a et al. (2016)

    This system is the revised H, G12 system by Penttil\"a et al. (2016,
    Planet. Space Sci. 123, 117-125) that has a continuous derivative with
    respect to parameter G12.  The original model as adopted by IAU as the
    "standard" model for disk-integrated phase functions of planetary objects
    is implemented in class `HG12`.

    Examples
    --------
    >>> # Define the phase function for Themis with
    >>> # H = 7.121, G12 = 0.68
    >>>
    >>> import astropy.units as u
    >>> from sbpy.calib import solar_fluxd
    >>> from sbpy.photometry import HG12_Pen16
    >>> themis = HG12_Pen16(7.121 * u.mag, 0.68, radius = 100 * u.km,
    ...     wfb = 'V')
    >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
    ...     print('geometric albedo = {0:.4f}'.format(themis.geomalb))
    ...     print('phase integral = {0:.4f}'.format(themis.phaseint))
    geometric albedo = 0.0622
    phase integral = 0.3804
    """

    @cite({'definition': '2016P&SS..123..117P'})
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def _G12_to_G1(g12):
        """Calculate G1 from G12"""
        return 0.84293649*g12

    @staticmethod
    def _G12_to_G2(g12):
        """Calculate G2 from G12"""
        return 0.53513350*(1-g12)

    @staticmethod
    def fit_deriv(ph, h, g):
        if hasattr(ph, '__iter__'):
            ddh = np.ones_like(ph)
        else:
            ddh = 1.
        g1 = HG12_Pen16._G12_to_G1(g)
        g2 = HG12_Pen16._G12_to_G2(g)
        phi1 = HG1G2._phi1(ph)
        phi2 = HG1G2._phi2(ph)
        phi3 = HG1G2._phi3(ph)
        dom = (g1*phi1+g2*phi2+(1-g1-g2)*phi3)
        p1 = 0.84293649
        p2 = -0.53513350
        ddg = 1.085736205*((phi3-phi1)*p1+(phi3-phi2)*p2)/dom
        return [ddh, ddg]
