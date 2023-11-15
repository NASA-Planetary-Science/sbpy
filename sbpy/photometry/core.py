# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Photometry Module

created on June 23, 2017

"""

__all__ = ['DiskIntegratedPhaseFunc', 'LinearPhaseFunc',
           'InvalidPhaseFunctionWarning']

__doctest_requires__ = {
    ("DiskIntegratedPhaseFunc",
     "DiskIntegratedPhaseFunc._phase_integral",
     "DiskIntegratedPhaseFunc.from_obs",
     "LinearPhaseFunc",
     "LinearPhaseFunc._phase_integral"
     ): ["scipy"],
    ("DiskIntegratedPhaseFunc.from_phys", "DiskIntegratedPhaseFunc.to_phys"): ["astroquery"],
}

from collections import OrderedDict
import numpy as np

try:
    from scipy.integrate import quad
except ImportError:
    quad = None

from astropy.modeling import (Fittable1DModel, Parameter)
import astropy.units as u
from astropy import log
from ..data import (Phys, Obs, Ephem, dataclass_input,
                    quantity_to_dataclass)
from ..units import reflectance
from ..exceptions import SbpyWarning
from ..utils.decorators import requires


class InvalidPhaseFunctionWarning(SbpyWarning):
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
    >>> print(phys['targetname','H','G'])    # doctest: +SKIP
    <QTable length=1>
        targetname       H       G
                        mag
          str17       float64 float64
    ----------------- ------- -------
    1 Ceres (A801 AA)    3.31    0.12
    >>> m = HG.from_phys(phys)                  # doctest: +REMOTE_DATA
    INFO: Model initialized for 1 Ceres (A801 AA). [sbpy.photometry.core]
    >>> print(m)                             # doctest: +SKIP
    Model: HG
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 1
    Parameters:
         H     G
        mag
        ---- ----
        3.31 0.12
    >>> print(m.meta['targetname'])          # doctest: +REMOTE_DATA
    1 Ceres (A801 AA)
    >>> print(m.radius)                      # doctest: +SKIP
    469.7 km
    >>>
    >>> # Initialize from orbital elements pulled from JPL Horizons that also
    >>> # contain the H and G parameters
    >>> elem = Orbit.from_horizons('Ceres')  # doctest: +REMOTE_DATA
    >>> print(elem['targetname','H','G'])    # doctest: +SKIP
    <QTable length=1>
        targetname       H       G
                        mag
          str17       float64 float64
    ----------------- ------- -------
    1 Ceres (A801 AA)    3.33    0.12
    >>> m = HG.from_phys(elem)                    # doctest: +REMOTE_DATA
    INFO: Model initialized for 1 Ceres (A801 AA). [sbpy.photometry.core]
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
        -------
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
        >>> print(phys['targetname','H','G'])   # doctest: +SKIP
        <QTable length=1>
            targetname       H       G
                            mag
              str17       float64 float64
        ----------------- ------- -------
        1 Ceres (A801 AA)    3.31    0.12
        >>> m = HG.from_phys(phys)              # doctest: +REMOTE_DATA
        INFO: Model initialized for 1 Ceres (A801 AA). [sbpy.photometry.core]
        >>> print(m)                            # doctest: +SKIP
        Model: HG
        Inputs: ('x',)
        Outputs: ('y',)
        Model set size: 1
        Parameters:
             H     G
            mag
            ---- ----
            3.31 0.12
        >>> print(m.meta['targetname'])         # doctest: +REMOTE_DATA
        1 Ceres (A801 AA)
        >>> print(m.radius)                     # doctest: +REMOTE_DATA
        469.7 km
        >>>
        >>> # Failed initialization due to the lack of field 'G'
        >>> phys = Phys.from_sbdb('12893')      # doctest: +REMOTE_DATA
        >>> print('G' in phys.field_names)      # doctest: +REMOTE_DATA
        False
        >>> m = HG.from_phys(phys)              # doctest: +SKIP
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
        >>> print(phys['targetname','radius','H','G']) # doctest: +SKIP
        <QTable length=1>
            targetname     radius    H       G
                             km     mag
              str17       float64 float64 float64
        ----------------- ------- ------- -------
        1 Ceres (A801 AA)   469.7    3.31    0.12
        >>> m = HG.from_phys(phys)   # doctest: +REMOTE_DATA
        INFO: Model initialized for 1 Ceres (A801 AA). [sbpy.photometry.core]
        >>> m.wfb = 'V'              # doctest: +REMOTE_DATA
        >>> with solar_fluxd.set({'V': -26.77 * u.mag}):
        ...     p = m.to_phys()      # doctest: +REMOTE_DATA
        >>> print(type(p))           # doctest: +REMOTE_DATA
        <class 'sbpy.data.phys.Phys'>
        >>> p.table.pprint(max_width=-1)  # doctest: +SKIP
            targetname    diameter  H    G            pv                  A
                             km    mag
        ----------------- -------- ---- ---- ------------------- --------------------
        1 Ceres (A801 AA)    939.4 3.31 0.12 0.07624470768627523 0.027779803126557152
        """  # noqa: E501
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
        >>> from astropy.modeling.fitting import SLSQPLSQFitter
        >>> fitter = SLSQPLSQFitter()
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
                dist_corr = -2.5 * np.log10(dist_corr)
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
                    dist_corr1 = -2.5 * np.log10(dist_corr)
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

        When ``append_results == False``: The calculated magnitude will be
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
        >>> eph = Ephem.from_dict({'alpha': np.linspace(
        ...                                     0, np.pi * 0.9, 200) * u.rad,
        ...                        'r': np.repeat(2.7 * u.au, 200),
        ...                        'delta': np.repeat(1.8 * u.au, 200)})
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
                    'Cannot calculate phase function in magnitude because the'
                    ' size of object is unknown.')
            if self.wfb is None:
                raise ValueError('Wavelength/Frequency/Band is unknown.')
            out = out.to(
                unit,
                reflectance(self.wfb, cross_section=np.pi * self.radius**2)
            )
        dist_corr = self._distance_module(eph)
        dist_corr = u.Quantity(dist_corr).to(u.mag, u.logarithmic())
        out = out - dist_corr
        if append_results:
            name = 'mag'
            i = 1
            while name in eph.field_names:
                name = 'mag'+str(i)
                i += 1
            eph.apply(out, name=name)
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

        When ``append_results == False``: The calculated reflectance will be
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
        >>> eph = Ephem.from_dict({'alpha': np.linspace(
        ...                                     0, np.pi * 0.9, 200) * u.rad,
        ...                        'r': np.repeat(2.7 * u.au, 200),
        ...                        'delta': np.repeat(1.8 * u.au, 200)})
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
                out = out.to(
                    '1/sr',
                    reflectance(self.wfb, cross_section=np.pi*self.radius**2)
                )
            else:
                out = out - norm
                out = out.to('', u.logarithmic())
        if append_results:
            name = 'ref'
            i = 1
            while name in eph.field_names:
                name = 'ref'+str(i)
                i += 1
            eph.apply(out, name=name)
            return eph
        else:
            return out

    @requires("scipy")
    def _phase_integral(self, integrator=quad):
        """Calculate phase integral with numerical integration

        Parameters
        ----------
        integrator : function, optional
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
