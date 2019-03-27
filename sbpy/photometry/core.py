# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Photometry Module

created on June 23, 2017

"""

__all__ = ['ref2mag', 'mag2ref', 'spline',
           'DiskIntegratedPhaseFunc', 'LinearPhaseFunc', 'HG', 'HG12', 'HG1G2', 'HG12_Pen16',
           'DiskFunctionModel', 'LommelSeeliger', 'Lambert', 'LunarLambert',
           'PhaseFunctionModel', 'ROLOPhase',
           'ResolvedPhotometricModelClass', 'ROLO']

from collections import OrderedDict
import numpy as np
from numbers import Number
from scipy.integrate import quad
from astropy.modeling import (FittableModel, Fittable1DModel,
                              Fittable2DModel, Parameter)
from astropy import units as u
from ..data import Ephem


def _process_ephem_input(eph, key=None):
    if eph is None:
        return None, None
    if not isinstance(eph, (Ephem, dict, np.ndarray, Number, u.Quantity)):
        raise ValueError('`~sbpy.data.Ephem`, `dict`, `numpy.ndarray`, `Number`, or `~astropy.units.Quantity expected, {0} received'.format(type(eph)))
    if isinstance(eph, dict):
        eph = Ephem(eph)
    if isinstance(eph, Ephem):
        if key is None:
            out = None
        else:
            out = eph[key]
    else:
        out = eph
        eph = None
    return eph, out


def ref2mag(ref, radius, M_sun=None):
    """Convert average bidirectional reflectance to reduced magnitude

    Parameters
    ----------
    ref : float, astropy.units.Quantity
        Average bidirectional reflectance
    radius : float, astropy.units.Quantity
        Radius of object
    M_sun : float, optional
        The magnitude of the Sun, default is -26.74

    Returns
    -------
    The same type as input, reduced magnitude (at r=delta=1 au)

    Examples
    --------
    >>> from astropy import units as u
    >>> from sbpy.photometry import ref2mag
    >>> mag = ref2mag(0.1, 460)
    >>> print('{0:.4}'.format(mag))
    2.078
    >>> mag = ref2mag(0.1, 460*u.km)
    >>> print('{0:.4}'.format(mag))
    2.078 mag
    """

    if M_sun is None:
        M_sun = -26.74
    Q = False
    if isinstance(ref, u.Quantity):
        ref = ref.value
        Q = True
    if isinstance(radius, u.Quantity):
        radius = radius.to('km').value
        Q = True
    if isinstance(M_sun, u.Quantity):
        M_sun = M_sun.to('mag').value
        Q = True

    mag = M_sun-2.5*np.log10(ref*np.pi*radius*radius*u.km.to('au')**2)
    if Q:
        return mag*u.mag
    else:
        return mag


def mag2ref(mag, radius, M_sun=None):
    """Convert reduced magnitude to average bidirectional reflectance

    Parameters
    ----------
    mag : float, astropy.units.Quantity
        Reduced magnitude
    radius : float, astropy.units.Quantity
        Radius of object
    M_sun : float, optional
        The magnitude of the Sun, default is -26.74

    Returns
    -------
    The same type as input, average bidirectional reflectance

    Examples
    --------
    >>> from astropy import units as u
    >>> from sbpy.photometry import mag2ref
    >>> ref = mag2ref(2.08, 460)
    >>> print('{0:.4}'.format(ref))
    0.09981
    >>> ref = mag2ref(2.08, 460*u.km)
    >>> print('{0:.4}'.format(ref))
    0.09981 1 / sr
    """

    if M_sun is None:
        M_sun = -26.74
    Q = False
    if isinstance(mag, u.Quantity):
        mag = mag.value
        Q = True
    if isinstance(radius, u.Quantity):
        radius = radius.to('km').value
        Q = True
    if isinstance(M_sun, u.Quantity):
        M_sun = M_sun.to('mag').value
        Q = True

    ref = 10**((M_sun-mag)*0.4)/(np.pi*radius*radius*u.km.to('au')**2)
    if Q:
        return ref/u.sr
    else:
        return ref


class spline(object):
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


class DiskIntegratedPhaseFunc(Fittable1DModel):
    """Base class for disk-integrated phase function model

    Examples
    --------
    Define a linear phase function with phase slope 0.04 mag/deg, and
    study its properties

    >>> # Define a disk-integrated phase function model
    >>> import numpy as np
    >>> from astropy.modeling import Parameter
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
    >>> linear_phasefunc = LinearPhaseFunc(5, 2.29, radius=300)
    >>> pha = np.linspace(0, 180, 200)
    >>> mag = linear_phasefunc.mag(np.deg2rad(pha))
    >>> ref = linear_phasefunc.ref(np.deg2rad(pha))
    >>> geoalb = linear_phasefunc.geoalb
    >>> phaseint = linear_phasefunc.phaseint
    >>> bondalb = linear_phasefunc.bondalb
    >>> print('Geometric albedo is {0:.3}'.format(geoalb))
    Geometric albedo is 0.0501
    >>> print('Bond albedo is {0:.3}'.format(bondalb))
    Bond albedo is 0.0184
    >>> print('Phase integral is {0:.3}'.format(phaseint))
    Phase integral is 0.368
    """

    _unit = None
    input_units = {'x': u.rad}
    input_units_allow_dimensionless = {'x': True}

    def __init__(self, *args, radius=None, M_sun=None, **kwargs):
        """Initialize DiskIntegratedPhaseFunc

        Parameters
        ----------
        radius : float or astropy.units.Quantity
            Radius of object, in km if a float.  Needed for magnitude and
            reflectance conversion
        M_sun : float or astropy.units.Quantity
            Solar magnitude.  Needed for magnitude and reflectance conversion.
            If not supplied, the V-band solar magnitude assumed.  See
            `~ref2mag` and `~mag2ref`

        """
        super().__init__(*args, **kwargs)
        self.radius = radius
        self.M_sun = M_sun

    def _check_unit(self):
        if self._unit is None:
            raise ValueError('the unit of phase function is unknown')

    @property
    def geoalb(self):
        """Geometric albedo"""
        alb = np.pi*self.ref(0)['ref'][0]
        if hasattr(alb, 'unit'):
            alb = alb*u.sr
        return alb

    @property
    def bondalb(self):
        """Bond albedo"""
        return self.geoalb*self.phaseint

    @property
    def phaseint(self):
        """Phase integral"""
        return self._phase_integral()

    def fit(self, eph, mag, fitter=None, return_fitter=False, **kwargs):
        """Fit disk-integrated phase function model to magnitude data and the
        corresponding ephemerides data.

        This is a wrapper that provides a consistent interface that is
        compatible with the `sbpy` guideline.

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`, dict_like, or array_like of float
            If `~sbpy.data.Ephem` or dict_like, must contain a column with a
            name that's compatible with the alternative names for phase angle
            (see `~sbpy.data.DataClass`).  If array_like, then it is the phase
            angles of object.  If any distance (heliocentric and geocentric)
            is provided, then they will be used to correct magnitude to
            1 au before fitting.  If no unit is provided via type
            `~astropy.units.Quantity`, then radiance is assumed for phase
            angle, and au is assumed for distances.
        mag : array_like, `~astropy.units.Quantity`
            Magnitude data to be fitted
        fitter : `~astropy.modeling.fitting.Fitter` class, optional
            The fitter class to be used for fitting.  Default:
            `~astropy.modeling.fitting.LevMarLSQFitter`
        return_fitter : bool, optional
            If `True`, then the fitter class is returned in a tuple.
            Default: `False`
        **kwargs : Keyword parameters accepted by `fitter.__call__()`.

        Returns
        -------
        Object of `DiskIntegratedPhaseFunc` subclass
            The best-fit model class object.

        Examples
        --------
        >>> from sbpy.photometry import HG # doctest: +SKIP
        >>> from sbpy.data import Misc # doctest: +SKIP
        >>> eph = Misc.mpc_observations('Bennu') # doctest: +SKIP
        >>> hg = HG() # doctest: +SKIP
        >>> best_hg = hg.fit(eph, eph['mag']) # doctest: +SKIP
        """
        eph, pha = _process_ephem_input(eph, 'alpha')

        mag = np.asanyarray(mag)
        dist_corr = self._distance_module(eph)
        if isinstance(mag, u.Quantity):
            dist_corr = dist_corr * u.mag
        mag0 = mag + dist_corr

        if fitter is None:
            from astropy.modeling.fitting import LevMarLSQFitter, SLSQPLSQFitter
            fitter = LevMarLSQFitter
        fit = fitter()

        model = fit(self, pha, mag0, **kwargs)
        if return_fitter:
            return model, fit
        else:
            return model

    def _distance_module(self, eph):
        """Return the correction magnitude or factor for heliocentric distance
        and observer distance

        Parameters
        ----------
        eph : any type
            If `~sbpy.data.Ephem` or dict_like, then it's used to calculate
            distanc correction.  If no unit is provided via type
            `~astropy.units.Quantity`, then the distance is assumed to be in
            unit of au.  For any other data type, a factor 1 a magnitude of 0
            will be returned

        Returns
        -------
        float or numpy array
            Magnitude correction to be added or scaling factor to be
            multiplied to the apparent magnitude or flux, respectively, in
            order to correct to heliocentric distance and observer distance of
            both 1 au.
        """
        eph, dummy = _process_ephem_input(eph)
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
        module = np.asarray(module)
        if self._unit == 'mag':
            return -2.5 * np.log10(module)
        else:
            return module

    @classmethod
    def from_data(cls, eph, mag, init=None, fitter=None, **kwargs):
        """Instantiate a photometric model class object from data

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`, dict_like, or array_like of float
            Ephemerides of data.  See `.fit()` for more details.
        mag : array_like, `~astropy.units.Quantity`
            Magnitude data to be fitted
        init : array_like, `~astropy.units.Quantity`
            The intial guess parameters
        fitter : `~astropy.modeling.fitting.Fitter` class, optional
            The fitter class to be used for fitting.  See `.fit()` for
            details.
        **kwargs : Keyword parameters accepted by `fitter.__call__()`.

        Returns
        -------
        `~DiskIntegratedPhaseFunc` class object

        Examples
        --------
        >>> from sbpy.photometry import HG # doctest: +SKIP
        >>> from sbpy.data import Misc # doctest: +SKIP
        >>> eph = Misc.mpc_observations('Bennu') # doctest: +SKIP
        >>> hg = HG.from_data(eph, eph['mag']) # doctest: +SKIP
        """
        if init is None:
            m0 = cls()
        else:
            m0 = cls(*init)
        return m0.fit(eph, mag, fitter=fitter, **kwargs)

    def mag(self, eph, **kwargs):
        """Calculate phase function in magnitude

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`, dict_like, float, or array_like of float
            If `~sbpy.data.Ephem` or dict_like, ephemerides of the object that
            can include phase angle, heliocentric and geocentric distances via
            keywords `phase`, `r` and `delta`.  If float or array_like, then
            the phase angle of object.  If any distance (heliocentric and
            geocentric) is not provided, then it will be assumed to be 1 au.
            If no unit is provided via type `~astropy.units.Quantity`, then
            radiance is assumed for phase angle, and au is assumed for
            distances.

        Returns
        -------
        `~sbpy.data.Ephem` object, with the calculated magnitudes added in
        a new column.

        Examples
        --------
        >>> import numpy as np
        >>> from astropy import units as u
        >>> from sbpy.photometry import HG
        >>> from sbpy.data import Ephem
        >>> ceres_hg = HG(3.34, 0.12)
        >>> # parameter `eph` as `~sbpy.data.Ephem` type
        >>> eph = Ephem({'alpha': np.linspace(0,np.pi*0.9,200),
        ...              'r': np.repeat(2.7*u.au, 200),
        ...              'delta': np.repeat(1.8*u.au, 200)})
        >>> mag1 = ceres_hg.mag(eph)
        >>> # parameter `eph` as numpy array
        >>> pha = np.linspace(0, 180, 200)
        >>> mag2 = ceres_hg.mag(np.deg2rad(pha))
        """
        self._check_unit()
        eph, pha = _process_ephem_input(eph, 'alpha')
        out = self(pha, **kwargs)
        if self._unit != 'mag':
            if self.radius is None:
                raise ValueError(
                    'cannot calculate phase funciton in magnitude because the size of object is unknown')
            out = ref2mag(out, self.radius, M_sun=self.M_sun)
        else:
            dist_corr = self._distance_module(eph)
            if isinstance(out, u.Quantity):
                dist_corr = dist_corr*u.mag
            out = out - dist_corr
        if eph is None:
            eph = Ephem({'alpha': pha, 'mag': out})
        else:
            eph.add_column(out, name='mag')
        return eph

    def ref(self, eph, normalized=None, **kwargs):
        """Calculate phase function in average bidirectional reflectance

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`, dict_like, float, or array_like of float
            If `~sbpy.data.Ephem` or dict_like, ephemerides of the object that
            can include phase angle, heliocentric and geocentric distances via
            keywords `phase`, `r` and `delta`.  If float or array_like, then
            the phase angle of object.  If any distance (heliocentric and
            geocentric) is not provided, then it will be assumed to be 1 au.
            If no unit is provided via type `~astropy.units.Quantity`, then
            radiance is assumed for phase angle, and au is assumed for
            distances.

        Returns
        -------
        `~sbpy.data.Ephem` object, with the calculated reflectances added in a
        new column.

        Examples
        --------
        >>> import numpy as np
        >>> from astropy import units as u
        >>> from sbpy.photometry import HG
        >>> from sbpy.data import Ephem
        >>> ceres_hg = HG(3.34, 0.12, radius=480)
        >>> # parameter `eph` as `~sbpy.data.Ephem` type
        >>> eph = Ephem({'alpha': np.linspace(0,np.pi*0.9,200),
        ...              'r': np.repeat(2.7*u.au, 200),
        ...              'delta': np.repeat(1.8*u.au, 200)})
        >>> ref1 = ceres_hg.ref(eph)
        >>> # parameter `eph` as numpy array
        >>> pha = np.linspace(0, 180, 200)
        >>> ref2 = ceres_hg.mag(np.deg2rad(pha))
        """
        self._check_unit()
        eph, pha = _process_ephem_input(eph, 'alpha')
        out = self(pha, **kwargs)
        if normalized is not None:
            norm = self(normalized, **kwargs)
        if self._unit == 'ref':
            if normalized is not None:
                out /= norm
        else:
            if self.radius is None:
                raise ValueError(
                    'cannot calculate phase function in reflectance unit because the size of object is unknown')
            out = mag2ref(out, self.radius, M_sun=self.M_sun)
            if normalized is not None:
                out /= mag2ref(norm, self.radius, M_sun=self.M_sun)
        if eph is None:
            eph = Ephem({'alpha': pha, 'ref': out})
        else:
            eph.add_column(out, name='ref')
        return eph

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
        >>> from sbpy.photometry import HG
        >>> ceres_hg = HG(3.34, 0.12, radius=480)
        >>> print('{0:.3}'.format(ceres_hg._phase_integral()))
        0.364

        """
        def integrand(x): return 2*self.ref(x, normalized=0.)['ref']*np.sin(x)
        return integrator(integrand, 0, np.pi)[0]


class LinearPhaseFunc(DiskIntegratedPhaseFunc):
    """Linear phase function model

    Examples
    --------
    >>> # Define a linear phase function model with absolute magnitude
    >>> # H = 5 and slope = 0.04 mag/deg = 2.29 mag/rad
    >>> from sbpy.photometry import LinearPhaseFunc
    >>>
    >>> linear_phasefunc = LinearPhaseFunc(5, 2.29, radius=300)
    >>> pha = np.linspace(0, 180, 200)
    >>> mag = linear_phasefunc.mag(np.deg2rad(pha))
    >>> ref = linear_phasefunc.ref(np.deg2rad(pha))
    >>> geoalb = linear_phasefunc.geoalb
    >>> phaseint = linear_phasefunc.phaseint
    >>> bondalb = linear_phasefunc.bondalb
    >>> print('Geometric albedo is {0:.3}'.format(geoalb))
    Geometric albedo is 0.0501
    >>> print('Bond albedo is {0:.3}'.format(bondalb))
    Bond albedo is 0.0184
    >>> print('Phase integral is {0:.3}'.format(phaseint))
    Phase integral is 0.368

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
    >>> from sbpy.photometry import HG
    >>> ceres = HG(3.34, 0.12, radius=480)
    >>> print('{0:.4f}'.format(ceres.geoalb))
    0.0902
    >>> print('{0:.4f}'.format(ceres.phaseint))
    0.3644

    """

    _unit = 'mag'
    H = Parameter(description='H parameter', default=8)
    G = Parameter(description='G parameter', default=0.4)

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
            func = func * u.mag
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

    class _spline_positive(spline):
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
    >>> from sbpy.photometry import HG1G2
    >>> themis = HG1G2(7.063,0.62,0.14,radius=100)
    >>> print('{0:.4f}'.format(themis.geoalb))
    0.0674
    >>> print('{0:.4f}'.format(themis.phaseint))
    0.3742

    """

    H = Parameter(description='H parameter', default=8)
    G1 = Parameter(description='G1 parameter', default=0.2)
    G2 = Parameter(description='G2 parameter', default=0.2)

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
            func = func * u.mag
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

    Examples
    --------

    >>> # Define the phase function for Themis with
    >>> # H = 7.121, G12 = 0.68
    >>>
    >>> from sbpy.photometry import HG12
    >>> themis = HG12(7.121, 0.68, radius=100)
    >>> print('{0:.4f}'.format(themis.geoalb))
    0.0639
    >>> print('{0:.4f}'.format(themis.phaseint))
    0.3949

    """

    H = Parameter(description='H parameter', default=8)
    G12 = Parameter(description='G12 parameter', default=0.3)

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

    Examples
    --------
    >>> # Define the phase function for Themis with
    >>> # H = 7.121, G12 = 0.68
    >>>
    >>> from sbpy.photometry import HG12
    >>> themis = HG12_Pen16(7.121, 0.68, radius=100)
    >>> print('{0:.4f}'.format(themis.geoalb))
    0.0639
    >>> print('{0:.4f}'.format(themis.phaseint))
    0.3804

    """

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


class DiskFunctionModel(FittableModel):
    """Base class for disk-function model"""
    pass


class LommelSeeliger(DiskFunctionModel):
    """Lommel-Seeliger model class"""

    inputs = ('i', 'e')
    outputs = ('d',)

    @staticmethod
    def evaluate(i, e):
        mu0 = np.cos(i)
        mu = np.cos(e)
        return mu0/(mu0+mu)


class Lambert(DiskFunctionModel):
    """Lambert model class"""

    inputs = ('i',)
    outputs = ('d',)

    @staticmethod
    def evaluate(i):
        return np.cos(i)


class LunarLambert(DiskFunctionModel):
    """Lunar-Lambert model, or McEwen model class"""
    inputs = ('i', 'e')
    outputs = ('d',)

    L = Parameter(default=0.1, description='Partition parameter')

    @staticmethod
    def evaluate(i, e, L):
        return (1-L) * LommelSeeliger.evaluate(i, e) + L * Lambert.evaluate(i)


class PhaseFunctionModel(FittableModel):
    """Base class for phase function model"""
    inputs = ('a',)
    outputs = ('f',)


class ROLOPhase(PhaseFunctionModel):
    """ROLO phase function model class"""
    A0 = Parameter(default=0.1, min=0., description='ROLO A0 parameter')
    A1 = Parameter(default=0.1, min=0., description='ROLO A1 parameter')
    C0 = Parameter(default=0.1, min=0., description='ROLO C0 parameter')
    C1 = Parameter(default=0.1, description='ROLO C1 parameter')
    C2 = Parameter(default=0.1, description='ROLO C2 parameter')
    C3 = Parameter(default=0.1, description='ROLO C3 parameter')
    C4 = Parameter(default=0.1, description='ROLO C4 parameter')

    @staticmethod
    def evaluate(pha, c0, c1, a0, a1, a2, a3, a4):
        pha2 = pha*pha
        return c0*np.exp(-c1*pha)+a0+a1*pha+a2*pha2+a3*pha*pha2+a4*pha2*pha2

    @staticmethod
    def fit_deriv(pha, c0, c1, a0, a1, a2, a3, a4):
        pha2 = pha*pha
        dc0 = np.exp(-c1*pha)
        if hasattr(pha, '__iter__'):
            dda = np.ones(len(pha))
        else:
            dda = 1.
        return [dc0, -c0*c1*dc0, dda, pha, pha2, pha*pha2, pha2*pha2]


class ResolvedPhotometricModelClass(object):
    """Base class for disk-resolved photometric model"""
    # composite model as the product of a disk function and a phase function
    pass


class ROLO(ResolvedPhotometricModelClass):
    """ROLO disk-resolved photometric model"""
    pass
