# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy photometry module
IAU HG series phase function

created on July, 2022

"""


__all__ = ['HG', 'HG12', 'HG1G2', 'HG12_Pen16']

__doctest_requires__ = {
    "HG*": ["scipy"],
}

import warnings
from collections import OrderedDict
import numpy as np
from astropy.modeling import Parameter
import astropy.units as u
from .core import DiskIntegratedPhaseFunc, InvalidPhaseFunctionWarning
from ..bib import cite


# define the bounds of various model parameters
_hg_g_bounds = [-0.253, 1.194]
_hg1g2_g1_bounds = [0, 1]
_hg1g2_g2_bounds = [0, 1]
_hg12_g12_bounds = [-0.0818919, 0.909714]
_hg12_pen16_g12_bounds = [0, 1]


def _check_bounds(value, bounds, exception, msg=''):
    """Check if all `value` is within `bounds`

    Parameters
    ----------
    value : number, number array
        Value to be checked
    bounds : 2-element array-like representing [lower, upper]
        Bounds for value.  Lower or upper could be None, representing
        no bounds.
    exception : Warning or Exception
        If out of bounds, issue a warning or exception.
    msg : str, optional
        Exception or warning string
    """
    value = np.asanyarray(value)
    if ((bounds[0] is not None and (value < bounds[0]).any())
        or (bounds[1] is not None and (value > bounds[1]).any())):
        if issubclass(exception, Warning):
            warnings.warn(msg, exception)
        elif issubclass(exception, Exception):
            raise exception(msg)


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
    G = Parameter(description='G parameter', default=0.4,
                  bounds=_hg_g_bounds)

    @cite({'definition': '1989aste.conf..524B'})
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @G.validator
    def G(self, value):
        _check_bounds(value, _hg_g_bounds, InvalidPhaseFunctionWarning,
            'G parameter could result in an invalid phase function')

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
    G1 = Parameter(description='G1 parameter', default=0.2,
                   bounds=_hg1g2_g1_bounds)
    G2 = Parameter(description='G2 parameter', default=0.2,
                   bounds=_hg1g2_g2_bounds)

    def __init__(self, *args, **kwargs):
        ineqcons = kwargs.pop('ineqcons', [])
        ineqcons.append(lambda x, *args: 1 - x[1] - x[2])
        super().__init__(*args, ineqcons=ineqcons, **kwargs)

    @G1.validator
    def G1(self, value):
        _check_bounds(value, _hg1g2_g1_bounds, InvalidPhaseFunctionWarning,
            'G1 results in an invalid phase function.')
        if (value + self.G2 > 1).any():
            warnings.warn('G1, G2 combination results in an invalid '
                'phase function.')

    @G2.validator
    def G2(self, value):
        _check_bounds(value, _hg1g2_g2_bounds, InvalidPhaseFunctionWarning,
            'G2 results in an invalid phase function.')
        if (value + self.G1 > 1).any():
            warnings.warn('G1, G2 combination results in an invalid '
                'phase function.')

    @property
    def _G1(self):
        return self.G1.value

    @property
    def _G2(self):
        return self.G2.value

    @staticmethod
    def evaluate(ph, h, g1, g2):
        ph = u.Quantity(ph, 'rad').to_value('rad')
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
    G12 = Parameter(description='G12 parameter', default=0.3,
                    bounds=_hg12_g12_bounds)

    @G12.validator
    def G12(self, value):
        _check_bounds(value, _hg12_g12_bounds, InvalidPhaseFunctionWarning,
            'G12 parameter could result in an invalid phase function')

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
    def evaluate(ph, h, g12):
        g1 = HG12._G12_to_G1(g12)
        g2 = HG12._G12_to_G2(g12)
        return HG1G2.evaluate(ph, h, g1, g2)

    @staticmethod
    def fit_deriv(ph, h, g12):
        if hasattr(ph, '__iter__'):
            ddh = np.ones_like(ph)
        else:
            ddh = 1.
        g1 = HG12._G12_to_G1(g12)
        g2 = HG12._G12_to_G2(g12)
        phi1 = HG1G2._phi1(ph)
        phi2 = HG1G2._phi2(ph)
        phi3 = HG1G2._phi3(ph)
        dom = (g1*phi1+g2*phi2+(1-g1-g2)*phi3)
        if g12 < 0.2:
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

    G12 = Parameter(description='G12 parameter', default=0.3,
                    bounds=_hg12_pen16_g12_bounds)

    @G12.validator
    def G12(self, value):
        _check_bounds(value, _hg12_pen16_g12_bounds,
                      InvalidPhaseFunctionWarning,
                      'G12 parameter could result in an invalid phsae function')

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
    def evaluate(ph, h, g12):
        g1 = HG12_Pen16._G12_to_G1(g12)
        g2 = HG12_Pen16._G12_to_G2(g12)
        return HG1G2.evaluate(ph, h, g1, g2)

    @staticmethod
    def fit_deriv(ph, h, g12):
        if hasattr(ph, '__iter__'):
            ddh = np.ones_like(ph)
        else:
            ddh = 1.
        g1 = HG12_Pen16._G12_to_G1(g12)
        g2 = HG12_Pen16._G12_to_G2(g12)
        phi1 = HG1G2._phi1(ph)
        phi2 = HG1G2._phi2(ph)
        phi3 = HG1G2._phi3(ph)
        dom = (g1*phi1+g2*phi2+(1-g1-g2)*phi3)
        p1 = 0.84293649
        p2 = -0.53513350
        ddg = 1.085736205*((phi3-phi1)*p1+(phi3-phi2)*p2)/dom
        return [ddh, ddg]
