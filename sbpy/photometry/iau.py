# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy photometry module
IAU HG series phase function

created on July, 2022

"""

import warnings
import numpy as np
from scipy.integrate import quad
from astropy.modeling import Parameter
from .core import DiskIntegratedPhaseFunc, NonmonotonicPhaseFunctionWarning
from ..bib import cite


__all__ = ['HG', 'HG12BaseClass', 'HG12', 'HG1G2', 'HG12_Pen16']


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
