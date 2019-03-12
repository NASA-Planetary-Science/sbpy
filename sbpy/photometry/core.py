# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Photometry Module

created on June 23, 2017

"""

__all__ = ['ref2mag', 'mag2ref', 'spline', 'SpectralGradient',
           'DiskIntegratedModelClass', 'LinearPhaseFunc', 'HG',
           'HG12', 'HG1G2', 'DiskFunctionModel', 'LommelSeeliger',
           'Lambert', 'LunarLambert', 'PhaseFunctionModel', 'ROLOPhase',
           'ResolvedPhotometricModelClass', 'ROLO']

from warnings import warn

import numpy as np
from scipy.integrate import quad
from astropy.modeling import (FittableModel, Fittable1DModel,
                              Fittable2DModel, Parameter)
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from astropy import log

try:
    import synphot
except ImportError:
    synphot = None
    raise warn(AstropyWarning('synphot is not present, photometric filter'
                              ' conversions will be limited.'))

from ..data import Ephem
from ..spectroscopy.sun import Sun
from ..spectroscopy.vega import Vega
from ..units import spectral_density_vega


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
    mag : float, astropy.u.Quantity
        Reduced magnitude
    radius : float, astropy.u.Quantity
        Radius of object
    M_sun : float, optional
        The magnitude of the Sun, default is -26.74

    Returns
    -------
    The same type as input, average bidirectional reflectance

    Examples
    --------
    >>> from astropy import u as u
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


class SpectralGradient(u.SpecificTypeQuantity):
    """Convert between magnitude and spectral gradient.

    ``SpectralGradient`` is a `~astropy.units.Quantity` object with
    units of percent change per wavelength.

    Spectral gradient is with respect to the flux density at a
    particular wavelength.  By convention this is the mean flux
    between the two filters, assuming the reflectance spectrum is
    linear with wavelength.  Eq. 1 of A'Hearn et al. (1984, AJ 89,
    579):

                R(λ1) - R(λ0)  2
        slope = ------------- ---
                R(λ1) + R(λ0) Δλ

    Δλ is typically measured in units of 100 nm.


    Parameters
    ----------
    value : number, `~astropy.units.Quantity`
        The value(s).

    unit : string, `~astropy.units.Unit`, optional
        The unit of the input value.  Strings must be parseable by
        :mod:`~astropy.units` package.

    wave : `~astropy.units.Quantity`, optional
        Effective wavelengths of the measurement for a solar spectral
        energy distribution.  Required for conversion to other
        wavelengths.

    wave0 : `~astropy.units.Quantity`, optional
        Normalization point.  If ``None``, the mean of ``wave``
        will be assumed.

    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`.

    copy : bool, optional
        See `~astropy.units.Quantity`.


    Examples
    --------
    >>> import astropy.units as u
    >>> from sbpy.units import hundred_nm
    >>> S = SpectralGradient(10 * u.percent / hundred_nm, wave0=5500 * u.AA)
    >>> print(S)
    10.0 % / 100 nm

    >>> from sbpy.units import VEGAMAG
    >>> bp = ('johnson_v', 'cousins_r')
    >>> VmR = 15.8 * VEGAMAG - 15.3 * VEGAMAG
    >>> VmR_sun = 0.37 * u.mag
    >>> S = SpectralGradient.from_color(bp, VmR - VmR_sun)
    ...                              # doctest: +REMOTE_DATA +IGNORE_OUTPUT
    >>> S.value                      # doctest: +REMOTE_DATA +FLOAT_CMP
    12.29185986266534
    >>> print(S.unit)                # doctest: +REMOTE_DATA
    % / 100 nm
    """

    _equivalent_unit = u.meter**-1
    _include_easy_conversion_members = False

    def __new__(cls, value, unit=None, wave=None, wave0=None,
                dtype=None, copy=None):
        S = super().__new__(cls, value, unit=unit, dtype=dtype,
                            copy=copy)

        if wave is not None:
            if np.size(wave) != 2:
                raise ValueError(
                    'Two wavelengths must be provided, got {}'
                    .format(np.size(wave)))
            S.wave = S._eff_wave(wave)

        if wave0 is None and wave is not None:
            S.wave0 = S.wave.mean()
        else:
            S.wave0 = wave0

        return S

    @classmethod
    def _eff_wave(cls, wfb):
        """Wavelength/frequency/bandpass to wavelength.

        Bandpass is converted to effective wavelength using a solar
        spectrum.

        """

        eff_wave = (0, 0) * u.um
        sun = Sun.from_default()
        for i in range(2):
            if isinstance(wfb, u.Quantity):
                eff_wave[i] = wfb[i].to(u.um)
            else:
                eff_wave[i] = sun.filt(wfb[i])[0]
                log.info('Using λ_eff = {}'.format(eff_wave[i]))

        return eff_wave

    @classmethod
    def from_color(cls, wfb, color):
        """Initialize from observed color.


        Parameters
        ----------
        wfb : two-element `~astropy.units.Quantity` or tuple
            Wavelengths, frequencies, or bandpasses of the
            measurement.  If a bandpass, the effective wavelength of a
            solar spectrum will be used.  Bandpasses may be a string
            (name) or `~synphot.SpectralElement` (see
            :func:`~sbpy.spectroscopy.sun.Sun.filt`).

        color : `~astropy.units.Quantity`, optional
            Observed color, ``blue - red`` for magnitudes, ``red
            / blue`` for linear units.  Must be dimensionless and have
            the solar color removed.


        Notes
        -----
        Computes spectral gradient from ``color_index``.
        ``wfb[0]`` is the blue-ward of the two measurements:

                  R(λ1) - R(λ0)  2    α - 1  2
          slope = ------------- --- = ----- ---
                  R(λ1) + R(λ0) Δλ    α + 1 Δλ

        where R(λ) is the reflectivity, and

          α = R(λ1) / R(λ0) = 10**(0.4 * color_index)
          color_index = Δm - C_sun

        Δλ is typically expressed in units of 100 nm.


        Examples
        --------
        >>> import astropy.units as u
        >>> w = [0.4719, 0.6185] * u.um
        >>> S = SpectralGradient.from_color(w, 0.10 * u.mag)
        >>> S.value                             # doctest: +FLOAT_CMP
        6.27819572
        >>> print(S.unit)
        % / 100 nm

        """

        eff_wave = SpectralGradient._eff_wave(wfb)

        try:
            # works for u.Magnitudes and dimensionless u.Quantity
            alpha = u.Quantity(color, u.dimensionless_unscaled)
        except u.UnitConversionError:
            # works for u.mag
            alpha = color.to(u.dimensionless_unscaled, u.logarithmic())

        dw = eff_wave[0] - eff_wave[1]
        S = ((2 / dw * (alpha - 1) / (alpha + 1))
             .to(u.percent / hundred_nm))

        return SpectralGradient(S, wave=eff_wave)

    def to_color(self, wfb):
        """Express as a color index.


        Parameters
        ----------
        wfb : two-element `~astropy.units.Quantity` or tuple
            Wavelengths, frequencies, or bandpasses of the
            measurement.  If a bandpass, the effective wavelength of a
            solar spectrum will be used.  Bandpasses may be a string
            (name) or `~synphot.SpectralElement` (see
            :func:`~sbpy.spectroscopy.sun.Sun.filt`).

        Notes
        -----
              1 + S * Δλ / 2
          α = -------------- 
              1 - S * Δλ / 2

        where S is the spectral gradient at the mean of λ0 and λ1, and

          α = R(λ1) / R(λ0) = 10**(0.4 * color_index)
          color_index = Δm - C_sun

        Δλ is typically expressed in units of 100 nm.


        Returns
        -------
        color : `~astropy.units.Quantity`
            ``blue - red`` color in magnitudes, dimensionless and
            excludes the solar color.


        Examples
        --------
        >>> import astropy.units as u
        >>> from sbpy.units import hundred_nm
        >>> S = SpectralGradient(10 * u.percent / hundred_nm,
        ...                      wave0=0.55 * u.um)
        >>> C = S.to_color((525, 575) * u.nm)
        >>> C.value    # doctest: +FLOAT_CMP
        0.05429812423309064
        >>> print(C.unit)
        mag

        """

        eff_wave = self._eff_wave(wfb)

        S = self.renormalize(eff_wave.mean())
        dw = eff_wave[0] - eff_wave[1]
        beta = (S * dw / 2).decompose()  # dimensionless
        color = ((1 + beta) / (1 - beta)).to(u.mag, u.logarithmic())

        return color

    def renormalize(self, wave0):
        """Re-normalize to another wavelength.

        The slope is linearly extrapolated to the new normalization
        point.  Requires the `wave` attribute to be defined, see
        `~SpectralGradient`.

        Parameters
        ----------
        wave0 : `~astropy.units.Quantity`
            Wavelength.

        Returns
        -------
        S : ``SpectralGradient``
            Renormalized gradient.

        """

        if self.wave0 is None:
            raise ValueError('wave0 attribute must be defined.')

        delta = wave0 - self.wave0
        S0 = 1 + self.to(delta.unit**-1) * delta
        S = self / S0
        S.wave0 = wave0
        return S


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


class DiskIntegratedModelClass(Fittable1DModel):
    """Base class for disk-integrated phase function model

    Examples
    --------
    Define a linear phase function with phase slope 0.04 mag/deg, and
    study its properties

    >>> # Define a disk-integrated phase function model
    >>> import numpy as np
    >>> from astropy.modeling import Parameter
    >>> from sbpy.photometry import DiskIntegratedModelClass
    >>>
    >>> class LinearPhaseFunc(DiskIntegratedModelClass):
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

    def __init__(self, *args, radius=None, M_sun=None, **kwargs):
        """Initialize DiskIntegratedModelClass

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
        alb = np.pi*self.ref(0)
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

    def fit(self, eph):
        """Fit photometric model to photometric data stored in sbpy.data.Ephem
        object

        Parameters
        ----------
        eph : `~sbpy.data.Ephem` instance, mandatory
            photometric data; must contain `phase` (phase angle) and `mag`
            (apparent magnitude) columns; `mag_err` optional

        Returns
        -------
        fit Chi2

        Examples
        --------
        >>> from sbpy.photometry import HG # doctest: +SKIP
        >>> from sbpy.data import Misc # doctest: +SKIP
        >>> eph = Misc.mpc_observations('Bennu') # doctest: +SKIP
        >>> hg = HG() # doctest: +SKIP
        >>> chi2 = hg.fit(eph) # doctest: +SKIP

        not yet implemented

        """

    def distance_module(self, eph):
        """Account magnitudes for distance module (distance from observer,
        distance to the Sun); return modified magnitudes

        Parameters
        ----------
        eph : list or array, mandatory
            phase angles

        Returns
        -------
        sbpy.data.Ephem instance

        Examples
        --------
        TBD

        not yet implemented

        """

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
        Numpy array of magnitude based on the phase function model at specified
        phase angle, heliocentric and geocentric distances

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
        if isinstance(eph, Ephem):
            pass
        elif isinstance(eph, dict):
            eph = Ephem(eph)
        else:
            eph = Ephem({'alpha': eph})
        pha = eph['alpha']
        if isinstance(pha, u.Quantity):
            pha = pha.to('rad').value
        out = self(pha, **kwargs)
        if self._unit != 'mag':
            if self.radius is None:
                raise ValueError(
                    'cannot calculate phase funciton in magnitude because the size of object is unknown')
            out = ref2mag(out, self.radius, M_sun=self.M_sun)
        if 'r' in eph.column_names:
            rh = eph['r']
            if isinstance(rh, u.Quantity):
                rh = rh.to('au').value
            out += 5*np.log10(rh)
        if 'delta' in eph.column_names:
            delta = eph['delta']
            if isinstance(delta, u.Quantity):
                delta = delta.to('au').value
            out += 5*np.log10(delta)
        return out

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
        Numpy array of average bidirectional reflectance based on the phase
        function model at specified phase angle

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
        if isinstance(eph, (dict, Ephem)):
            pha = eph['alpha']
        else:
            pha = eph
        if isinstance(pha, u.Quantity):
            pha = pha.to('rad').value
        out = self(pha, **kwargs)
        if normalized is not None:
            if isinstance(normalized, u.Quantity):
                normalized = normalized.to('rad').value
            norm = self(normalized, **kwargs)
        if self._unit == 'ref':
            if normalized is not None:
                out /= norm
            return out
        else:
            if self.radius is None:
                raise ValueError(
                    'cannot calculate phase function in reflectance unit because the size of object is unknown')
            out = mag2ref(out, self.radius, M_sun=self.M_sun)
            if normalized is not None:
                out /= mag2ref(norm, self.radius, M_sun=self.M_sun)
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
        >>> from sbpy.photometry import HG
        >>> ceres_hg = HG(3.34, 0.12, radius=480)
        >>> print('{0:.3}'.format(ceres_hg._phase_integral()))
        0.364

        """
        def integrand(x): return 2*self.ref(x, normalized=0.)*np.sin(x)
        return integrator(integrand, 0, np.pi)[0]


class LinearPhaseFunc(DiskIntegratedModelClass):
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
    S = Parameter(description='Linear slope (mag/rad)')

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


class HG(DiskIntegratedModelClass):
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
    H = Parameter(description='H parameter')
    G = Parameter(description='G parameter')

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
        return hh-2.5*np.log10((1-gg)*HG._hgphi(pha, 1)+gg*HG._hgphi(pha, 2))

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


class HG12BaseClass(DiskIntegratedModelClass):
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

    H = Parameter(description='H parameter')
    G1 = Parameter(description='G1 parameter')
    G2 = Parameter(description='G2 parameter')

    @property
    def _G1(self):
        return self.G1.value

    @property
    def _G2(self):
        return self.G2.value

    @staticmethod
    def evaluate(ph, h, g1, g2):
        return h-2.5*np.log10(g1*HG1G2._phi1(ph)+g2*HG1G2._phi2(ph)+(1-g1-g2)*HG1G2._phi3(ph))

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

    H = Parameter(description='H parameter')
    G12 = Parameter(description='G12 parameter')

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
        ddg = 1.085736205*((phi3-phi1)*p1+(phi3-phi1)*p2)/dom
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
