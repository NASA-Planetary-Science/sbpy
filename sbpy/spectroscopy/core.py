# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Spectroscopy Module

created on June 23, 2017
"""

import numpy as np
import astropy.units as u


__all__ = ["Spectrum", "SpectralModel", "SpectralGradient"]

__doctest_requires__ = {
    "SpectralGradient": ["synphot"],
    "SpectralGradient.from_color": ["synphot"],
    "SpectralGradient.to_color": ["synphot"],
    "SpectralGradient.renormalize": ["synphot"],
}


class SpectralGradient(u.SpecificTypeQuantity):
    r"""Convert between magnitude and spectral gradient.

    ``SpectralGradient`` is a `~astropy.units.Quantity` object with
    units of percent change per wavelength.

    Spectral gradient is with respect to the flux density at a
    particular wavelength.  By convention this is the mean flux
    between the two filters, assuming the reflectance spectrum is
    linear with wavelength.  Eq. 1 of A'Hearn et al. [ADT84]_:

    .. math::

        S = \frac{R(λ1) - R(λ0)}{R(λ1) + R(λ0)} \frac{2}{Δλ}

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
    10.0 % / (100 nm)

    >>> from sbpy.units import VEGAmag
    >>> from sbpy.photometry import bandpass
    >>> V = bandpass('Johnson V')
    >>> R = bandpass('Cousins R')
    >>> VmR = 15.8 * VEGAmag - 15.3 * VEGAmag
    >>> VmR_sun = 0.37 * u.mag
    >>> S = SpectralGradient.from_color((V, R), VmR - VmR_sun)
    >>> print(S)    # doctest: +FLOAT_CMP
    12.29185986266534 % / (100 nm)


    References
    ----------
    .. [ADT84] A'Hearn, Dwek & Tokunaga 1984. Infrared Photometry of
        Comet Bowell and Other Comets. ApJ 282, 803-806.

    """

    _equivalent_unit = u.meter**-1
    _include_easy_conversion_members = False

    def __new__(cls, value, unit=None, wave=None, wave0=None, dtype=None, copy=True):
        S = super().__new__(cls, value, unit=unit, dtype=dtype, copy=copy)

        if wave is not None:
            if np.size(wave) != 2:
                raise ValueError(
                    "Two wavelengths must be provided, got {}".format(np.size(wave))
                )
            S.wave = S._lambda_eff(wave)

        if wave0 is None and wave is not None:
            S.wave0 = S.wave.mean()
        else:
            S.wave0 = wave0

        return S

    @classmethod
    def _lambda_eff(cls, wfb):
        """Wavelength/frequency/bandpass to wavelength.

        Bandpass is converted to effective wavelength using a solar
        spectrum.

        """
        from ..calib import Sun

        lambda_eff, ci = Sun.from_default().color_index(wfb, u.ABmag)

        return lambda_eff

    @classmethod
    def from_color(cls, wfb, color):
        r"""Initialize from observed color.


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
        ``wfb[0]`` is the blue-ward of the two measurements

        .. math::

           S &= \frac{R(λ1) - R(λ0)}{R(λ1) + R(λ0)} \frac{2}{Δλ} \\
           &= \frac{α - 1}{α + 1} \frac{2}{Δλ}

        where R(λ) is the reflectivity, and:

        .. math::

            α = R(λ1) / R(λ0) = 10^{0.4 color_index}

            color_index = Δm - C_{sun}

        Δλ is typically expressed in units of 100 nm.


        Examples
        --------
        >>> import astropy.units as u
        >>> w = [0.4719, 0.6185] * u.um
        >>> S = SpectralGradient.from_color(w, 0.10 * u.mag)
        >>> print(S)                            # doctest: +FLOAT_CMP
        6.27819572 % / (100 nm)

        """
        from ..units import hundred_nm

        lambda_eff = SpectralGradient._lambda_eff(wfb)

        if isinstance(color, u.Magnitude):
            alpha = u.Quantity(-1 * color, u.dimensionless_unscaled)
        elif u.Quantity(color).unit.is_equivalent(u.mag):
            alpha = (-color).to(u.dimensionless_unscaled, u.logarithmic())
        else:
            alpha = u.Quantity(color, u.dimensionless_unscaled)

        dw = lambda_eff[1] - lambda_eff[0]
        S = (2 / dw * (alpha - 1) / (alpha + 1)).to(u.percent / hundred_nm)

        return SpectralGradient(S, wave=lambda_eff)

    def to_color(self, wfb):
        r"""Express as a color index.


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
        Color index is computed from:

        .. math::

            α = \frac{1 + S Δλ / 2}{1 - S * Δλ / 2}

        where S is the spectral gradient at the mean of λ0 and λ1, and:

        .. math::

            α = R(λ1) / R(λ0) = 10^{0.4 color_index}

            color_index = Δm - C_{sun}

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
        >>> print(C)    # doctest: +FLOAT_CMP
        0.05429812423309064 mag

        """

        lambda_eff = self._lambda_eff(wfb)

        S = self.renormalize(lambda_eff.mean())
        dw = lambda_eff[0] - lambda_eff[1]
        beta = (S * dw / 2).decompose()  # dimensionless
        color = ((1 + beta) / (1 - beta)).to(u.mag, u.logarithmic())

        return color

    def renormalize(self, wave0):
        """Re-normalize to another wavelength.

        The slope is linearly extrapolated to the new normalization
        point.  Requires the `wave0` attribute to be defined, see
        `~SpectralGradient`.


        Parameters
        ----------
        wave0 : `~astropy.units.Quantity`
            Wavelength.


        Returns
        -------
        S : ``SpectralGradient``
            Renormalized gradient.


        Examples
        --------
        >>> import astropy.units as u
        >>> from sbpy.spectroscopy import SpectralGradient
        >>> from sbpy.units import hundred_nm
        >>> S1 = SpectralGradient(10 * u.percent / hundred_nm,
        ...                      wave0=0.55 * u.um)
        >>> S2 = S1.renormalize(3.6 * u.um)
        >>> print(S2.renormalize(3.6 * u.um))    # doctest: +FLOAT_CMP
        2.469135802469136 % / (100 nm)
        >>> print(S2.wave0)
        3.6 um
        """

        if self.wave0 is None:
            raise ValueError("wave0 attribute must be defined.")

        delta = wave0 - self.wave0
        S0 = 1 + self.to(delta.unit**-1) * delta
        S = self / S0
        S.wave0 = wave0
        return S


class SpectralModel:
    """Range of spectral models"""

    def haser():
        """Haser model

        should allow direct creation of a `sbpy.actvity.Haser` instance"""
        pass

    def emission_lines():
        """Emission lines"""
        pass

    def reflectance():
        """Reflectance spectrum (asteroids)"""


class Spectrum:

    def __init__(self, flux, dispersionaxis, unit):
        self.flux = flux
        self.axis = dispersionaxis
        self.unit = unit

    @classmethod
    def read(cls, filename, columns="auto"):
        """Read spectrum from file

        Parameters
        ----------
        filename : str, mandatory
            data file name
        columns : str or list-like, optional, default: 'auto'
            file format, `auto` will try to recognize the format
            automatically

        Returns
        -------
        `Spectrum` instance

        Examples
        --------
        >>> spec = Spectrum.read('2012_XY.dat') # doctest: +SKIP

        not yet implemented

        """

    def write(self, filename, columns="all"):
        """Write spectrum to file

        Parameters
        ----------
        filename : str, mandatory
            data file name
        columns : str or list-like, optional: default: 'all'
            file format; `all` will write all fields to the file

        Examples
        --------
        >>> spec = Spectrum.read('2012_XY.dat') # doctest: +SKIP
        >>> spec.write('2012_XY.dat.bak') # doctest: +SKIP

        not yet implemented

        """

    def convert_units(self, **kwargs):
        """Convert Spectrum units as provided by user

        Examples
        --------
        >>> spec.convert_units(flux_unit=u.K) # doctest: +SKIP
        >>> spec.convert_units(dispersion_unit=u.km/u.s) # doctest: +SKIP

        not yet implemented

        """

    def baseline(self, subtract=False):
        """fit baseline to `Spectrum` instance

        Parameters
        ----------
        subtract : bool, optional, default=False
            if `True`, subtract the baseline

        Returns
        -------
        float

        Examples
        --------
        >>> baseline = spec.baseline() # doctest: +SKIP
        >>> spec.baseline(subtract=True) # doctest: +SKIP

        not yet implemented

        """

    def slope(self, subtract=False):
        """fit slope to `Spectrum` instance

        Parameters
        ----------
        subtract : bool, optional, default=False
            if `True`, subtract the slope

        Returns
        -------
        float

        Examples
        --------
        >>> slope = spec.slope() # doctest: +SKIP
        >>> spec.slope(subtract=True) # doctest: +SKIP

        not yet implemented

        """

    def integrated_flux(self, frequency, interval=1 * u.km / u.s):
        """
        Calculate integrated flux of emission line.

        Parameters
        ----------
        frequency : `~astropy.units.Quantity`
            Transition frequency
        interval : `~astropy.units.Quantity`
            line width

        Examples
        --------
        >>> flux = spec.integrated_flux(frequency=556.9*u.GHz, # doctest: +SKIP
        >>>                             interval=1.7*u.km/u.s) # doctest: +SKIP

        not yet implemented

        """

    def fit(self, spec):
        """Fit `SpectralModel` to different model types

        Parameters
        ----------
        spec : str, mandatory
            `SpectralModel` instance to fit

        Examples
        --------
        >>> spec_model = SpectralModel(type='Haser', molecule='H2O')  # doctest: +SKIP

        >>> spec.fit(spec_model) # doctest: +SKIP
        >>> print(spec.fit_info) # doctest: +SKIP

        not yet implemented

        """

    def plot(self):
        """Plot spectrum

        not yet implemented
        """
