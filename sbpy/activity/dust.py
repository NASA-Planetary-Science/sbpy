# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy activity.dust.core
=======================

All things dust coma related.


"""

__all__ = ["phase_HalleyMarcus", "DustComaQuantity", "Afrho", "Efrho"]

__doctest_requires__ = {
    "Afrho.to_fluxd": ["synphot"],
    "Efrho.to_fluxd": ["synphot"],
}

import abc

import numpy as np
import astropy.units as u

try:
    from scipy.interpolate import splrep, splev
except ImportError:
    pass

from .. import bib
from ..calib import Sun
from ..spectroscopy import BlackbodySource
from ..utils import optional_packages
from .. import data as sbd
from .. import units as sbu
from ..spectroscopy.sources import SinglePointSpectrumError
from .core import CircularAperture


@bib.cite(
    {
        "Halley-Marcus phase function": "2011AJ....141..177S",
        "Halley phase function": "1998Icar..132..397S",
        "Marcus phase function": "2007ICQ....29...39M",
    }
)
def phase_HalleyMarcus(phase):
    """Halley-Marcus composite dust phase function.

    Uses `~scipy.interpolate` for spline interpolation, otherwise uses
    linear interpolation from `~numpy.interp`.


    Parameters
    ----------
    phase : `~astropy.units.Quantity`, `~astropy.coordinate.Angle`
        Phase angle.


    Returns
    -------
    Phi : float, `~np.ndarray`


    Notes
    -----
    The Halley-Marcus phase function was first used by Schleicher and
    Bair (2011), but only described in detail by Schleicher and Marcus
    (May 2010) online at:

        https://asteroid.lowell.edu/comet/dustphase.html

        "To distinguish this curve from others, we designate this as
        the HM phase function, for the sources of the two components:
        Halley and Marcus, where the Halley curve for smaller phase
        angles comes from our previous work (Schleicher et al. 1998)
        while Joe Marcus has fit a Henyey-Greenstein function to a
        variety of mid- and large-phase angle data sets (Marcus 2007);
        see here for details. Note that we do not consider our
        composite curve to be a definitive result, but rather
        appropriate for performing first-order adjustments to dust
        measurements for changing phase angle."


    References
    ----------
    Schleicher & Bair 2011, AJ 141, 177.
    Schleicher, Millis, & Birch 1998, Icarus 132, 397-417.
    Marcus 2007, International Comets Quarterly 29, 39-66.


    Examples
    --------
    >>> from sbpy.activity import phase_HalleyMarcus
    >>> import astropy.units as u
    >>> phase_HalleyMarcus(0 * u.deg)                 # doctest: +FLOAT_CMP
    1.0
    >>> phase_HalleyMarcus(15 * u.deg)                # doctest: +FLOAT_CMP
    5.8720e-01

    """

    th = np.arange(181)
    ph = np.array(
        [
            1.0000e00,
            9.5960e-01,
            9.2170e-01,
            8.8590e-01,
            8.5220e-01,
            8.2050e-01,
            7.9060e-01,
            7.6240e-01,
            7.3580e-01,
            7.1070e-01,
            6.8710e-01,
            6.6470e-01,
            6.4360e-01,
            6.2370e-01,
            6.0490e-01,
            5.8720e-01,
            5.7040e-01,
            5.5460e-01,
            5.3960e-01,
            5.2550e-01,
            5.1220e-01,
            4.9960e-01,
            4.8770e-01,
            4.7650e-01,
            4.6590e-01,
            4.5590e-01,
            4.4650e-01,
            4.3770e-01,
            4.2930e-01,
            4.2150e-01,
            4.1420e-01,
            4.0730e-01,
            4.0090e-01,
            3.9490e-01,
            3.8930e-01,
            3.8400e-01,
            3.7920e-01,
            3.7470e-01,
            3.7060e-01,
            3.6680e-01,
            3.6340e-01,
            3.6030e-01,
            3.5750e-01,
            3.5400e-01,
            3.5090e-01,
            3.4820e-01,
            3.4580e-01,
            3.4380e-01,
            3.4210e-01,
            3.4070e-01,
            3.3970e-01,
            3.3890e-01,
            3.3850e-01,
            3.3830e-01,
            3.3850e-01,
            3.3890e-01,
            3.3960e-01,
            3.4050e-01,
            3.4180e-01,
            3.4320e-01,
            3.4500e-01,
            3.4700e-01,
            3.4930e-01,
            3.5180e-01,
            3.5460e-01,
            3.5760e-01,
            3.6090e-01,
            3.6450e-01,
            3.6830e-01,
            3.7240e-01,
            3.7680e-01,
            3.8150e-01,
            3.8650e-01,
            3.9170e-01,
            3.9730e-01,
            4.0320e-01,
            4.0940e-01,
            4.1590e-01,
            4.2280e-01,
            4.3000e-01,
            4.3760e-01,
            4.4560e-01,
            4.5400e-01,
            4.6270e-01,
            4.7200e-01,
            4.8160e-01,
            4.9180e-01,
            5.0240e-01,
            5.1360e-01,
            5.2530e-01,
            5.3750e-01,
            5.5040e-01,
            5.6380e-01,
            5.7800e-01,
            5.9280e-01,
            6.0840e-01,
            6.2470e-01,
            6.4190e-01,
            6.5990e-01,
            6.7880e-01,
            6.9870e-01,
            7.1960e-01,
            7.4160e-01,
            7.6480e-01,
            7.8920e-01,
            8.1490e-01,
            8.4200e-01,
            8.7060e-01,
            9.0080e-01,
            9.3270e-01,
            9.6640e-01,
            1.0021e00,
            1.0399e00,
            1.0799e00,
            1.1223e00,
            1.1673e00,
            1.2151e00,
            1.2659e00,
            1.3200e00,
            1.3776e00,
            1.4389e00,
            1.5045e00,
            1.5744e00,
            1.6493e00,
            1.7294e00,
            1.8153e00,
            1.9075e00,
            2.0066e00,
            2.1132e00,
            2.2281e00,
            2.3521e00,
            2.4861e00,
            2.6312e00,
            2.7884e00,
            2.9592e00,
            3.1450e00,
            3.3474e00,
            3.5685e00,
            3.8104e00,
            4.0755e00,
            4.3669e00,
            4.6877e00,
            5.0418e00,
            5.4336e00,
            5.8682e00,
            6.3518e00,
            6.8912e00,
            7.4948e00,
            8.1724e00,
            8.9355e00,
            9.7981e00,
            1.0777e01,
            1.1891e01,
            1.3166e01,
            1.4631e01,
            1.6322e01,
            1.8283e01,
            2.0570e01,
            2.3252e01,
            2.6418e01,
            3.0177e01,
            3.4672e01,
            4.0086e01,
            4.6659e01,
            5.4704e01,
            6.4637e01,
            7.7015e01,
            9.2587e01,
            1.1237e02,
            1.3775e02,
            1.7060e02,
            2.1348e02,
            2.6973e02,
            3.4359e02,
            4.3989e02,
            5.6292e02,
            7.1363e02,
            8.8448e02,
            1.0533e03,
            1.1822e03,
            1.2312e03,
        ]
    )

    _phase = np.abs(u.Quantity(phase, "deg").value)

    if optional_packages("scipy"):
        Phi = splev(_phase, splrep(th, ph))
    else:
        Phi = np.interp(_phase, th, ph)

    if np.iterable(_phase):
        Phi = np.array(Phi).reshape(np.shape(_phase))
    else:
        Phi = float(Phi)

    return Phi


class DustComaQuantity(u.SpecificTypeQuantity, metaclass=abc.ABCMeta):
    """Abstract base class for dust coma photometric models: Afrho, Efrho."""

    _equivalent_unit = u.meter
    _include_easy_conversion_members = False

    def __new__(cls, value, unit=None, dtype=None, copy=True):
        return super().__new__(cls, value, unit=unit, dtype=dtype, copy=copy)

    @classmethod
    def from_fluxd(cls, wfb, fluxd, aper, eph, **kwargs):
        """Initialize from spectral flux density.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, list
            Wavelengths, frequencies, bandpass, or list of
            bandpasses of the observation.  Bandpasses require
            `~synphot`.

        fluxd : `~astropy.units.Quantity`
            Flux density per unit wavelength or frequency.

        aper : `~astropy.units.Quantity` or `~sbpy.activity.Aperture`
            Aperture of the observation as a circular radius (length
            or angular units), or as an `~sbpy.activity.Aperture`.

        eph: dictionary-like, `~sbpy.data.Ephem`
            Ephemerides of the comet.  Required fields: 'rh', 'delta'.
            Optional: 'phase'.

        **kwargs
            Keyword arguments for `~to_fluxd`.

        """

        fluxd1cm = cls(1 * u.cm).to_fluxd(
            wfb,
            aper,
            eph,
            unit=fluxd.unit,
            **kwargs,
        )

        if isinstance(fluxd1cm, u.Magnitude):
            coma = cls((fluxd - fluxd1cm).physical * u.cm)
        else:
            coma = cls((fluxd / fluxd1cm).decompose() * u.cm)

        return coma

    @sbd.dataclass_input(eph=sbd.Ephem)
    def to_fluxd(self, wfb, aper, eph, unit=None, **kwargs):
        """Express as spectral flux density in an observation.

        Assumes the small angle approximation.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`, list
            Wavelengths, frequencies, bandpass, or list of
            bandpasses of the observation.  Bandpasses require
            `~synphot`.  Ignored if ``S`` is provided.

        aper: `~astropy.units.Quantity`, `~sbpy.activity.Aperture`
            Aperture of the observation as a circular radius (length
            or angular units), or as an sbpy `~sbpy.activity.Aperture`.

        eph: dictionary-like, `~sbpy.data.Ephem`
            Ephemerides of the comet.  Required fields: 'rh', 'delta'.
            Optional: 'phase'.

        unit : `~astropy.units.Unit`, string, optional
            The flux density unit for the output.

        """

        # This method handles the geometric quantities.  Sub-classes
        # will handle the photometric quantities in `_source_fluxd`.

        # rho = effective circular aperture radius at the distance of
        # the comet.  Keep track of array dimensionality as Ephem
        # objects can needlessly increase the number of dimensions.
        rho = CircularAperture.from_coma_equivalent(aper).dim
        ndim = np.ndim(rho)
        rho = rho.to("km", sbu.projected_size(eph))

        ndim = max(ndim, np.ndim(self))

        # validate unit
        if unit is not None:
            unit = u.Unit(unit)

        # get source spectral flux density
        # * sunlight for Afrho,
        # * blackbody emission for Efrho
        # quantity = (delta**2 * F / rho) / source
        # must have spectral flux density units
        source = self._source_fluxd(wfb, eph, unit=unit, **kwargs)

        if isinstance(source, u.Magnitude):
            _source = source.physical
        else:
            _source = source
        fluxd = self * rho / eph["delta"] ** 2 * _source

        # using Ephem can unnecessarily promote fluxd to an array
        if np.ndim(fluxd) > ndim:
            fluxd = np.squeeze(fluxd)

        # and back to magnitudes, as needed
        return fluxd.to(source.unit)

    @abc.abstractmethod
    def _source_fluxd(self, wfb, eph, unit=None, **kwargs):
        """Photometric calibration of dust coma quantity.

        quantity = delta**2 * F / rho / source

        delta - observer-comet distance
        F - observed spectral flux density
        rho - photometric aperture radius at the distance of the comet
        source - source function flux density (this method)

        For Afrho, source = S / rh**2 / 4 * Phi(phase).

        For Efrho, source = 1 / pi / B(T).

        Must respect requested units.

        """


class Afrho(DustComaQuantity):
    """Coma dust quantity for scattered light.

    ``Afrho`` objects behave like `~astropy.units.Quantity` objects
    with units of length.


    Parameters
    ----------
    value : number, `~astropy.units.Quantity`
        The value(s).

    unit : string, `~astropy.units.Unit`, optional
        The unit of the input value.  Strings must be parseable by
        :mod:`~astropy.units` package.

    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`.

    copy : bool, optional
        See `~astropy.units.Quantity`.


    Notes
    -----
    Afρ is the product of dust albedo, dust filling factor, and
    circular aperture radius.  It is nominally a constant for a
    steady-state coma in free expansion.  See A'Hearn et al. (1984)
    for details.


    References
    ----------
    A'Hearn et al. 1984, AJ 89, 579-591.


    Examples
    --------
    >>> from sbpy.activity import Afrho
    >>> print(Afrho(1000.0, 'cm'))
    1000.0 cm

    """

    @classmethod
    def from_fluxd(cls, wfb, fluxd, aper, eph, **kwargs):
        return super().from_fluxd(wfb, fluxd, aper, eph, **kwargs)

    from_fluxd.__doc__ = (
        DustComaQuantity.from_fluxd.__doc__
        + """
        Examples
        --------
        Convert observed flux density to Afρ, with a user-provided
        solar flux density for the V-band:
        >>> from sbpy.activity import Afrho
        >>> import astropy.units as u
        >>> from sbpy.calib import solar_fluxd
        >>>
        >>> solar_fluxd.set({'V': 1869 * u.W / u.m**2 / u.um})
        >>>
        >>> fluxd = 6.730018324465526e-14 * u.W / u.m**2 / u.um
        >>> aper = 1 * u.arcsec
        >>> eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        >>> afrho = Afrho.from_fluxd('V', fluxd, aper, eph=eph)
        >>> print(afrho)                        # doctest: +FLOAT_CMP
        999.9999999999999 cm

        """
    )

    def to_fluxd(self, wfb, aper, eph, unit=None, phasecor=False, Phi=None):
        return super().to_fluxd(wfb, aper, eph, unit=unit, phasecor=phasecor, Phi=Phi)

    to_fluxd.__doc__ = (
        DustComaQuantity.to_fluxd.__doc__
        + """
        phasecor: bool, optional
            Scale the result by the phase function ``Phi``, assuming
            ``Afrho`` is quoted for 0° phase.  Requires phase angle in
            ``eph``.

        Phi : callable, optional
            Phase function, see :func:`~Afrho.to_phase`.

        **kwargs
            Keyword arguments for `~Sun.observe`.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`
            Spectral flux density.


        Examples
        --------
        >>> from sbpy.activity import Afrho
        >>> import astropy.units as u
        >>> afrho = Afrho(1000 * u.cm)
        >>> wave = 0.55 * u.um
        >>> aper = 1 * u.arcsec
        >>> eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        >>> fluxd = afrho.to_fluxd(wave, aper, eph)
        >>> print(fluxd)    # doctest: +FLOAT_CMP
        6.730018324465526e-14 W / (um m2)

        With a phase correction:
        >>> eph['phase'] = 30 * u.deg
        >>> fluxd = afrho.to_fluxd(wave, aper, eph, phasecor=True)
        >>> print(fluxd)    # doctest: +FLOAT_CMP
        2.8017202649540757e-14 W / (um m2)

        In magnitudes through the Johnson V filter:
        >>> import sbpy.units as sbu
        >>> from sbpy.photometry import bandpass
        >>> bp = bandpass('Johnson V')
        >>> fluxd = afrho.to_fluxd(bp, aper, eph, unit=sbu.JMmag,
        ...                     phasecor=True)
        >>> print(fluxd)    # doctest: +FLOAT_CMP
        15.321242371548918 mag(JM)

        """
    )

    @bib.cite({"model": "1984AJ.....89..579A"})
    def _source_fluxd(self, wfb, eph, unit=None, phasecor=False, Phi=None, **kwargs):
        # get solar flux density
        sun = Sun.from_default()
        try:
            S = sun.observe(wfb, unit=unit, **kwargs)
        except SinglePointSpectrumError:
            S = sun(wfb, unit=unit)

        if not (
            S.unit.is_equivalent(u.W / u.m**2 / u.um)
            or S.unit.is_equivalent(u.W / u.m**2 / u.Hz)
            or isinstance(S, u.Magnitude)
        ):
            raise ValueError(
                "Solar flux density must have units of spectral flux "
                "density, e.g., W/m2/μm or W/m2/Hz"
            )

        if phasecor:
            Phi = phase_HalleyMarcus if Phi is None else Phi
            _Phi = Phi(eph["phase"]) / Phi(0 * u.deg)
        else:
            _Phi = 1

        # compute
        _S = S.physical if isinstance(S, u.Magnitude) else S
        source = _S * _Phi / 4 * u.au**2 / eph["rh"] ** 2

        return source.to(S.unit)

    def to_phase(self, to_phase, from_phase, Phi=None):
        """Scale to another phase angle.


        Parameters
        ----------
        to_phase : `~astropy.units.Quantity`
            New target phase angle.

        from_phase : `~astropy.units.Quantity`
            Current target phase angle.

        Phi : callable, optional
            Phase function, a callable object that takes a single
            parameter, phase angle as a `~astropy.units.Quantity`, and
            returns a scale factor.  Default is `~phase_HalleyMarcus`.


        Returns
        -------
        afrho : `~Afrho`
            The scaled Afρ quantity.


        Examples
        --------
        >>> from sbpy.activity import Afrho
        >>> afrho = Afrho(10 * u.cm).to_phase(15 * u.deg, 0 * u.deg)
        >>> print(afrho)                       # doctest: +FLOAT_CMP
        5.87201 cm

        """

        if Phi is None:
            Phi = phase_HalleyMarcus

        return self * Phi(to_phase) / Phi(from_phase)


class Efrho(DustComaQuantity):
    """Coma dust quantity for thermal emission.

    ``Efrho`` behave like `~astropy.units.Quantity` objects with units
    of length.


    Parameters
    ----------
    value : number, `~astropy.units.Quantity`
        The value(s).

    unit : str, `~astropy.units.Unit`, optional
        The unit of the input value.  Strings must be parseable by
        :mod:`~astropy.units` package.

    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`.

    copy : bool, optional
        See `~astropy.units.Quantity`.


    Notes
    -----
    εfρ is the product of dust emissivity, dust filling factor, and
    circular aperture radius.  It is nominally a constant for a
    steady-state coma in free expansion, and is the thermal emission
    equivalent for the Afρ quantity.  See Kelley et al. (2013) for
    details.


    References
    ----------
    A'Hearn et al. 1984, AJ 89, 579-591.
    Kelley et al. 2013, Icarus 225, 475-494.


    Examples
    --------
    >>> from sbpy.activity import Efrho
    >>> print(Efrho(1000.0, 'cm'))
    1000.0 cm

    """

    @classmethod
    def from_fluxd(cls, wfb, fluxd, aper, eph, **kwargs):
        return super().from_fluxd(wfb, fluxd, aper, eph, **kwargs)

    from_fluxd.__doc__ = (
        DustComaQuantity.from_fluxd.__doc__
        + """
        Examples
        --------
        >>> from sbpy.activity import Efrho
        >>> import astropy.units as u
        >>> wave = 15.8 * u.um
        >>> fluxd = 6.52 * u.mJy
        >>> aper =  11.1 * u.arcsec
        >>> eph = {'rh': 4.42 * u.au, 'delta': 4.01 * u.au}
        >>> efrho = Efrho.from_fluxd(wave, fluxd, aper, eph=eph)
        >>> print(efrho)                          # doctest: +FLOAT_CMP
        120.00836963059808 cm

    """
    )

    def to_fluxd(self, wfb, aper, eph, unit=None, Tscale=1.1, T=None, B=None):
        return super().to_fluxd(wfb, aper, eph, unit=unit, Tscale=Tscale, T=T, B=B)

    to_fluxd.__doc__ = (
        DustComaQuantity.to_fluxd.__doc__
        + """
        Tscale : float, optional
            Scale factor for blackbody in LTE with sunlight.  Ignored
            if ``T`` or ``B`` is provided.

        T : `~astropy.units.Quantity`, optional
            Blackbody temperature.  Ignored if ``B`` is provided.

        B : `~astropy.units.Quantity`, optional
            Observed spectral flux density from a blackbody sphere,
            i.e., pi * Planck function.  Overrides ``T`` and
            ``Tscale``.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`
            Spectral flux density.


        Examples
        --------
        >>> from sbpy.activity import Efrho
        >>> import astropy.units as u
        >>> efrho = Efrho(120.0, 'cm')
        >>> freq = 15.8 * u.um
        >>> aper = 11.1 * u.arcsec
        >>> eph = {'rh': 4.42 * u.au, 'delta': 4.01 * u.au}
        >>> fluxd = efrho.to_fluxd(freq, aper, eph=eph, unit='Jy')
        >>> print(fluxd)                           # doctest: +FLOAT_CMP
        0.006519545281786034 Jy

    """
    )

    @bib.cite({"model": "2013Icar..225..475K"})
    def _source_fluxd(self, wfb, eph, unit=None, Tscale=1.1, T=None, B=None):
        if T is None:
            T = Tscale * 278 / np.sqrt(eph["rh"].to("au").value)

        if B is None:
            BB = BlackbodySource(T)
            try:
                B = BB.observe(wfb, unit=unit)
            except SinglePointSpectrumError:
                B = BB(wfb, unit=unit)
        else:
            if not (
                B.unit.is_equivalent(u.W / u.m**2 / u.um)
                or B.unit.is_equivalent(u.W / u.m**2 / u.Hz)
                or isinstance(B, u.Magnitude)
            ):
                raise ValueError(
                    "B must be a magnitude or have units of spectral "
                    "flux density, e.g., W/m2/μm or W/m2/Hz"
                )

        return B
