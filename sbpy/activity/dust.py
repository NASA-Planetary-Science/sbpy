# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
==========================
SBPy Activity: Dust Module
==========================

All things dust coma related.

Functions
---------
phase_HalleyMarcus - Halley-Marcus composite dust phase function.

Classes
-------
Afrho     - Coma dust quantity of A'Hearn et al. (1984).
Syndynes  - Dust dynamical model for zero-ejection velocities.


"""

import numpy as np
import astropy.units as u


__all__ = [
    'phase_HalleyMarcus',
    'Afrho',
    'Syndynes'
]

def phase_HalleyMarcus(phase):
    """Halley-Marcus composite dust phase function.

    Uses `~scipy.interpolate` for spline interpolation, otherwise uses
    linear interpolation from `~numpy.interp`.


    Parameters
    ----------
    phase : `~astropy.units.Quantity`
      Phase angle.


    Returns
    -------
    Phi : float or `~np.ndarray`


    Notes
    -----
    The Halley phase function was first used by Schleicher and Bair
    (2011), but only described in detail by Schleicher and Marcus (May
    2010) online at:

      http://asteroid.lowell.edu/comet/dustphase.html

      "To distinguish this curve from others, we designate this as the
      HM phase function, for the sources of the two components: Halley
      and Marcus, where the Halley curve for smaller phase angles
      comes from our previous work (Schleicher et al. 1998) while Joe
      Marcus has fit a Henyey-Greenstein function to a variety of mid-
      and large-phase angle data sets (Marcus 2007); see here for
      details. Note that we do not consider our composite curve to be
      a definitive result, but rather appropriate for performing
      first-order adjustments to dust measurements for changing phase
      angle."


    References
    ----------
    Schleicher & Bair 2011, AJ 141, 177.
    Schleicher, Millis, & Birch 1998, Icarus 132, 397-417.
    Marcus 2007, International Comets Quarterly 29, 39-66.


    Examples
    --------
    >>> from sbpy.activity import phase_HalleyMarcus
    >>> import astropy.units as u
    >>> phase_HalleyMarcus(0 * u.deg)                    # doctest: +FLOAT_CMP
    1.0
    >>> phase_HalleyMarcus(15 * u.deg)                   # doctest: +FLOAT_CMP
    5.8720e-01

    """

    from .. import bib
    bib.register('activity.dust.phase_HalleyMarcus',
                 ('1998Icar..132..397S', '2007ICQ....29...39M'))

    th = np.arange(181)
    ph = np.array(
        [  1.0000e+00,   9.5960e-01,   9.2170e-01,   8.8590e-01,
           8.5220e-01,   8.2050e-01,   7.9060e-01,   7.6240e-01,
           7.3580e-01,   7.1070e-01,   6.8710e-01,   6.6470e-01,
           6.4360e-01,   6.2370e-01,   6.0490e-01,   5.8720e-01,
           5.7040e-01,   5.5460e-01,   5.3960e-01,   5.2550e-01,
           5.1220e-01,   4.9960e-01,   4.8770e-01,   4.7650e-01,
           4.6590e-01,   4.5590e-01,   4.4650e-01,   4.3770e-01,
           4.2930e-01,   4.2150e-01,   4.1420e-01,   4.0730e-01,
           4.0090e-01,   3.9490e-01,   3.8930e-01,   3.8400e-01,
           3.7920e-01,   3.7470e-01,   3.7060e-01,   3.6680e-01,
           3.6340e-01,   3.6030e-01,   3.5750e-01,   3.5400e-01,
           3.5090e-01,   3.4820e-01,   3.4580e-01,   3.4380e-01,
           3.4210e-01,   3.4070e-01,   3.3970e-01,   3.3890e-01,
           3.3850e-01,   3.3830e-01,   3.3850e-01,   3.3890e-01,
           3.3960e-01,   3.4050e-01,   3.4180e-01,   3.4320e-01,
           3.4500e-01,   3.4700e-01,   3.4930e-01,   3.5180e-01,
           3.5460e-01,   3.5760e-01,   3.6090e-01,   3.6450e-01,
           3.6830e-01,   3.7240e-01,   3.7680e-01,   3.8150e-01,
           3.8650e-01,   3.9170e-01,   3.9730e-01,   4.0320e-01,
           4.0940e-01,   4.1590e-01,   4.2280e-01,   4.3000e-01,
           4.3760e-01,   4.4560e-01,   4.5400e-01,   4.6270e-01,
           4.7200e-01,   4.8160e-01,   4.9180e-01,   5.0240e-01,
           5.1360e-01,   5.2530e-01,   5.3750e-01,   5.5040e-01,
           5.6380e-01,   5.7800e-01,   5.9280e-01,   6.0840e-01,
           6.2470e-01,   6.4190e-01,   6.5990e-01,   6.7880e-01,
           6.9870e-01,   7.1960e-01,   7.4160e-01,   7.6480e-01,
           7.8920e-01,   8.1490e-01,   8.4200e-01,   8.7060e-01,
           9.0080e-01,   9.3270e-01,   9.6640e-01,   1.0021e+00,
           1.0399e+00,   1.0799e+00,   1.1223e+00,   1.1673e+00,
           1.2151e+00,   1.2659e+00,   1.3200e+00,   1.3776e+00,
           1.4389e+00,   1.5045e+00,   1.5744e+00,   1.6493e+00,
           1.7294e+00,   1.8153e+00,   1.9075e+00,   2.0066e+00,
           2.1132e+00,   2.2281e+00,   2.3521e+00,   2.4861e+00,
           2.6312e+00,   2.7884e+00,   2.9592e+00,   3.1450e+00,
           3.3474e+00,   3.5685e+00,   3.8104e+00,   4.0755e+00,
           4.3669e+00,   4.6877e+00,   5.0418e+00,   5.4336e+00,
           5.8682e+00,   6.3518e+00,   6.8912e+00,   7.4948e+00,
           8.1724e+00,   8.9355e+00,   9.7981e+00,   1.0777e+01,
           1.1891e+01,   1.3166e+01,   1.4631e+01,   1.6322e+01,
           1.8283e+01,   2.0570e+01,   2.3252e+01,   2.6418e+01,
           3.0177e+01,   3.4672e+01,   4.0086e+01,   4.6659e+01,
           5.4704e+01,   6.4637e+01,   7.7015e+01,   9.2587e+01,
           1.1237e+02,   1.3775e+02,   1.7060e+02,   2.1348e+02,
           2.6973e+02,   3.4359e+02,   4.3989e+02,   5.6292e+02,
           7.1363e+02,   8.8448e+02,   1.0533e+03,   1.1822e+03,
           1.2312e+03])

    try:
        from scipy.interpolate import splrep, splev
        Phi = splev(np.abs(phase), splrep(th, ph))
    except ImportError as e:
        from astropy.utils.exceptions import AstropyWarning
        from warnings import warn
        warn(AstropyWarning('scipy is not present, using linear interpolation.'))
        Phi = np.interp(np.abs(phase), th, ph)

    if np.iterable(phase):
        Phi = np.array(Phi).reshape(phase)
    else:
        Phi = float(Phi)

    return Phi


class Afrho(u.SpecificTypeQuantity):
    """Coma dust quantity of A'Hearn et al. (1984).

    `Afrho` objects behave like astropy `~astropy.units.Quantity`
    objects with units of length.


    Parameters
    ----------
    value : number, astropy `~astropy.units.Quantity`
      The value(s).

    unit : string, astropy `~Unit`
      The unit of the input value, if `value` is a number.  If a
      string, it must be parseable by :mod:`~astropy.units` package.

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
    >>> import astropy.units as u
    >>> Afrho(1000 * u.cm)
    <Afrho 1000.0 cm>

    """

    _equivalent_unit = u.meter
    _include_easy_conversion_members = True
    
    def __new__(cls, value, unit=None, dtype=None, copy=None):
        return super().__new__(cls, value, unit=unit, dtype=dtype, copy=copy)

    @classmethod
    def from_fluxd(cls, wave_or_freq, fluxd, aper, eph=None, S=None, **kwargs):
        """Initialize from flux density.

        Phase angle is not considered.

        Parameters
        ----------
        wave_or_freq : `~astropy.units.Quantity`
          Wavelengths or frequencies of the observation.

        fluxd : `~astropy.units.Quantity`
          Flux density per unit wavelength or frequency, depending on
          the units of `wave_or_freq`.

        aper : `~astropy.units.Quantity`, `~Aperture`
          Aperture of the observation as a circular radius (length
          or angular units), or as an sbpy `Aperture` class.

        eph : dictionary-like or `~Ephem`, optional
          Ephemerides of the comet, describing heliocentric and
          geocentric distances as `~astropy.units.Quantity` via
          keywords `rh` and `delta`.

        S : `~astropy.units.Quantity`, optional
          Solar flux density per unit wavelength at 1 au and `wave`.
          If `None`, then it will be retrieved from
          `~sbpy.data.solar_fluxd`.

        **kwargs
          Additional keyword arguments are passed to
          `~sbpy.data.solar_fluxd`.

        
        Example
        -------
        >>> from sbpy.activity import Afrho
        >>> import astropy.units as u
        >>> fluxd = 6.730018324465526e-14 * u.W / u.m**2 / u.um
        >>> aper = 1 * u.arcsec
        >>> eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        >>> S = 1869 * u.W / u.m**2 / u.um
        >>> afrho = Afrho.from_fluxd(None, fluxd, aper, eph=eph, S=S)
        >>> afrho.cm                              # doctest: +FLOAT_CMP
        999.9999999999999

        """

        fluxd1cm = Afrho(1 * u.cm).fluxd(wave_or_freq, aper, eph=eph,
                                         S=S, **kwargs)
        return Afrho((fluxd / fluxd1cm).decompose() * u.cm)

    def from_mag(self, filt, mag, aper, eph=None, S=None, **kwargs):
        """Initialize from apparent magnitude.

        Phase angle is not considered.

        Parameters
        ----------
        filt : string
          Name of the filter/bandpass.

        mag : `~astropy.units.Quantity` ???
          Apparent magntiude.

        aper : `~astropy.units.Quantity`, `~Aperture`
          Aperture of the observation as a circular radius (length
          or angular units), or as an sbpy `Aperture` class.

        eph : dictionary-like or `~Ephem`, optional
          Ephemerides of the comet, describing heliocentric and
          geocentric distances as `~astropy.units.Quantity` via
          keywords `rh` and `delta`.

        S : `~astropy.units.Quantity`, optional
          Absolute magnitude of the Sun at `filter`.  If `None`, then
          it will be retrieved from `~sbpy.data.solar_absmag`.

        **kwargs
          Additional keyword arguments are passed to
          `~sbpy.data.solar_absmag`.

        
        Example
        -------

        """
        
        raise NotImplemented


    def fluxd(self, wave_or_freq, aper, eph=None, S=None, **kwargs):
        """Coma flux density.

        Phase angle is not considered.

        Parameters
        ----------
        wave_or_freq : `~astropy.units.Quantity`
          Wavelengths or frequencies of the observation, or `None` if
          `S` is provided.

        aper : `~astropy.units.Quantity`, `~Aperture`
          Aperture of the observation as a circular radius (length
          or angular units), or as an sbpy `~Aperture` class.

        eph : dictionary-like or `~Ephem`, optional
          Ephemerides of the comet, describing heliocentric and
          geocentric distances as `~astropy.units.Quantity` via
          keywords `rh` and `delta`.

        S : `~astropy.units.Quantity`, optional 
          Solar flux density at 1 au.  If `None`, then it will be
          retrieved from `~sbpy.data.solar_fluxd`.

        **kwargs
          Additional keyword arguments are passed to
          `~sbpy.data.solar_fluxd`.

        Returns
        -------
        fluxd : `~astropy.units.Quantity`
          Flux density per unit wavelength or frequency, depending on
          the units of `wave_or_freq` or `S`.

        Example
        -------
        >>> from sbpy.activity import Afrho
        >>> import astropy.units as u
        >>> afrho = Afrho(1000, 'cm')
        >>> aper = 1 * u.arcsec
        >>> eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        >>> S = 1869 * u.W / u.m**2 / u.um
        >>> fluxd = afrho.fluxd(None, aper, eph=eph, S=S)
        >>> fluxd.value                                  # doctest: +FLOAT_CMP
        6.730018324465526e-14

        """

        from . import core
        #from ..data import solar_fluxd
        from .. import bib

        bib.register('activity.dust.Afrho', '1984AJ.....89..579A')

        # check aperture radius
        if isinstance(aper, core.Aperture):
            rho = aper.coma_equivalent_radius()
        else:
            rho = aper

        rho = core.rho_as_distance(rho, eph)

        # check solar flux density
        if S is None:
            S = solar_fluxd(wave_or_freq)
        else:
            assert (S.unit.is_equivalent(u.W / u.m**2 / u.um)
                    or S.unit.is_equivalent(u.W / u.m**2 / u.Hz))

        # compute
        fluxd = self * rho * S / 4 / eph['delta']**2 * u.au**2 / eph['rh']**2

        return fluxd.to(S.unit)

    def mag(self, filt, aper, eph=None, S=None, **kwargs):
        """Coma apparent magnitude.

        Phase angle is not considered.

        Parameters
        ----------
        filt : string
          The name of the filter / bandpass.

        aper : `~astropy.units.Quantity`, `~Aperture`
          Aperture of the observation as a circular radius (length
          or angular units), or as an sbpy `~Aperture` class.

        eph : dictionary-like or `~Ephem`, optional
          Ephemerides of the comet, describing heliocentric and
          geocentric distances as `~astropy.units.Quantity` via
          keywords `rh` and `delta`.

        S : `~astropy.units.Quantity`, optional 
          Absolute magnitude of the Sun at `filt`.  If `None`, then it
          will be retrieved from `~sbpy.data.solar_absmag`.

        **kwargs
          Additional keyword arguments are passed to
          `~sbpy.data.solar_absmag`.

        Returns
        -------
        mag : `~astropy.units.Quantity` ???


        Example
        -------

        """
        raise NotImplemented

    def to_phase(self, to_phase, from_phase, Phi=None):
        """Scale Afρ to another phase angle.

        Parameters
        ----------
        to_phase : `~astropy.units.Quantity`
          New target phase angle.

        from_phase : `~astropy.units.Quantity`
          Current target phase angle.

        Phi : callable or `None`
          Phase function, a callable object that takes a single
          parameter, phase angle as a `~astropy.units.Quantity`, and
          returns a scale factor.  If `None`, `phase_HalleyMarcus` is
          used.  The phase function is expected to be 1.0 at 0 deg.


        Returns
        -------
        afrho : `~Afrho`
          The scaled Afρ quantity.

        
        Examples
        --------
        >>> from sbpy.activity import Afrho
        >>> afrho = Afrho(10 * u.cm).to_phase(15 * u.deg, 0 * u.deg)
        >>> afrho.cm                                     # doctest: +FLOAT_CMP
        5.87201
        
        """

        if Phi is None:
            Phi = phase_HalleyMarcus
        
        return self * Phi(to_phase) / Phi(from_phase)
        

class Syndynes:
    """Syndynes and Synchrones"""

    def __init__(self, orb, date):
        self.orb = orb
        self.date = date
        raise NotImplemented

    def plot_syndynes(self, betas, ages, location='500'):
        """Parameters
        ----------
        betas : array, mandatory
            beta values
        ages : array, mandatory
            synchrone ages
        location : str, optional, default: '500' (geocentric)
            observer location MPC code

        Returns
        -------
        `matplotlib.pyplot` instance

        Examples
        --------
        TBD

        not yet implemented

        """
        
        pass
