# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""activity.gas core

"""

__all__ = [
    'photo_lengthscale',
    'photo_timescale',
    'fluorescence_band_strength',

    'Haser'
]

from warnings import warn
from abc import ABC, abstractmethod

import numpy as np
import astropy.units as u

try:
    import scipy
    from scipy import special
    from scipy.integrate import quad, dblquad
except ImportError:
    scipy = None

from astropy.table import Table
from ... import bib
from ... import data as sbd
from ...exceptions import RequiredPackageUnavailable
from .. core import (Aperture, RectangularAperture, GaussianAperture,
                     AnnularAperture, CircularAperture)


def photo_lengthscale(species, source=None):
    """Photodissociation lengthscale for a gas species.


    Parameters
    ----------
    species : string
        The species to look up.

    source : string, optional
        Retrieve values from this source (case insensitive).  See
        references for keys.


    Returns
    -------
    gamma : `~astropy.units.Quantity`
      The lengthscale at 1 au.


    Examples
    --------
    >>> from sbpy.activity import photo_lengthscale
    >>> gamma = photo_lengthscale('OH')


    References
    ----------
    [CS93] H2O and OH from Table IV of Cochran & Schleicher 1993,
    Icarus 105, 235-253.  Quoted for intermediate solar activity.

    """

    from .data import photo_lengthscale as data

    default_sources = {
        'H2O': 'CS93',
        'OH': 'CS93',
    }

    if species not in data:
        summary = ''
        for k, v in sorted(data.items()):
            summary += '\n{} [{}]'.format(k, ', '.join(v.keys()))

        raise ValueError(
            'Invalid species {}.  Choose from:{}'
            .format(species, summary))

    gas = data[species]
    source = default_sources[species] if source is None else source

    if source not in gas:
        raise ValueError(
            'Source key {} not available for {}.  Choose from: {}'
            .format(source, species, ', '.join(gas.keys())))

    gamma, bibcode = gas[source]
    bib.register(photo_lengthscale, bibcode)

    return gamma


def photo_timescale(species, source=None):
    """Photodissociation timescale for a gas species.


    Parameters
    ----------
    species : string
        Species to look up.

    source : string, optional
        Retrieve values from this source.  See references for keys.


    Returns
    -------
    tau : `~astropy.units.Quantity`
      The timescale at 1 au.  May be a two-element array: (quiet Sun,
      active Sun).


    Examples
    --------
    >>> from sbpy.activity import photo_timescale
    >>> tau = photo_timescale('OH')


    References
    ----------
    [CS93] Table IV of Cochran & Schleicher 1993, Icarus 105, 235-253.
    Quoted for intermediate solar activity.

    [C94] Crovisier 1994, JGR 99, 3777-3781.

    [CE83] Crovisier & Encrenaz 1983, A&A 126, 170-182.

    [H92] Huebner et al. 1992, Astroph. & Space Sci. 195, 1-294.

    """

    from .data import photo_timescale as data

    default_sources = {
        'H2O': 'CS93',
        'OH': 'CS93',
        'HCN': 'C94',
        'CH3OH': 'C94',
        'H2CO': 'C94',
        'CO2': 'CE83',
        'CO': 'CE83',
        'CN': 'H92'
    }

    if species not in data:
        summary = ''
        for k, v in sorted(data.items()):
            summary += '\n{} [{}]'.format(k, ', '.join(v.keys()))

        raise ValueError(
            "Invalid species {}.  Choose from:{}"
            .format(species, summary))

    gas = data[species]
    source = default_sources[species] if source is None else source

    if source not in gas:
        raise ValueError(
            'Source key {} not available for {}.  Choose from: {}'
            .format(source, species, ', '.join(gas.keys())))

    tau, bibcode = gas[source]
    bib.register(photo_timescale, bibcode)

    return tau


@sbd.dataclass_input(eph=sbd.Ephem)
def fluorescence_band_strength(species, eph=None, source=None):
    """Fluorescence band strength.


    Parameters
    ----------
    species : string
        Species to look up.

    eph : `~astropy.units.Quantity`, `~sbpy.data.Ephem` or `dict` optional
        The target ephemeris.  The strength is scaled to the given
        heliocentric distance, if present.  Some species require
        heliocentric radial velocity ('rdot').

    source : string, optional
        Retrieve values from this source (case insensitive).  See
        references for keys.


    Returns
    -------
    LN : `~astropy.units.Quantity`
        Luminosity per molecule, scaled to rh, if provided.


    Examples
    --------
    >>> import astropy.units as u
    >>> from sbpy.activity import fluorescence_band_strength
    >>>
    >>> eph = {'rh': 1 * u.au, 'rdot': -1 * u.km / u.s}
    >>> LN = fluorescence_band_strength('OH 0-0', eph, 'SA88')
    >>> print(LN)    # doctest: +FLOAT_CMP
    [1.54e-15] erg / s

    """

    from .data import fluorescence_band_strength as data

    default_sources = {
        'OH 0-0': 'SA88',
        'OH 1-0': 'SA88',
        'OH 1-1': 'SA88',
        'OH 2-2': 'SA88',
        'OH 0-1': 'SA88',
        'OH 0-2': 'SA88',
        'OH 2-0': 'SA88',
        'OH 2-1': 'SA88',
    }

    if species not in data:
        raise ValueError(
            'No data available for {}.  Choose one of: {}'
            .format(species, ', '.join(data.keys())))

    band = data[species]
    source = default_sources[species] if source is None else source

    if source not in band:
        raise ValueError(
            'No source {} for {}.  Choose one of: {}'
            .format(source, species, ', '.join(band.keys())))

    LN, note, bibcode = band[source]
    if bibcode is not None:
        bib.register(fluorescence_band_strength, bibcode)

    return LN(eph)


class GasComa(ABC):
    """Abstract base class for gas coma models.


    Parameters
    ----------
    Q : `~astropy.units.Quantity`
        Production rate, number per time.

    v : `~astropy.units.Quantity`
        Radial outflow speed, distance per time.

    """

    @u.quantity_input(Q=(u.s**-1, u.mol / u.s), v=u.m / u.s)
    def __init__(self, Q, v):
        self.Q = Q
        self.v = v

    @u.quantity_input(r=u.m)
    def volume_density(self, r):
        """Coma volume density.


        Parameters
        ----------
        r : `~astropy.units.Quantity`
            Linear distance to the nucleus.


        Returns
        -------
        n : `~astropy.units.Quantity`
            Local number density.

        """

        return self._volume_density(r.to('m').value) / u.m**3

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
    def column_density(self, rho, eph=None):
        """Coma column density at a projected distance from nucleus.


        Parameters
        ----------
        rho : `~astropy.units.Quantity`
            Projected distance to the region of interest on the plane
            of the sky in units of length or angle.

        eph : dictionary-like, `~sbpy.data.Ephem`, `~astropy.units.Quantity`. optional
            Target-observer distance, or ephemeris with ``'delta'``
            field.  Required to convert rho to a projected size.


        Returns
        -------
        sigma : `~astropy.units.Quantity`
            Coma column density along the line of sight at a distance
            rho.

        """

        equiv = []
        if eph is not None:
            equiv = sbu.projected_size(eph)

        rho = rho.to('m', equiv).value
        return self._column_density(rho) / u.m**2

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
    def total_number(self, aper, eph=None):
        """Total number of molecules in aperture.


        Parameters
        ----------
        aper : `~astropy.units.Quantity`, `~sbpy.activity.Aperture`
            Observation aperture.  May be a circular aperture radius
            with units of length or angle.

        eph : dictionary-like, `~sbpy.data.Ephem`, `~astropy.units.Quantity`
            Target-observer distance, or ephemeris with `delta`.
            Required if the aperture has angular units.


        Returns
        -------
        N : float
            Total number of molecules within the aperture.

        """

        if eph is not None:
            aper = aper.as_length(eph)
        return self._integrate_column_density(aper)[0]

    @abstractmethod
    def _volume_density(self, r):
        """Unitless volumne density function.


        Parameters
        ----------
        r : float
            Linear distance to the nucleus in meters.


        Returns
        -------
        n : float
            Local number density in inverse cubic-meters.

        """

    @abstractmethod
    def _column_density(self, rho):
        """Unitless column density function.


        Parameters
        ----------
        rho : float
            Projected distance of the region of interest on the plane
            of the sky in units of meters.


        Returns
        -------
        sigma : float
            Coma column density along the line of sight at a distance
            rho in units of inverse square-meters.

        """

    def _integrate_volume_density(self, rho, epsabs=1.49e-8):
        """Integrate volume density along the line of sight.


        Parameters
        ----------
        rho : float
            Projected distance of the region of interest on the plane of
            the sky in units of meters

        epsabs : float, int, optional
            Absolute and relative error tolerance for integrals.  See
            `scipy.integrate.quad`.


        Returns
        -------
        sigma : float
            Coma column density along ``rho`` in units of inverse
            square-meters.

        err : float
            Estimated integration error.

        """

        if not scipy:
            raise RequiredPackageUnavailable('scipy')

        def f(s, rho2):
            r = np.sqrt(rho2 + s**2)
            return self._volume_density(r)

        # quad diverges integrating to infinity, but 1e6 × rho is good
        # enough
        limit = 30
        points = rho * np.logspace(-4, 4, limit / 2)
        sigma, err = quad(f, 0, 1e6 * rho, args=(rho**2,),
                          limit=limit, points=points, epsabs=epsabs)

        # spherical symmetry
        sigma *= 2
        err *= 2

        return sigma, err

    def _integrate_column_density(self, aper, epsabs=1.49e-8):
        """Integrate column density over an aperture.


        Parameters
        ----------
        aper : `~sbpy.activity.Aperture`
            Aperture, in units of length.

        epsabs : float, int, optional
            Absolute and relative error tolerance for integrals.  See
            `scipy.integrate.quad` (circular, annular, Gaussian) and
            `~scipy.integrate.dblquad` (rectangular) for details.


        Returns
        -------
        N : float
            Total number.

        err : float
            Estimated integration error.

        """

        if not scipy:
            raise RequiredPackageUnavailable('scipy')

        if isinstance(aper, (CircularAperture, AnnularAperture)):
            if isinstance(aper, CircularAperture):
                limits = (0, aper.radius.to('m').value)
            else:
                limits = aper.shape.to('m').value

            # integrate in polar coordinates
            def f(rho):
                """Column density integration in polar coordinates.

                rho in m, column_density in m**-2

                """
                return rho * self._column_density(rho)

            N, err = quad(f, *limits, epsabs=epsabs)
            N *= 2 * np.pi
            err *= 2 * np.pi
        elif isinstance(aper, RectangularAperture):
            shape = aper.shape.to('m').value

            def f(rho, th):
                """Column density integration in polar coordinates.

                rho in m, column_density in m**-2

                th is ignored (azimuthal symmetry)

                """
                return rho * self._column_density(rho)

            # first "octant"; rho1 and rho2 are the limits of the
            # integration
            def rho1(th):
                "Lower limit"
                return 0

            def rho2(th):
                "Upper limit (a line)"
                return shape[0] / 2 / np.cos(th)

            th = np.arctan(shape[1] / shape[0])
            N1, err1 = dblquad(f, 0, th, rho1, rho2, epsabs=epsabs)

            # second "octant"
            def rho2(th):
                "Upper limit (a line)"
                return shape[1] / 2 / np.cos(th)

            th = np.arctan(shape[0] / shape[1])
            N2, err2 = dblquad(f, 0, th, rho1, rho2, epsabs=epsabs)

            # N1 + N2 constitute 1/4th of the rectangle
            N = 4 * (N1 + N2)
            err = 4 * (err1 + err2)
        elif isinstance(aper, GaussianAperture):
            # integrate in polar coordinates
            def f(rho, sigma):
                """Column density integration in polar coordinates.

                rho and sigma in m, column_density in m**-2

                """
                return (rho * np.exp(-rho**2 / sigma**2 / 2)
                        * self._column_density(rho))

            sigma = aper.sigma.to('m').value
            N, err = quad(f, 0, np.inf, args=(sigma,), epsabs=epsabs)
            N *= 2 * np.pi
            err *= 2 * np.pi

        return N, err


class Haser(GasComa):
    """Haser coma model.

    Some functions require `scipy`.


    Parameters
    ----------
    Q : `~astropy.units.Quantity`
        Production rate, per time.

    v : `~astropy.units.Quantity`
        Radial outflow speed, distance per time.

    parent : `~astropy.units.Quantity`
        Coma lengthscale of the parent species.

    daughter : `~astropy.units.Quantity`, optional
        Coma lengthscale of the daughter species.


    References
    ----------
    Haser 1957, Bulletin de la Societe Royale des Sciences de Liege
    43, 740.

    Newburn and Johnson 1978, Icarus 35, 360-368.

    """

    @bib.cite({'model': '1957BSRSL..43..740H'})
    @u.quantity_input(parent=u.m, daughter=u.m)
    def __init__(self, Q, v, parent, daughter=None):
        super().__init__(Q, v)
        self.parent = parent
        self.daughter = daughter

    def _volume_density(self, r):
        n = (self.Q / self.v).to('1/m').value / r**2 / 4 / np.pi
        parent = self.parent.to('m').value
        if self.daughter is None or self.daughter == 0:
            # parent only
            n *= np.exp(-r / parent)
        else:
            daughter = self.daughter.to('m').value
            n *= (daughter / (parent - daughter)
                  * (np.exp(-r / parent) - np.exp(-r / daughter)))

        return n

    def _iK0(self, x):
        """Integral of the modified Bessel function of 2nd kind, 0th order."""
        if not scipy:
            raise RequiredPackageUnavailable('scipy')
        return special.iti0k0(x)[1]

    def _K1(self, x):
        """Modified Bessel function of 2nd kind, 1st order."""
        if not scipy:
            raise RequiredPackageUnavailable('scipy')
        return special.k1(x)

    @bib.cite({'model': '1978Icar...35..360N'})
    def _column_density(self, rho):
        sigma = (self.Q / self.v).to('1/m').value / rho / 2 / np.pi
        parent = self.parent.to('m').value
        if self.daughter is None or self.daughter == 0:
            sigma *= np.pi / 2 - self._iK0(rho / parent)
        else:
            daughter = self.daughter.to('m').value
            sigma *= (daughter / (parent - daughter)
                      * (self._iK0(rho / daughter) - self._iK0(rho / parent)))
        return sigma

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
    def total_number(self, aper, eph=None):
        if isinstance(aper, u.Quantity):
            aper = CircularAperture(aper)

        if eph is not None:
            aper = aper.as_length(eph)

        # Inspect aper and handle as appropriate
        if isinstance(aper, (RectangularAperture, GaussianAperture)):
            return super().total_number(aper)
        elif isinstance(aper, AnnularAperture):
            N0 = self.total_number(aper.shape[0])
            N1 = self.total_number(aper.shape[1])
            return N1 - N0

        # Solution for the circular aperture of radius rho:
        bib.register(self.total_number, {'model': '1978Icar...35..360N'})

        rho = aper.radius
        parent = self.parent.to(rho.unit)
        x = (rho / parent).to('').value
        N = (self.Q * rho / self.v).to(u.dimensionless_unscaled).value
        if self.daughter is None or self.daughter == 0:
            N *= 1 / x - self._K1(x) + np.pi / 2 - self._iK0(x)
        else:
            daughter = self.daughter.to(rho.unit)
            y = (rho / daughter).to('').value
            N *= (daughter / (parent - daughter)
                  * (self._iK0(y) - self._iK0(x) + x**-1 - y**-1
                     + self._K1(y) - self._K1(x)))

        return N
    total_number.__doc__ = GasComa.total_number.__doc__
