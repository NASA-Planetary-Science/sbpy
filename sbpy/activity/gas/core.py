# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
==================
SBPy Activity: Gas
==================


Functions
---------
photo_lengthscale          - Photodissociation lengthscale.
photo_timescale            - Photodissociation timescale.
fluorescence_band_strength - Fluorescence band efficiency of a specific
                             species and transition.


Classes
-------
GasComa             - Abstract base class for gas coma models.
Haser               - Haser coma model for gas (Haser 1957).


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
from astropy.utils.exceptions import AstropyWarning
from ... import bib
from ... import data as sbd
from .. core import (Aperture, RectangularAperture, GaussianAperture,
                     AnnularAperture, CircularAperture)
from .. core import rho_as_length


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

    data = {   # (value, {key feature: ADS bibcode})
        'H2O': {
            'CS93': (2.4e4 * u.km,
                     {'H2O photodissociation lengthscale':
                      '1993Icar..105..235C'})
        },
        'OH': {
            'CS93': (1.6e5 * u.km,
                     {'OH photodissociation lengthscale':
                      '1993Icar..105..235C'})
        },
    }

    default_sources = {
        'H2O': 'CS93',
        'OH': 'CS93',
    }

    if species.upper() not in data:
        raise ValueError(
            'No timescale available for {}.  Choose from: {}'
            .format(species, ', '.join(data.keys())))

    gas = data[species.upper()]
    source = default_sources[species.upper()] if source is None else source

    if source.upper() not in gas:
        raise ValueError(
            'Source key {} not available for {}.  Choose from: {}'
            .format(source, species, ', '.join(gas.keys())))

    gamma, bibcode = gas[source.upper()]
    bib.register(photo_lengthscale, bibcode)

    return gamma


def photo_timescale(species, source=None):
    """Photodissociation timescale for a gas species.


    Parameters
    ----------
    species : string, ``None``
      The species to look up, or ``None`` to summarize available
      species.

    source : string, optional
      Retrieve values from this source (case insensitive).  See
      references for keys.


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

    data = {   # (value, {key feature: ADS bibcode})
        'H2O': {
            'CS93': (5.2e4 * u.s,
                     {'H2O photodissociation timescale':
                      '1993Icar..105..235C'})
        },
        'OH': {
            'CS93': (1.6e5 * u.s,
                     {'OH photodissociation timescale':
                      '1993Icar..105..235C'})
        },
        'HCN': {
            'C94': (6.7e4 * u.s,
                    {'HCN photodissociation timescale':
                     '1994JGR....99.3777C'})
        },
        'CH3OH': {
            'C94': (7.7e4 * u.s,
                    {'CH3OH photodissociation timescale':
                     '1994JGR....99.3777C'})
        },
        'H2CO': {
            'C94': (5.0e3 * u.s,
                    {'H2CO photodissociation timescale':
                     '1994JGR....99.3777C'})
        },
        'CO': {
            'CE83': (1.5e6 * u.s,
                     {'CO photodissociation timescale':
                      '1983A%26A...126..170C'})
        },
        'CO2': {
            'CE83': (5.0e5 * u.s,
                     {'CO2 photodissociation timescale':
                      '1983A%26A...126..170C'})
        },
        'CN': {
            'H92': ([3.15e5, 1.35e5] * u.s,
                    {'CN photodissociation timescale':
                     '1992Ap%26SS.195....1H'})
        },
    }

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

    if species is None:
        tab = Table(
            names=('Species', 'Source', 'Default', 'Lifetime_1 (s)',
                   'Lifetime_2 (s)', 'Bibcode'),
            dtype=('S6', 'S6', bool, float, float, 'S128'),
            masked=True)
        tab['Lifetime_2 (s)'].masked = True

        for species, sources in data.items():
            for source, (tau, bibcode) in sources.items():
                if np.size(tau) == 2:
                    tau1, tau2 = tau
                    mask = None
                else:
                    tau1 = tau
                    tau2 = 0
                    mask = [False, False, False, False, True, False]

                default = default_sources[species] == source
                tab.add_row((species, source, default, tau1, tau2, bibcode),
                            mask=mask)

        tab.pprint(max_lines=-1, max_width=-1)
        return

    if species.upper() not in data:
        raise ValueError(
            "No timescale available for {}.  Choose from: {}"
            .format(species, ', '.join(data.keys())))

    gas = data[species.upper()]
    source = default_sources[species.upper()] if source is None else source

    if source.upper() not in gas:
        raise ValueError(
            'Source key {} not available for {}.  Choose from: {}'
            .format(source, species, ', '.join(gas.keys())))

    tau, bibcode = gas[source.upper()]
    bib.register(photo_timescale, bibcode)

    return tau


def fluorescence_band_strength(species, eph=None, source=None):
    """Fluorescence band strength.


    Parameters
    ----------
    species : string
        The species to look up.

    rdot : `~astropy.units.Quantity`, optional
        Heliocentric radial speed, required for some species.

    eph : `~sbpy.data.Ephem`, optional
        The target ephemeris for species that require heliocentric
        radial velocity ('rdot').

    source : string, optional
        Retrieve values from this source (case insensitive).  See
        references for keys.


    Returns
    -------
    tau : `~astropy.units.Quantity`
        The timescale, scaled to `rh` or `eph['rh']`.


    Examples
    --------
    >>> from sbpy.activity import fluorescence_band_strength
    >>> LN = fluorescence_band_strength('OH')  # doctest: +SKIP

    References
    ----------
    [SA88] OH from Schleicher & A'Hearn 1988, ApJ 331, 1058-1077.
    Requires `rdot`.

    """

    raise NotImplemented

    # implement list treatment

    data = {   # (value, {key feature: bibcode})
        'OH 0-0': {
            'SA88': (func0_0,
                     {'OH 0-0 fluorescence band efficiency':
                      '1988ApJ...331.1058S'})
        },
        'OH 1-0': {
            'SA88': (func1_0,
                     {'OH 1-0 fluorescence band efficiency':
                      '1988ApJ...331.1058S'})
        },
        'OH 1-1': {
            'SA88': (func1_1,
                     {'OH 1-1 fluorescence band efficiency':
                      '1988ApJ...331.1058S'})
        },
        'OH 2-2': {
            'SA88': (func2_2,
                     {'OH 2-2 fluorescence band efficiency':
                      '1988ApJ...331.1058S'})
        },
    }

    default_sources = {
        'OH 0-0': 'SA88',
        'OH 1-0': 'SA88',
        'OH 1-1': 'SA88',
        'OH 2-2': 'SA88',
    }

    if species.upper() not in data:
        raise ValueError(
            'No data available for {}.  Choose one of: {}'
            .format(species, ', '.join(data.keys())))

    band = data[species.upper()]

    if source.upper() not in band:
        raise ValueError(
            'No source {} for {}.  Choose one of: {}'
            .format(source, species, ', '.join(band.keys())))

    LN, bibcode = band[source.upper()]
    bib.register(fluorescence_band_strength, bibcode)

    something_about_rdot_here

    return LN


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

        eph : dictionary-like, `~sbpy.data.Ephem`, `~astropy.units.Quantity`
            Target-observer distance, or ephemeris with `delta`.
            Required if the aperture has angular units.


        Returns
        -------
        sigma : `~astropy.units.Quantity`
            Coma column density along the line of sight at a distance
            rho.

        """

        rho_m = rho_as_length(rho, eph=eph).to('m').value
        return self._column_density(rho_m) * u.m**2

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

        return self._integrate_column_density(aper.as_length(eph))[0]

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
        pass

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
        pass

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

        """

        if not scipy:
            raise AstropyWarning(
                'scipy is required for integrating volume density.')

        if not rho.unit.is_equivalent(u.m):
            raise ValueError('rho must have units of length.')

        def f(s, rho2):
            r = np.sqrt(rho2 + s**2)
            return self._volume_density(r)

        # Without points, quad diverges.
        points = rho * np.logspace(-4, 4)
        sigma, err = quad(f, 0, np.inf, args=(rho**2,), points=points,
                          epsabs=epsabs)

        # spherical symmetry
        sigma *= 2

        return sigma

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
            raise AstropyWarning(
                'scipy is required for integrating column density')

        if not aper.dim.unit.is_equivalent(u.m):
            raise ValueError('aper must have units of length')

        if isinstance(aper, (CircularAperture, AnnularAperture)):
            if isinstance(aper, CircularAperture):
                limits = (0, aper.radius.to('m').value)
            else:
                limits = aper.shape.to('m').value

            # integrate in polar coordinates
            def f(rho):
                # rho in m, column_density in m**-2
                return rho * self._column_density(rho)

            N, err = quad(f, *limits, epsabs=epsabs)
            N *= 2 * np.pi
            err *= 2 * np.pi
        elif isinstance(aper, RectangularAperture):
            shape = aper.shape.to('m').value

            # integrate in polar coordinates
            def f(rho, th):
                return rho * self._column_density(rho)

            # first "octant"; g and h are the limits of the
            # integration of rho
            def g(th):
                return 0

            def h(th):
                return shape[0] / 2 / np.cos(th)

            th = np.arctan(shape[1] / shape[0])
            N1, err1 = dblquad(f, 0, th, g, h, epsabs=epsabs)

            # second "octant"
            def g(th):
                return 0

            def h(th):
                return shape[1] / 2 / np.cos(th)

            th = np.arctan(shape[0] / shape[1])
            N2, err2 = dblquad(f, 0, th, g, h, epsabs=epsabs)

            # N1 + N2 constitute 1/4th of the rectangle
            N = 4 * (N1 + N2)
            err = 4 * (err1 + err2)
        elif isinstance(aper, GaussianAperture):
            # integrate in polar coordinates
            def f(rho, sigma):
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
            raise AstropyWarning('scipy is not present, cannot continue.')
        return special.iti0k0(x.decompose().value)[1]

    def _K1(self, x):
        """Modified Bessel function of 2nd kind, 1st order."""
        if not scipy:
            raise AstropyWarning('scipy is not present, cannot continue.')
        return special.k1(x.decompose().value)

    @bib.cite({'model': '1978Icar...35..360N'})
    def _column_density(self, rho):
        sigma = (self.Q / self.v).to('1/m').value / r / 2 / np.pi
        parent = self.parent.to('m').value
        if self.daughter is None or self.daughter == 0:
            sigma *= np.pi / 2 - self._iK0(r / parent)
        else:
            daughter = self.daughter.to('m').value
            sigma *= (daughter / (parent - daughter)
                      * (self._iK0(r / daughter) - self._iK0(r / parent)))
        return sigma

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
    def total_number(self, aper, eph=None):
        if isinstance(aper, u.Quantity):
            aper = CircularAperture(aper)

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
        N = (self.Q * rho / self.v).to('dimensionless_unscaled').value
        if self.daughter is None or self.daughter == 0:
            x = rho / parent
            N *= 1 / x - self._K1(x) + np.pi / 2 - self._iK0(x)
        else:
            daughter = self.daughter.to(rho.unit)
            y = rho / daughter
            N *= (daughter / (parent - daughter)
                  * (self._iK0(y) - self._iK0(x) + x**-1 - y**-1
                     + self._K1(y) - self._K1(x)))

        return N
    total_number.__doc__ = GasComa.total_number.__doc__
