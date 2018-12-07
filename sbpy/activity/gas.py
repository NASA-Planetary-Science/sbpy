# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
========================
SBPy Activity Gas Module
========================

Functions
---------
photo_lengthscale          - Photodissociation lengthscale.
photo_timescale            - Photodissociation timescale.
fluorescence_band_strength - Fluorescence band efficiency of a specific
                             species and transition.

Classes
-------
Activity            - Abstract base class for gas coma models.
Haser               - Haser coma model for gas (Haser 1957).
Vectorial           - Vectorial coma model for gas (Festou 1981).


"""

from abc import ABC, abstractmethod
import numpy as np
import astropy.units as u
from .. import bib


__all__ = [
    'photo_lengthscale',
    'photo_timescale',
    'fluorescence_band_strength',

    'Haser',
    'Vectorial',
]


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
    gamma : astropy Quantity
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
        'H2O': {'CS93': (2.4e4 * u.km, {'H2O photodissociation lengthscale': '1993Icar..105..235C'})},
        'OH': {'CS93': (1.6e5 * u.km, {'OH photodissociation lengthscale': '1993Icar..105..235C'})},

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
    bib.register('activity.gas.photo_lengthscale', bibcode)

    return gamma


def photo_timescale(species, source=None):
    """Photodissociation timescale for a gas species.

    Parameters
    ----------
    species : string or None
      The species to look up, or `None` to summarize available
      species.

    source : string, optional
      Retrieve values from this source (case insensitive).  See
      references for keys.


    Returns
    -------
    tau : astropy Quantity
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
        'H2O': {'CS93': (5.2e4 * u.s, {'H2O photodissociation timescale': '1993Icar..105..235C'})},
        'OH': {'CS93': (1.6e5 * u.s, {'OH photodissociation timescale': '1993Icar..105..235C'})},
        'HCN': {'C94': (6.7e4 * u.s, {'HCN photodissociation timescale': '1994JGR....99.3777C'})},
        'CH3OH': {'C94': (7.7e4 * u.s, {'CH3OH photodissociation timescale': '1994JGR....99.3777C'})},
        'H2CO': {'C94': (5.0e3 * u.s, {'H2CO photodissociation timescale': '1994JGR....99.3777C'})},
        'CO': {'CE83': (1.5e6 * u.s, {'CO photodissociation timescale': '1983A%26A...126..170C'})},
        'CO2': {'CE83': (5.0e5 * u.s, {'CO2 photodissociation timescale': '1983A%26A...126..170C'})},
        'CN': {'H92': ([3.15e5, 1.35e5] * u.s, {'CN photodissociation timescale': '1992Ap%26SS.195....1H'})},
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
        from astropy.table import Table

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

                default = True if default_sources[species] == source else False
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
    bib.register('activity.gas.photo_timescale', bibcode)

    return tau


def fluorescence_band_strength(species, rdot=0 * u.km / u.s,
                               eph=None, source=None):
    """Fluorescence band efficiency of a specific species and transition.

    Parameters
    ----------
    species : string
      The species to look up.
    rdot : astropy Quantity, optional
      Heliocentric radial speed, required for some species.
    eph : sbpy Ephem, optional
      The target ephemeris.  Must include heliocentric radial
      velocity.
    source : string, optional
      Retrieve values from this source (case insensitive).  See
      references for keys.

    Returns
    -------
    tau : astropy Quantity
      The timescale, scaled to `rh` or `eph['rh']`.

    Notes
    -----
    One of `rdot` or `eph` is required for some species.

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
        'OH 0-0': {'SA88': ('XXX', {'OH 0-0 fluorescence band efficiency': '1988ApJ...331.1058S'})},
        'OH 1-0': {'SA88': ('XXX', {'OH 1-0 fluorescence band efficiency': '1988ApJ...331.1058S'})},
        'OH 1-1': {'SA88': ('XXX', {'OH 1-1 fluorescence band efficiency': '1988ApJ...331.1058S'})},
        'OH 2-2': {'SA88': ('XXX', {'OH 2-2 fluorescence band efficiency': '1988ApJ...331.1058S'})},
    }

    default_sources = {
        'OH 0-0': ('model', 'SA88'),
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
    bib.register('activity.gas.fluorescence_band_strength', bibcode)

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

    def __init__(self, Q, v):
        if not Q.unit.is_equivalent((u.s**-1, u.mol / u.s)):
            raise ValueError('Q must have units equivalent to 1/s or mol/s')
        self.Q = Q

        if not v.unit.is_equivalent(u.m / u.s):
            raise ValueError('v must have units of length/time')
        self.v = v

    @abstractmethod
    def volume_density(self, r):
        """Coma volume density.

        Parameters
        ----------
        r : `~astropy.units.Quantity`
          Linear distance to the nucleus.


        Returns
        -------
        n : float

        """
        pass

    @abstractmethod
    def column_density(self, rho, eph=None):
        """Coma column density at a projected distance from nucleus.

        Parameters
        ----------
        rho : `~astropy.units.Quantity`
          Projected distance of the region of interest on the plane of
          the sky in units of length or angle.

        eph : dictionary-like or `~sbpy.data.Ephem`
          Ephemerides at epoch; requires geocentric distance as
          `delta` keyword if aperture has angular units.


        Returns
        -------
        sigma : float
          Coma column density along the line of sight at a distance rho.

        """
        pass

    def _integrate_volume_density(self, rho, epsabs=1.49e-8):
        """
        Integrate volume density along the line of sight.

        Parameters
        ----------
        rho : `~astropy.units.Quantity`
            Projected distance of the region of interest on the plane of
            the sky in units of length.
        epsabs : float or int, optional
            Absolute and relative error tolerance for integrals.  See
            `scipy.integrate.quad`.

        Returns
        -------
        sigma : float
          Coma column density along the line of sight at a distance rho.

        """

        try:
            from scipy.integrate import quad
        except ImportError:
            from astropy.utils.exceptions import AstropyWarning
            from warnings import warn
            warn(AstropyWarning(
                'scipy is not present, cannot integrate volume density.'))
            return None

        if not rho.unit.is_equivalent(u.m):
            raise ValueError('rho must have units of length.')

        def f(s):
            r = np.sqrt(rho.to(u.km).value**2 + s**2)
            n = self.volume_density(r*u.km) * u.km
            return n.decompose().value

        # Using an upper limit of integration than 1e9 m makes the
        # integral divergent
        # sigma, err = quad(f, 0, np.inf, epsabs=epsabs)
        sigma, err = quad(f, 0, np.max(
            (1.e6, 10*rho.to(u.km).value)), epsabs=epsabs)

        # spherically symmetric coma
        sigma *= 2

        return sigma

    @abstractmethod
    def total_number(self, aper, eph=None):
        """Total number of molecules in aperture.

        Parameters
        ----------
        aper : `~astropy.units.Quantity` or `~sbpy.activity.Aperture`
          Observation aperture as a radius for a circular aperture
          (projected length, or angle) or an `Aperture` instance.

        eph : dictionary-like or `~sbpy.data.Ephem`, optional
          Ephemerides at epoch; requires geocentric distance as
          `delta` keyword if aperture has angular units.


        Returns
        -------
        N : float
          Total number of molecules within the aperture.

        """
        pass

    def _integrate_column_density(self, aper, epsabs=1.49e-8):
        """Integrate column density over an aperture.

        Parameters
        ----------
        aper : `~sbpy.activity.Aperture`
          Aperture, in units of length.

        epsabs : float or int, optional
          Absolute and relative error tolerance for integrals.  See
          `scipy.integrate.quad` (circular, annular, Gaussian) and
          `scipy.integrate.dblquad` (rectangular) for details.

        """

        from .core import RectangularAperture, GaussianAperture, AnnularAperture, CircularAperture

        try:
            from scipy.integrate import quad, dblquad
        except ImportError as e:
            from astropy.utils.exceptions import AstropyWarning
            from warnings import warn
            warn(AstropyWarning(
                'scipy is not present, cannot integrate column density.'))
            return None

        if not aper.dim.unit.is_equivalent(u.m):
            raise ValueError('aper must have units of length')

        if isinstance(aper, CircularAperture):
            # integrate in polar coordinates
            def f(rho):
                x = rho * self.column_density(rho * u.km) * u.km**2
                return x.decompose().value

            N, err = quad(f, 0, aper.radius.to(u.km).value, epsabs=epsabs)
            N *= 2 * np.pi
        elif isinstance(aper, AnnularAperture):
            # integrate in polar coordinates
            def f(rho):
                x = rho * self.column_density(rho * u.km) * u.km**2
                return x.decompose().value

            N, err = quad(f, aper.shape[0].to(u.km).value,
                          aper.shape[1].to(u.km).value, epsabs=epsabs)
            N *= 2 * np.pi
        elif isinstance(aper, RectangularAperture):
            # integrate in polar coordinates
            def f(rho, th):
                x = rho * self.column_density(rho * u.km) * u.km**2
                return x.decompose().value

            shape = aper.shape.to(u.km).value

            # first "octant"; g and h are the limits of the
            # integration of rho
            def g(th): return 0

            def h(th): return shape[0] / 2 / np.cos(th)
            th = np.arctan(shape[1] / shape[0])
            N1, err1 = dblquad(f, 0, th, g, h, epsabs=epsabs)

            # second "octant"
            def g(th): return 0

            def h(th): return shape[1] / 2 / np.cos(th)
            th = np.arctan(shape[0] / shape[1])
            N2, err2 = dblquad(f, 0, th, g, h, epsabs=epsabs)

            # N1 + N2 constitute 1/4th of the rectangle
            N = 4 * (N1 + N2)
        elif isinstance(aper, GaussianAperture):
            # integrate in polar coordinates
            def f(rho): return (rho * aper(rho * u.km).value
                                * self.column_density(rho * u.km).to(u.km**-2).value)
            N, err = quad(f, 0, np.inf, epsabs=epsabs)
            N *= 2 * np.pi

        return N


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

    def __init__(self, Q, v, parent, daughter=None):
        super(Haser, self).__init__(Q, v)

        bib.register('activity.gas.Haser', {'model': '1957BSRSL..43..740H'})

        if not parent.unit.is_equivalent(u.m):
            raise ValueError('parent must have units of length')
        self.parent = parent

        if daughter is None:
            self.daughter = None
        else:
            if not daughter.unit.is_equivalent(u.m):
                raise ValueError('daugher must have units of length')
            self.daughter = daughter

    def volume_density(self, r):
        if not r.unit.is_equivalent(u.m):
            raise ValueError('r must have units of length')

        n = self.Q / 4 / np.pi / r**2 / self.v
        if self.daughter is None or self.daughter == 0:
            # parent only
            n *= np.exp(-r / self.parent)
        else:
            n *= (self.daughter / (self.parent - self.daughter)
                  * (np.exp(-r / self.parent) - np.exp(-r / self.daughter)))

        return n.decompose()

    def _iK0(self, x):
        """Integral of the modified Bessel function of 2nd kind, 0th order."""
        try:
            from scipy.special import iti0k0
        except ImportError as e:
            from astropy.utils.exceptions import AstropyWarning
            from warnings import warn
            warn(AstropyWarning('scipy is not present, cannot continue.'))
            return None

        return iti0k0(x.decompose().value)[1]

    def _K1(self, x):
        """Modified Bessel function of 2nd kind, 1st order."""
        try:
            from scipy.special import k1
        except ImportError as e:
            from astropy.utils.exceptions import AstropyWarning
            from warnings import warn
            warn(AstropyWarning('scipy is not present, cannot continue.'))
            return None

        return k1(x.decompose().value)

    def column_density(self, rho, eph=None):
        from .core import rho_as_length

        bib.register('activity.gas.Haser.column_density',
                     {'model': '1978Icar...35..360N'})

        r = rho_as_length(rho, eph=eph)
        x = 0 if self.parent is None else (r / self.parent).decompose()
        y = 0 if self.daughter is None else (r / self.daughter).decompose()
        sigma = self.Q / 2 / np.pi / r / self.v
        if self.daughter is None or self.daughter == 0:
            sigma *= np.pi / 2 - self._iK0(x)
        elif self.parent is None or self.parent == 0:
            sigma *= np.pi / 2 - self._iK0(y)
        else:
            sigma *= (self.daughter / (self.parent - self.daughter)
                      * (self._iK0(y) - self._iK0(x)))

        return sigma.decompose()

    def total_number(self, aper, eph=None):
        from .core import rho_as_length, Aperture
        from .core import RectangularAperture, GaussianAperture, AnnularAperture, CircularAperture

        bib.register('activity.gas.Haser.total_number',
                     {'model': '1978Icar...35..360N'})

        # Inspect aper and handle as appropriate
        if isinstance(aper, Aperture):
            aper = aper.as_length(eph)
            if isinstance(aper, (RectangularAperture, GaussianAperture)):
                return self._integrate_column_density(aper)
            elif isinstance(aper, AnnularAperture):
                return self.total_number(aper.shape[1]) - self.total_number(aper.shape[0])
            elif isinstance(aper, CircularAperture):
                rho = aper.radius
            else:
                raise NotImplemented(
                    "Integration of {} apertures is not implemented.".format(type(aper)))
        else:
            rho = rho_as_length(aper, eph)

        # Solution for the circular aperture of radius rho:
        x = 0 if self.parent is None else (rho / self.parent).decompose()
        y = 0 if self.daughter is None else (rho / self.daughter).decompose()

        N = self.Q * rho / self.v
        if self.daughter is None or self.daughter == 0:
            N *= 1 / x - self._K1(x) + np.pi / 2 - self._iK0(x)
        elif self.parent is None or self.parent == 0:
            N *= 1 / y - self._K1(y) + np.pi / 2 - self._iK0(y)
        else:
            N *= (self.daughter / (self.parent - self.daughter)
                  * (self._iK0(y) - self._iK0(x) + x**-1 - y**-1
                     + self._K1(y) - self._K1(x)))

        return N.decompose().value


class Vectorial(GasComa):
    """Vectorial model implementation"""

    def __init__(self, Q, species):
        """Parameters
        ----------
        Q : `Astropy.units` quantity or iterable, mandatory
            production rate usually in units of `u.molecule / u.s`
        species : dictionary or list of dictionaries, mandatory
            defines gas velocity, lifetimes, disassociative lifetimes

        Returns
        -------
        Vectorial instance

        Examples
        --------
        TBD

        not yet implemented

        """

        pass
