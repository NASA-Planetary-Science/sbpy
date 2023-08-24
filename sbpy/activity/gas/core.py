# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""activity.gas core

"""

__all__ = [
    "photo_lengthscale",
    "photo_timescale",
    "fluorescence_band_strength",
    "Haser",
    "VectorialModel",
]

from abc import ABC, abstractmethod

# distutils is deprecated in python 3.10 and will be removed in 3.12 (PEP 632).
# Migration from distutils.log -> logging
from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np
import astropy.units as u

try:
    import scipy
    from scipy import special
    from scipy.integrate import quad, dblquad, romberg
    from scipy.interpolate import CubicSpline, PPoly
except ImportError:
    scipy = None
    PPoly = None

from ... import bib
from ... import data as sbd
from ... import units as sbu
from ...utils.decorators import requires
from ..core import (
    Aperture,
    RectangularAperture,
    GaussianAperture,
    AnnularAperture,
    CircularAperture,
)


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
        "H2O": "CS93",
        "OH": "CS93",
    }

    if species not in data:
        summary = ""
        for k, v in sorted(data.items()):
            summary += "\n{} [{}]".format(k, ", ".join(v.keys()))

        raise ValueError(
            "Invalid species {}.  Choose from:{}".format(species, summary))

    gas = data[species]
    source = default_sources[species] if source is None else source

    if source not in gas:
        raise ValueError(
            "Source key {} not available for {}.  Choose from: {}".format(
                source, species, ", ".join(gas.keys())
            )
        )

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
        "H2O": "CS93",
        "OH": "CS93",
        "HCN": "C94",
        "CH3OH": "C94",
        "H2CO": "C94",
        "CO2": "CE83",
        "CO": "CE83",
        "CN": "H92",
    }

    if species not in data:
        summary = ""
        for k, v in sorted(data.items()):
            summary += "\n{} [{}]".format(k, ", ".join(v.keys()))

        raise ValueError(
            "Invalid species {}.  Choose from:{}".format(species, summary))

    gas = data[species]
    source = default_sources[species] if source is None else source

    if source not in gas:
        raise ValueError(
            "Source key {} not available for {}.  Choose from: {}".format(
                source, species, ", ".join(gas.keys())
            )
        )

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
        "OH 0-0": "SA88",
        "OH 1-0": "SA88",
        "OH 1-1": "SA88",
        "OH 2-2": "SA88",
        "OH 0-1": "SA88",
        "OH 0-2": "SA88",
        "OH 2-0": "SA88",
        "OH 2-1": "SA88",
    }

    if species not in data:
        raise ValueError(
            "No data available for {}.  Choose one of: {}".format(
                species, ", ".join(data.keys())
            )
        )

    band = data[species]
    source = default_sources[species] if source is None else source

    if source not in band:
        raise ValueError(
            "No source {} for {}.  Choose one of: {}".format(
                source, species, ", ".join(band.keys())
            )
        )

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

        return self._volume_density(r.to_value("m")) / u.m**3

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, "delta"))
    def column_density(self, rho, eph=None):
        """Coma column density at a projected distance from nucleus.


        Parameters
        ----------
        rho : `~astropy.units.Quantity`
            Projected distance to the region of interest on the plane
            of the sky in units of length or angle.

        eph : dictionary-like, `~sbpy.data.Ephem`, `~astropy.units.Quantity`, optional
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

        rho = rho.to_value("m", equiv)
        return self._column_density(rho) / u.m**2

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, "delta"))
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

        if not isinstance(aper, Aperture):
            aper = CircularAperture(aper)

        if eph is not None:
            aper = aper.as_length(eph)

        return self._total_number(aper)

    def _total_number(self, aper):
        """Total number of molecules in aperture.

        Sub-classes of ``GasComa`` may override this method, instead of
        ``total_number``, which avoids reusing the boiler plate aperture
        conversions.


        Parameters
        ----------
        aper : `~sbpy.activity.Aperture`
            Observation aperture in units of length.

        """

        return self._integrate_column_density(aper)[0]

    @abstractmethod
    def _volume_density(self, r):
        """Unitless volume density function.


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

    @requires("scipy")
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

        def f(s, rho2):
            r = np.sqrt(rho2 + s**2)
            return self._volume_density(r)

        # quad diverges integrating to infinity, but 1e6 × rho is good
        # enough
        limit = 30
        points = rho * np.logspace(-4, 4, limit // 2)
        sigma, err = quad(
            f, 0, 1e6 * rho, args=(rho**2,), limit=limit, points=points, epsabs=epsabs
        )

        # spherical symmetry
        sigma *= 2
        err *= 2

        return sigma, err

    @requires("scipy")
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

        if isinstance(aper, (CircularAperture, AnnularAperture)):
            if isinstance(aper, CircularAperture):
                limits = (0, aper.radius.to_value("m"))
            else:
                limits = aper.shape.to_value("m")

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
            shape = aper.shape.to_value("m")

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
                return (
                    rho
                    * np.exp(-(rho**2) / sigma**2 / 2)
                    * self._column_density(rho)
                )

            sigma = aper.sigma.to_value("m")
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

    @bib.cite({"model": "1957BSRSL..43..740H"})
    @u.quantity_input(parent=u.m, daughter=u.m)
    def __init__(self, Q, v, parent, daughter=None):
        super().__init__(Q, v)
        self.parent = parent
        self.daughter = daughter

    def _volume_density(self, r):
        n = (self.Q / self.v).to_value("1/m") / r**2 / 4 / np.pi
        parent = self.parent.to_value("m")
        if self.daughter is None or self.daughter == 0:
            # parent only
            n *= np.exp(-r / parent)
        else:
            daughter = self.daughter.to_value("m")
            n *= (
                daughter
                / (parent - daughter)
                * (np.exp(-r / parent) - np.exp(-r / daughter))
            )

        return n

    def _iK0(self, x):
        """Integral of the modified Bessel function of 2nd kind, 0th order."""
        return special.iti0k0(x)[1]

    def _K1(self, x):
        """Modified Bessel function of 2nd kind, 1st order."""
        return special.k1(x)

    @requires("scipy")
    @bib.cite({"model": "1978Icar...35..360N"})
    def _column_density(self, rho):
        sigma = (self.Q / self.v).to_value("1/m") / rho / 2 / np.pi
        parent = self.parent.to_value("m")
        if self.daughter is None or self.daughter == 0:
            sigma *= np.pi / 2 - self._iK0(rho / parent)
        else:
            daughter = self.daughter.to_value("m")
            sigma *= (
                daughter
                / (parent - daughter)
                * (self._iK0(rho / daughter) - self._iK0(rho / parent))
            )
        return sigma

    @requires("scipy")
    def _total_number(self, aper):
        # Inspect aper and handle as appropriate
        if isinstance(aper, (RectangularAperture, GaussianAperture)):
            return self._integrate_column_density(aper)[0]
        elif isinstance(aper, AnnularAperture):
            N0 = self._total_number(CircularAperture(aper.shape[0]))
            N1 = self._total_number(CircularAperture(aper.shape[1]))
            return N1 - N0

        # Solution for the circular aperture of radius rho:
        bib.register(self.total_number, {"model": "1978Icar...35..360N"})

        rho = aper.radius
        parent = self.parent.to(rho.unit)
        x = (rho / parent).to_value(u.dimensionless_unscaled)
        N = (self.Q * rho / self.v).to_value(u.dimensionless_unscaled)
        if self.daughter is None or self.daughter == 0:
            N *= 1 / x - self._K1(x) + np.pi / 2 - self._iK0(x)
        else:
            daughter = self.daughter.to(rho.unit)
            y = (rho / daughter).to_value("")
            N *= (daughter / (parent - daughter)).to_value("") * (
                self._iK0(y)
                - self._iK0(x)
                + x**-1
                - y**-1
                + self._K1(y)
                - self._K1(x)
            )

        return N


@dataclass
class VMParent:
    """
    Physical information about the parent necessary for the vectorial model

    Parameters
    ----------
    v_outflow : float
        See VectorialModel documentation.

    tau_d : float
        See VectorialModel documentation.

    tau_T : float
        See VectorialModel documentation.

    sigma : float
        See VectorialModel documentation.
    """
    v_outflow: float
    tau_d: float
    tau_T: float
    sigma: float


@dataclass
class VMFragment:
    """
    Physical information about the fragment necessary for the vectorial model

    Parameters
    ----------
    v_photo : float
        See VectorialModel documentation.

    tau_T : float
        See VectorialModel documentation.
    """
    v_photo: float
    tau_T: float


@dataclass
class VMGridParams:
    """
    Vectorial model gridding parameters to control how finely the space around
    the comet should be gridded, and how detailed the outflow sampling should
    be.

    Parameters
    ----------
    radial_points : int
        See VectorialModel documentation.

    angular_points : int
        See VectorialModel documentation.

    radial_substeps : int
        See VectorialModel documentation.
    """
    radial_points: int
    angular_points: int
    radial_substeps: int


@dataclass
class VMParams:
    """
    Vectorial model parameters unrelated to physical inputs, dealing with the
    model's assumptions and detalis of the calculations.

    Parameters
    ----------
    parent_destrution_level : float
        See VectorialModel documentation.

    fragment_destruction_level : float
        See VectorialModel documentation.

    max_fragment_lifetimes : float
        See VectorialModel documentation.
    """
    parent_destruction_level: float
    fragment_destruction_level: float
    max_fragment_lifetimes: float


@dataclass
class VMFragmentSputterPolar:
    """
    Describes the distribution of fragment volume density
    (``fragment_density[i][j]``) as function of (``rs[i]``, ``thetas[j]``) in a
    spherical coordinate system, given a column of parent molecules flowing
    outward along the azimuthal (z) axis.

    Parameters
    ----------
    rs : `np.ndarray`
        List of radii (`astropy.units.Quantity`).

    thetas: `np.ndarray`
        List of polar angles theta (radians) in a spherical coordinate system.

    fragment_density: `np.ndarray`
        List of fragment densities at the corresponding ``rs[i]`` and
        ``thetas[j]``.
    """
    rs: np.ndarray
    thetas: np.ndarray
    fragment_density: np.ndarray


@dataclass
class VMResult:
    """
    Dataclass to hold a collection of vectorial model details and results in a
    language- and model-independent way.

    Parameters
    ----------
    volume_density_grid : `np.ndarray`
        List of radii (`~astropy.units.Quantity`) that the model used for
        volume density calculations.

    volume_density : `np.ndarray`
        List of volume densities (`~astropy.units.Quantity`) computed at the
        corresponding radius: at radius = ``volume_density_grid[i]``, the volume
        density is ``volume_density[i]``.

    column_density_grid : `np.ndarray`
        List of radii (`~astropy.units.Quantity`) that the model used for
        column density calculations.

    column_density : `np.ndarray`
        List of column densities (`~astropy.units.Quantity`), computed at the
        corresponding radius: at radius = ``column_density_grid[i]``, the column
        density is ``column_density[i]``.

    fragment_sputter : `VMFragmentSputterPolar`
        Describes the distribution of fragment volume density as function of
        (r, theta) in a spherical coordinate system, given a column of parent
        molecules flowing outward along the azimuthal (z) axis.

    solid_angle_sputter : `VMFragmentSputterPolar`
        Similar to ``fragment_sputter``, but adjusted by a factor of
        ``sin(theta)``.

    volume_density_interpolation : `scipy.interpolate.PPoly`
        Function that takes radius as a float in meters (no astropy units) and
        returns the volume density in 1/m^3.  Not reliable beyond
        ``max_grid_radius``, and questionable for radii less than
        ``collision_sphere_radius``.

    column_density_interpolation : `scipi.interpolate.PPoly`
        Function that takes radius as a float in meters (no astropy units) and
        returns the column density in 1/m^2.  Not relaible beyond
        ``max_grid_radius``, and questionable for radii less than
        ``collision_sphere_radius``.

    collision_sphere_radius : `astropy.units.Quantity`
        The radius of the collision sphere as described and calculated in
        Festou (1981): the radius beyond which a molecule can expect to see one
        collision over its lifetime.

    max_grid_radius : `astropy.units.Quantity`
        Cutoff radius for the model calculations.  Attempting to take results
        from the model beyond this radius will be wrong or unreliable at best.

    coma_radius : `astropy.units.Quantity`
        Cutoff radius beyond which we take the parents to be entirely
        dissociated.

    num_fragments_theory : `float`
        Total number of fragments that we expect based on a steady-state
        calculation.  Unreliable when time dependence is strong.

    num_fragments_grid : `float`
        Total number of fragments based on the actual results of the model.
        Agreement with ``num_fragments_theory`` is generally very good in the
        steady-state case.

    t_perm_flow : `float`
        Time for the comet to reach a steady state/permanent flow regime, where
        the number of fragments produced is equal to the number of fragments
        lost to dissociation.  In a time-dependent production context, this
        will also measure how long the effects of a single outburst can affect
        the density.
    """
    volume_density_grid: np.ndarray = None
    volume_density: np.ndarray = None
    column_density_grid: np.ndarray = None
    column_density: np.ndarray = None

    fragment_sputter: VMFragmentSputterPolar = None
    solid_angle_sputter: VMFragmentSputterPolar = None

    volume_density_interpolation: PPoly = None
    column_density_interpolation: PPoly = None

    collision_sphere_radius: u.Quantity = None
    max_grid_radius: u.Quantity = None
    coma_radius: u.Quantity = None

    num_fragments_theory: float = None
    num_fragments_grid: float = None
    t_perm_flow: u.Quantity = None


class VectorialModel(GasComa):
    """
    Vectorial model for fragments in a coma produced with a dissociative energy
    kick.


    Parameters
    ----------
    base_q : `~astropy.units.Quantity`
        Base production rate, per time

    parent: `~sbpy.data.Phys`
        Object with the following physical property fields:
            * ``tau_T``: total lifetime of the parent molecule
            * ``tau_d``: photodissociative lifetime of the parent molecule
            * ``v_outflow``: outflow velocity of the parent molecule
            * ``sigma``: cross-sectional area of the parent molecule

    fragment: `~sbpy.data.Phys`
        Object with the following physical property fields:
            * ``tau_T``: total lifetime of the fragment molecule
            * ``v_photo``: velocity of fragment resulting from
                photodissociation of the parent

    q_t: callable, optional
        Calculates the parent production rate as a function of time: ``q_t(t)``.
        The argument ``t`` is the look-back time as a float in units of seconds.
        The return value is the production rate in units of inverse seconds.  If
        provided, this value is added to ``base_q``.

        If no time-dependence function is given, the model will run with steady
        production at ``base_q`` stretching infinitely far into the past.

    radial_points: int, optional
        Number of radial grid points the model will use

    radial_substeps: int, optional
        Number of points along the outflow axis to integrate over

    angular_points: int, optional
        Number of angular grid points the model will use

    parent_destruction_level: float, optional
        Model will attempt to track parents until this percentage has
        dissociated

    fragment_destruction_level: float, optional
        Model will attempt to track fragments until this percentage has
        dissociated

    max_fragment_lifetimes: float, optional
        Fragments traveling through the coma will be ignored if they take longer
        than this to arrive and contribute to the density at any considered
        point.

    print_progress: bool, optional
        Print progress while calculating.


    References
    ----------
    The density distribution of neutral compounds in cometary atmospheres. I -
    Models and equations, Festou, M. C. 1981, Astronomy and Astrophysics, vol.
    95, no. 1, Feb. 1981, p. 69-79.

    """

    @bib.cite({"model": "1981A&A....95...69F"})
    @u.quantity_input(base_q=(u.s**-1, u.mol / u.s))
    def __init__(
        self,
        base_q,
        parent,
        fragment,
        q_t=None,
        radial_points=50,
        radial_substeps=80,
        angular_points=30,
        parent_destruction_level=0.99,
        fragment_destruction_level=0.95,
        max_fragment_lifetimes=8.0,
        print_progress=False,
    ):
        super().__init__(base_q, parent["v_outflow"][0])

        # Calculations are done internally in meters and seconds to match the
        # base GasComa class

        # Convert to unitless value of production per second
        self.base_q = base_q.to(1 / u.s).value

        # Copy time dependence or create a steady production function
        if q_t is None:
            self.q_t = self._make_steady_production()
        else:
            self.q_t = q_t

        # Copy parent info, stripping astropy units and converting to meters
        # and seconds
        self.parent = VMParent(
            tau_T=parent["tau_T"][0].to(u.s).value,
            tau_d=parent["tau_d"][0].to(u.s).value,
            v_outflow=parent["v_outflow"][0].to(u.m / u.s).value,
            sigma=parent["sigma"][0].to(u.m**2).value,
        )

        # Same for the fragment info
        self.fragment = VMFragment(
            tau_T=fragment["tau_T"][0].to(u.s).value,
            v_photo=fragment["v_photo"][0].to(u.m / u.s).value,
        )

        # Grid settings
        self.grid = VMGridParams(
            radial_points=radial_points,
            radial_substeps=radial_substeps,
            angular_points=angular_points,
        )

        self.model_params = VMParams(
            # Helps define cutoff for radial grid at this percentage of parents
            # lost to decay
            parent_destruction_level=parent_destruction_level,
            # Default here is lower than parents because they are born farther
            # from nucleus, tracking them too long will stretch the radial grid
            # a bit too much
            fragment_destruction_level=fragment_destruction_level,
            # If a fragment has to travel longer than this many lifetimes to
            # contribute to the density at a point, ignore it
            max_fragment_lifetimes=max_fragment_lifetimes,
        )

        self.print_progress = print_progress

        self.vmr = VMResult()

        # Calculate up a few things
        self._setup_calculations()

        # Build the radial grid
        self.fast_voldens_grid = self._make_radial_logspace_grid()
        self.vmr.volume_density_grid = self.fast_voldens_grid * (u.m)

        # Angular grid
        self.d_alpha = self.epsilon_max / self.grid.angular_points
        # Make array of angles adjusted up away from zero, to keep from
        # encroaching on the outflow axis
        self.angular_grid = np.linspace(
            0, self.epsilon_max, num=self.grid.angular_points, endpoint=False
        )
        self.angular_grid += self.d_alpha / 2

        # Makes a 2d array full of zero values
        self.fragment_sputter = np.zeros(
            (self.grid.radial_points, self.grid.angular_points)
        )

        # Do the main computation
        self._compute_fragment_density()

        # Turn our volume density into a function of r
        self._interpolate_volume_density()

        # Get the column density at our grid points and interpolate
        self._compute_column_density()

        # Count up the number of fragments in the grid versus theoretical value
        self.vmr.num_fragments_theory = self._calc_num_fragments_theory()
        self.vmr.num_fragments_grid = self._calc_num_fragments_grid()

        # Convert fragment sputters to group of one-dimensional arrays
        # of the form r_i, theta_i, fragment_density_i
        sputterlist = []
        for (i, j), frag_dens in np.ndenumerate(self.fragment_sputter):
            sputterlist.append(
                [
                    self.vmr.volume_density_grid[i].to(u.m).value,
                    self.angular_grid[j],
                    frag_dens,
                ]
            )
        sputter = np.array(sputterlist)
        # fill in the fragment sputter results
        self.vmr.fragment_sputter = VMFragmentSputterPolar(
            rs=sputter[:, 0] * u.m,
            thetas=sputter[:, 1],
            fragment_density=sputter[:, 2] / u.m**3,
        )
        normsputterlist = []
        for (i, j), norm_dens in np.ndenumerate(self.solid_angle_sputter):
            normsputterlist.append(
                [
                    self.vmr.volume_density_grid[i].to(u.m).value,
                    self.angular_grid[j],
                    norm_dens,
                ]
            )
        normsputter = np.array(normsputterlist)
        self.vmr.solid_angle_sputter = VMFragmentSputterPolar(
            rs=normsputter[:, 0] * u.m,
            thetas=normsputter[:, 1],
            fragment_density=normsputter[:, 2] / u.m**3,
        )

        if print_progress:
            print("Vectorial model calculations complete!")

    @classmethod
    def binned_production(cls, qs, fragment, parent, ts, **kwargs):
        """Alternate constructor for vectorial model


        Parameters
        ----------
        qs : `~astropy.units.Quantity`
            List of steady production rates, per time, with length equal to that
            of ``ts``.

        parent: `~sbpy.data.Phys`
            Same as __init__

        fragment: `~sbpy.data.Phys`
            Same as __init__

        ts : `~astropy.units.Quantity`
            List of times corresponding to when the production qs begin, with
            positive times indicating the past.

        kwargs: variable, optional
            Any additional parameters in kwargs are passed on to __init__, which
            are documented above and may be passed in here.


        Returns
        -------
        VectorialModel
            Instance of the VectorialModel class


        Examples
        --------
        This specifies that from 30 days ago to 7 days ago, the production was
        1.e27, changes to 3.e27 between 7 and 5 days ago, then falls to 2.e27
        from 5 days ago until now: >>> q_example = [1.e27, 3.e27, 1.e27] *
        (1/u.s) >>> t_example = [30, 7, 5] * u.day


        Notes
        -----
        Preserves Festou's original fortran method of describing time dependence
        in the model - time bins of steady production at specified intervals.

        The base production of the model is taken from the first element in the
        production array, which assumes the arrays are time-ordered from oldest
        to most recent.  The base production extends backward in time to
        infinity, so take care when using this method for time dependence if
        that is not what is intended.

        """

        return cls(
            base_q=qs[0],
            q_t=VectorialModel._make_binned_production(qs, ts),
            fragment=fragment,
            parent=parent,
            **kwargs
        )

    def _make_binned_production(qs, ts) -> Callable[[np.float64], np.float64]:
        """Produces a time dependence function out of lists given to
        binned_production constructor.


        Parameters
        ----------
        qs : `astropy.units.Quantity`
            See binned_production for description

        ts : `astropy.units.Quantity`
            See binned_production for description


        Returns
        -------
        q_t : function
            See __init__ for description


        Notes
        -----
        We create a model-compatible function for time dependence out of the
        information specified in the arrays qs and ts The resulting function
        gives a steady specified production within the given time windows.

        """

        base_q = qs[0]
        q_variations = qs - base_q

        q_invsecs = q_variations.to_value(1 / (u.s))
        t_at_p_secs = ts.to_value(u.s)

        # extend the arrays to simplify our comparisons in binned_q_t
        # last bin stops at observation time, t = 0
        t_at_p_secs = np.append(t_at_p_secs, [0])
        # junk value for production because this is in the future
        q_invsecs = np.append(q_invsecs, [0])

        # this function represents the deviation from base_q at any time t, so
        # we need to handle times outside of the range specified in 'ts'
        def q_t(t):
            # too long in the past?
            if t > t_at_p_secs[0]:
                return 0
            # in the future?
            if t < 0:
                return 0

            # find which bin the given time falls in, and return the
            # corresponding production for that interval
            for i in range(len(q_invsecs) - 1):
                if t < t_at_p_secs[i] and t > t_at_p_secs[i + 1]:
                    return q_invsecs[i]

        return q_t

    def _make_steady_production(self) -> Callable[[np.float64], np.float64]:
        """Produces a time dependence function that contributes no extra
        parents at any time.


        Returns
        -------
        q_t : function
            See __init__ for description


        Notes
        -----
        If no q_t is given, we use this as our time dependence as the
        model needs a q_t to run

        """

        # No additional production at any time
        def q_t(_):
            return 0.0

        return q_t

    def _setup_calculations(self) -> None:
        """Miscellaneous calculations to inform the model later.

        Notes
        -----
        Calculates the collision sphere radius, coma radius, time to
        permanent flow regime, the maximum radius our grid could possibly
        need to extend out to, and the maximum angle that a fragment's
        trajectory can deviate from its parent's trajectory (which is
        assumed to be radial).

        """

        """
            Calculate collision sphere radius based on the first production
            rate, Eq. (5) in Festou 1981

            Note that this is only calculated with the base production rate,
            because it is assumed that the base production rate has had roughly
            enough time to reach a steady state before letting production vary
            with time.
        """

        if self.print_progress:
            print("Performing setup calculations...")

        # This v_therm factor comes from molecular flux of ideal gas moving
        # through a surface, in our case the surface of the collision sphere
        # The gas is collisional inside this sphere and free beyond it, so we
        # can model this as effusion
        v_therm = self.parent.v_outflow * 0.25
        q = self.base_q
        vp = self.parent.v_outflow
        vf = self.fragment.v_photo

        # Eq. 5 of Festou 1981
        self.vmr.collision_sphere_radius = (
            (self.parent.sigma * q * v_therm) / (vp * vp)
        ) * u.m

        # Calculates the radius of the coma (defined as the space the parents
        # occupy)
        # NOTE: Equation (16) of Festou 1981 where alpha is the percent
        # destruction of molecules
        parent_beta_r = - \
            np.log(1.0 - self.model_params.parent_destruction_level)
        parent_r = parent_beta_r * vp * self.parent.tau_T
        self.vmr.coma_radius = parent_r * u.m

        fragment_beta_r = - \
            np.log(1.0 - self.model_params.fragment_destruction_level)
        # Calculate the time needed to hit a steady, permanent production of
        # fragments
        perm_flow_radius = self.vmr.coma_radius.value + (
            (vp + vf) * fragment_beta_r * self.fragment.tau_T
        )

        t_secs = self.vmr.coma_radius.value / vp + (
            perm_flow_radius - self.vmr.coma_radius.value
        ) / (vp + vf)
        self.vmr.t_perm_flow = (t_secs * u.s).to(u.day)

        # This is the total radial size that parents & fragments occupy, beyond
        # which we assume zero density
        self.vmr.max_grid_radius = perm_flow_radius * u.m

        # Two cases for angular range of ejection of fragment based on relative
        # velocities of parent and fragment species
        if vf < vp:
            self.epsilon_max = np.arcsin(vf / vp)
        else:
            self.epsilon_max = np.pi

    def production_at_time(self, t: np.float64) -> np.float64:
        """Get production rate at time t.


        Parameters
        ----------
        t : numpy.float64
            Time in seconds, with positive values representing the past


        Returns
        -------
        numpy.float64
            Production rate, unitless, at the specified time

        """

        return self.base_q + self.q_t(t)

    def _make_radial_logspace_grid(self) -> np.ndarray:
        """Create an appropriate radial grid based on the model parameters.


        Returns
        -------
        ndarray
            Logarithmically spaced samples of the radial space around the coma,
            out to a maximum distance.


        Notes
        -----
        Creates a grid (in meters) with numpy's logspace function that covers
        the expected radial size, stretching from 2 times the collision sphere
        radius (near the nucleus be dragons) out to the calculated max.  If we
        get too close to the nucleus things go very badly so don't do it, dear
        reader.

        """

        start_power = np.log10(self.vmr.collision_sphere_radius.value * 2)
        end_power = np.log10(self.vmr.max_grid_radius.value)
        return np.logspace(
            start_power, end_power, num=self.grid.radial_points, endpoint=True
        )

    def _outflow_axis_sampling(
        self, x: np.float64, y: np.float64, theta: np.float64
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Construct a list of points along the outflow axis, sampled to be
        more dense around the minimum distance to (x, y).


        Parameters
        ----------
        x : numpy.float64
            x-coordinate of the space where we are calculating the fragment
            density

        y : numpy.float64
            y-coordinate of the space where we are calculating the fragment
            density

        theta: numpy.float64
            polar spherical coordinate of (x, y) to avoid calculating it from
            (x, y)


        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray)
            Tuple with list of ejection sites in the first index and their
            extents in the second index


        Notes
        -----
        Returns list of radial points along outflow axis to sample for fragment
        ejection, with a list describing the size dr that each ejection point
        occupies along the outflow axis.

        """

        outflow_axis_edge = self.vmr.coma_radius.value

        # calc angle to edge of the grid, and use it if it's less than
        # epsilon max.  without this, when epsilon_max approaches pi we
        # get values of r outside the edge of the grid
        grid_edge_angle = np.arctan2(x, (y - outflow_axis_edge))
        if grid_edge_angle < self.epsilon_max:
            max_subangle = grid_edge_angle
        else:
            max_subangle = self.epsilon_max

        # Trace angles from our grid point at x, y and find where dissociations
        # would have to occur along the outflow axis. This samples shorter
        # trajectories more finely along the axis, which is necessary because
        # the exponential decay along their transit means they contribute the
        # most to the fragment density
        ces = np.vectorize(lambda epsilon: y - x / np.tan(epsilon))

        # array one element larger than radial_substeps because we will
        # use the midpoints defined by this array
        subangles = np.linspace(
            theta, max_subangle, num=self.grid.radial_substeps + 1, endpoint=True
        )

        # Find where this set of angles intersects the outflow axis
        ejection_grid = ces(subangles)
        # Find dr for radial integration
        drs = np.diff(ejection_grid)

        # This puts the sampling of the ejection sites at the midpoints
        # of our ejection grid
        ejection_sites = (ejection_grid[1:] + ejection_grid[:-1]) / 2

        return ejection_sites, drs

    def _fragment_sputter(self, r: np.float64, theta: np.float64) -> np.float64:
        """Compute the fragment density at (r, theta) in a spherical
        coordinate system where theta is the polar angle and the parents flow
        radially outward along the z axis.


        Parameters
        ----------
        r : numpy.float64
            radial coordinate in meters where we are calculating the fragment
            density.

        theta : numpy.float64
            polar spherical coordinate where we are calculating the fragment
            density.


        Returns
        -------
        np.float64
            Fragment density in 1/m**3, no astropy units attached.

        """

        sputter = 0.0
        vp = self.parent.v_outflow
        vf = self.fragment.v_photo
        p_tau_T = self.parent.tau_T
        f_tau_T = self.fragment.tau_T

        x = r * np.sin(theta)
        y = r * np.cos(theta)

        ejection_sites, drs = self._outflow_axis_sampling(x, y, theta)

        # Loop over these ejection sites that contribute to x, y
        for slice_r, dr in zip(ejection_sites, drs):
            # Distance from dissociation site to our grid point
            sep_dist = np.sqrt(x**2 + (slice_r - y) ** 2)

            cos_eject = (y - slice_r) / sep_dist
            sin_eject = x / sep_dist

            # Parent extinction when traveling along to the
            # dissociation site
            p_extinction = np.e ** (-slice_r / (p_tau_T * vp))

            # Calculate sqrt(vR^2 - u^2 sin^2 gamma)
            v_factor = np.sqrt(vf * vf - (vp * vp) * sin_eject**2)

            # The geometry of the problem can admit one or two solutions for
            # the velocity of the fragment
            v_one = vp * cos_eject + v_factor
            if vf > vp:
                velocities = [v_one]
            else:
                v_two = vp * cos_eject - v_factor
                velocities = [v_one, v_two]

            # TODO: this shouldn't be necessary if we only pick r, theta inside ejection
            if v_one == 0.0:
                continue

            for v in velocities:
                # Time taken to travel from the dissociation point at v, reject
                # if the time is too large and fragments have decayed beyond
                # self.fragment_destruction_level percent
                t_frag = sep_dist / v
                if t_frag > self.time_limit:
                    continue

                # total time between parent emission from nucleus and fragment
                # arriving at our point of interest, which we then use to look
                # up Q at that time in the past
                t_total = (slice_r / vp) + t_frag

                # Division by parent velocity makes this production per
                # unit distance for radial integration q(r, epsilon)
                # given by eq. 32, Festou 1981
                q = self.production_at_time(t_total) / vp
                q_r_eps = (v**2 * q) / (vf * np.abs(v - vp * cos_eject))

                # Fragment extinction when traveling at speed v from
                # dissociation site to x, y
                f_extinction = np.e ** (-t_frag / f_tau_T)

                # differential addition to the density integrating along dr,
                # similar to eq. (36) Festou 1981
                n_r = (p_extinction * f_extinction *
                       q_r_eps) / (sep_dist**2 * v)

                sputter += n_r * dr

        return sputter

    def _compute_fragment_density(self) -> None:
        """Computes the density of fragments as a function of radius.

        Notes
        -----
        Computes the density of fragments sputtered around the comet due to an
        outflow of parents along the z-axis, performing the integration in eq.
        (36), Festou 1981. The resulting radial fragment density will be in
        units of 1/(m^3) as we work in m, s, and m/s.

        We then interpolate the fragment density as a function of arbitrary
        radius.  We use our results to calculate the total number of fragments
        in the coma for comparison to the theoretical number we expect, to
        provide the user with a rough idea of how well the chosen radial and
        angular grid sizes have captured the appropriate amount of particles.
        Note that some level of disagreement is expected because the
        parent_destruction_level and fragment_destruction_level parameters cut
        the grid off before all parents can dissociate, and thus some escape the
        model and come up missing in the fragment count based on the grid.

        """

        if self.print_progress:
            print("Starting fragment density computations...")

        # Follow fragments until they have been totally destroyed
        self.time_limit = self.model_params.max_fragment_lifetimes * self.fragment.tau_T

        # More factors to fill out integral similar to eq. (36) Festou 1981
        integration_factor = (
            (1 / (4 * np.pi * self.parent.tau_d)) * self.d_alpha / (4.0 * np.pi)
        )

        # vectorize and apply to each combination of (r, theta)
        sputter_func = np.vectorize(self._fragment_sputter)
        rs, thetas = np.meshgrid(
            self.fast_voldens_grid, self.angular_grid, indexing="ij"
        )
        self.fragment_sputter = integration_factor * sputter_func(rs, thetas)

        # Make array to hold our data, no units
        self.fast_voldens = np.zeros(self.grid.radial_points)

        # Integration factors from angular part of integral, similar to
        # eq. (36) Festou 1981
        # integrate over theta to produce radial fragment volume density
        # axis = 1 sums over second column, theta
        # Equivalent to summing over j for sin(theta[j]) *
        # fragment_sputter[i][j] with numpy magic
        self.solid_angle_sputter = np.sin(thetas) * self.fragment_sputter
        self.fast_voldens = 2.0 * np.pi * \
            np.sum(self.solid_angle_sputter, axis=1)

        # Tag with proper units
        self.vmr.volume_density = self.fast_voldens / (u.m**3)

    @requires("scipy")
    def _interpolate_volume_density(self) -> None:
        """Interpolate the volume density as a function of radial distance
        from the nucleus.

        Takes our fragment density grid and constructs density as a function of
        arbitrary radius.

        """

        if self.print_progress:
            print("Interpolating radial fragment density...")

        # Interpolate this radial density grid with a cubic spline for lookup
        # at non-grid radii, input in m, output in 1/m^3
        self.vmr.volume_density_interpolation = CubicSpline(
            self.fast_voldens_grid, self.fast_voldens, bc_type="natural"
        )

    def _column_density_at_rho(self, rho: np.float64) -> np.float64:
        """
        Calculate the column density of fragments at a given impact parameter.


        Parameters
        ----------
        rho : np.float64
            Impact parameter of the column density integration, in meters


        Returns
        -------
        np.float64
            Column density at the given impact parameter in m^-2, no
            astropy units attached


        Notes
        -----
        We return zero column density beyond a certain distance past the
        grid edge, which can lead to strange graphing results and sharp
        cutoffs.

        """

        reach_factor = 2.0

        # _volume_density() approximates beyond the edge of the grid, so allow
        # an arbirtary limit beyond which we just return zero column density
        r_max = self.vmr.max_grid_radius.value * reach_factor

        if rho > r_max:
            return 0

        rhosq = rho**2
        z_max = np.sqrt(r_max**2 - rhosq)

        def column_density_integrand(z):
            return self._volume_density(np.sqrt(z**2 + rhosq))

        c_dens = 2 * romberg(column_density_integrand, 0,
                             z_max, rtol=0.0001, divmax=50)

        # result is in 1/m^2
        return c_dens

    @requires("scipy")
    def _compute_column_density(self) -> None:
        """Compute the column density on the grid and interpolate the results.

        Notes
        -----
        The interpolator returns column density in m^-2, no astropy units
        attached.

        """

        if self.print_progress:
            print("Computing column densities...")

        # make a grid for the column density values
        column_density_grid = self._make_radial_logspace_grid()
        # vectorize so we can hand the whole grid to our function
        cd_vectorized = np.vectorize(self._column_density_at_rho)
        # array now holds corresponding column density values
        column_densities = cd_vectorized(column_density_grid)

        self.fast_column_density_grid = column_density_grid
        self.vmr.column_density_grid = column_density_grid * u.m
        self.vmr.column_density = column_densities / (u.m**2)
        # Interpolation gives column density in m^-2
        self.vmr.column_density_interpolation = CubicSpline(
            column_density_grid, column_densities, bc_type="natural"
        )

    def _calc_num_fragments_theory(self) -> np.float64:
        """The total number of fragment species we expect in the coma.

        Returns
        -------
        np.float64
            Total number of fragment species we expect in the coma theoretically

        Notes
        -----
        Outbursts/time dependent production in general will make this result
        poor due to the grid being sized to capture a certain fraction of
        parents/fragments at the oldest (first) production rate.  The farther
        you get from this base production, the farther the model will deviate
        from capturing the requested percentage of particles.

        """

        vp = self.parent.v_outflow
        vf = self.fragment.v_photo
        p_tau_T = self.parent.tau_T
        f_tau_T = self.fragment.tau_T
        p_tau_d = self.parent.tau_d
        t_perm = self.vmr.t_perm_flow.to(u.s).value

        alpha = f_tau_T * p_tau_T / p_tau_d

        max_r = self.vmr.max_grid_radius.value
        edge_adjust = np.pi * max_r * max_r * (vf + vp) * self.fast_voldens[-1]

        num_time_slices = 1000

        # array starting at tperm (farthest back in time we expect something to
        # hang around), stepped down to zero (when the observation takes place)
        time_slices = np.linspace(t_perm, 0, num_time_slices, endpoint=True)

        # estimate based on discrete-time source and sink model
        theory_total = 0
        for i, t in enumerate(time_slices[:-1]):
            extinction_one = t / p_tau_T
            extinction_two = time_slices[i + 1] / p_tau_T
            mult_factor = -np.e ** (-extinction_one) + \
                np.e ** (-extinction_two)
            theory_total += self.production_at_time(t) * mult_factor

        return theory_total * alpha - edge_adjust

    def _calc_num_fragments_grid(self) -> np.float64:
        """Total number of fragments in the coma.

        Calculates the total number of fragments by integrating the density grid
        over its volume


        Returns
        -------
        np.float64
            Number of fragments in the coma based on our grid calculations


        Notes
        -----
        Outbursts/time dependent production in general will make this result
        poor due to the grid being sized to capture a certain fraction of
        parents/fragments at the oldest (first) production rate.  The farther
        you get from this base production, the farther the model will deviate
        from capturing the requested percentage of particles.

        """

        max_r = self.vmr.max_grid_radius.value

        def vol_integrand(r, r_func):
            return r_func(r) * r**2

        r_int = romberg(
            vol_integrand,
            0,
            max_r,
            args=(self.vmr.volume_density_interpolation,),
            rtol=0.0001,
            divmax=20,
        )
        return 4 * np.pi * r_int

    def _column_density(self, rho) -> np.float64:
        """Gives fragment column density at arbitrary impact parameter.

        Parameters
        ----------
        rho : np.float64
            Impact parameter, in meters, no astropy units attached.


        Returns
        -------
        np.float64
            Fragment column density at given impact parameter, in m^-2, no
            astropy units.

        """

        if rho < self.fast_column_density_grid[0]:
            return self.vmr.column_density[0].value
        if rho > self.vmr.max_grid_radius.value:
            return 0
        return self.vmr.column_density_interpolation(rho)

    def _volume_density(self, r) -> np.float64:
        """Gives fragment volume density at arbitrary radius.


        Parameters
        ----------
        r : np.float64
            Distance from nucleus, in meters, no astropy units attached.


        Returns
        -------
        np.float64
            Fragment volume density at specified radius, in m^-3, no astropy
            units.


        Notes
        -----
        When asked for a value at a radius smaller than the first grid point, we
        have two choices: return a value based on the interpolation of the
        volume density, or return the value of the closest grid point to the
        nucleus. This function returns the closest grid point to the nucleus,
        because the model does not say anything about what happens inside the
        collision sphere - so we approximate by clamping the value as constant.
        We can't trust the interpolation outside of the grid because the density
        values can get very large or become negative.

        Outside the radius of the grid, the volume density is approximated with
        exponential decay based on the total lifetime of the fragment.

        """

        if r < self.fast_voldens_grid[0]:
            return self.fast_voldens[0]
        if r > self.vmr.max_grid_radius.value:
            diff = r - self.vmr.max_grid_radius.value
            guess = self.fast_voldens[-1] * np.exp(
                -diff / (self.fragment.tau_T * self.fragment.v_photo)
            )
            return guess
        return self.vmr.volume_density_interpolation(r)
