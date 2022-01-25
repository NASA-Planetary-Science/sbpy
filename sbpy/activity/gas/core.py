# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""activity.gas core

"""

__all__ = [
    'photo_lengthscale',
    'photo_timescale',
    'fluorescence_band_strength',
    'Haser',
    'VectorialModel'
]

from abc import ABC, abstractmethod
from distutils.log import warn
import warnings

import numpy as np
import astropy.units as u

try:
    import scipy
    from scipy import special
    from scipy.integrate import quad, dblquad, romberg
    from scipy.interpolate import CubicSpline
except ImportError:
    scipy = None

from ... import bib
from ... import data as sbd
from ... import units as sbu
from ...exceptions import RequiredPackageUnavailable, TestingNeeded
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

    @u.quantity_input(Q=(u.s ** -1, u.mol / u.s), v=u.m / u.s)
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

        return self._volume_density(r.to_value('m')) / u.m ** 3

    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'delta'))
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

        rho = rho.to_value('m', equiv)
        return self._column_density(rho) / u.m ** 2

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
            r = np.sqrt(rho2 + s ** 2)
            return self._volume_density(r)

        # quad diverges integrating to infinity, but 1e6 × rho is good
        # enough
        limit = 30
        points = rho * np.logspace(-4, 4, limit // 2)
        sigma, err = quad(f, 0, 1e6 * rho, args=(rho ** 2,),
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
                limits = (0, aper.radius.to_value('m'))
            else:
                limits = aper.shape.to_value('m')

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
            shape = aper.shape.to_value('m')

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
                return (rho * np.exp(-rho ** 2 / sigma ** 2 / 2)
                        * self._column_density(rho))

            sigma = aper.sigma.to_value('m')
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
        n = (self.Q / self.v).to_value('1/m') / r ** 2 / 4 / np.pi
        parent = self.parent.to_value('m')
        if self.daughter is None or self.daughter == 0:
            # parent only
            n *= np.exp(-r / parent)
        else:
            daughter = self.daughter.to_value('m')
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
        sigma = (self.Q / self.v).to_value('1/m') / rho / 2 / np.pi
        parent = self.parent.to_value('m')
        if self.daughter is None or self.daughter == 0:
            sigma *= np.pi / 2 - self._iK0(rho / parent)
        else:
            daughter = self.daughter.to_value('m')
            sigma *= (daughter / (parent - daughter)
                      * (self._iK0(rho / daughter) - self._iK0(rho / parent)))
        return sigma

    def _total_number(self, aper):
        # Inspect aper and handle as appropriate
        if isinstance(aper, (RectangularAperture, GaussianAperture)):
            return self._integrate_column_density(aper)[0]
        elif isinstance(aper, AnnularAperture):
            N0 = self._total_number(CircularAperture(aper.shape[0]))
            N1 = self._total_number(CircularAperture(aper.shape[1]))
            return N1 - N0

        # Solution for the circular aperture of radius rho:
        bib.register(self.total_number, {'model': '1978Icar...35..360N'})

        rho = aper.radius
        parent = self.parent.to(rho.unit)
        x = (rho / parent).to_value(u.dimensionless_unscaled)
        N = (self.Q * rho / self.v).to_value(u.dimensionless_unscaled)
        if self.daughter is None or self.daughter == 0:
            N *= 1 / x - self._K1(x) + np.pi / 2 - self._iK0(x)
        else:
            daughter = self.daughter.to(rho.unit)
            y = (rho / daughter).to_value('')
            N *= ((daughter / (parent - daughter)).to_value('')
                  * (self._iK0(y) - self._iK0(x) + x ** -1 - y ** -1
                     + self._K1(y) - self._K1(x)))

        return N


class VectorialModel(GasComa):
    """ Vectorial model for fragments in a coma produced
         with a dissociative energy kick

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
            Calculates the parent production rate as a function of time:
            ``q_t(t)``. The argument ``t`` is the look-back time as a float
            in units of seconds.  The return value is the production rate in
            units of inverse seconds.  If provided, this value is added to
            ``base_q``.

            If no time-dependence function is given, the model will run with
            steady production at ``base_q`` stretching infinitely far into the
            past.

        radial_points: int, optional
            Number of radial grid points the model will use

        radial_substeps: int, optional
            Number of points along the contributing axis to integrate over

        angular_points: int, optional
            Number of angular grid points the model will use

        angular_substeps: int, optional
            Number of angular steps per radial substep to integrate over

        parent_destruction_level: float, optional
            Model will attempt to track parents until
            this percentage has dissociated

        fragment_destruction_level: float, optional
            Model will attempt to track fragments until
            this percentage has dissociated

        max_fragment_lifetimes: float, optional
            Fragments traveling through the coma will be ignored if they take
            longer than this to arrive and contribute to the density at any
            considered point

        print_progress: bool, optional
            Print progress percentage while calculating

        References:
            The density distribution of neutral compounds in cometary
            atmospheres. I - Models and equations,
            Festou, M. C. 1981, Astronomy and Astrophysics, vol. 95, no. 1,
            Feb. 1981, p. 69-79.
    """
    @bib.cite({'model': '1981A&A....95...69F'})
    @u.quantity_input(base_q=(u.s ** -1, u.mol / u.s))
    def __init__(self, base_q, parent, fragment, q_t=None, radial_points=50,
                 radial_substeps=12, angular_points=30, angular_substeps=7,
                 parent_destruction_level=0.99,
                 fragment_destruction_level=0.95,
                 max_fragment_lifetimes=8.0,
                 print_progress=False):

        warnings.warn("Literature tests with the Vectorial model are generally"
                      " in agreement at the 20% level or better.  The cause"
                      " for the differences with the Festou FORTRAN code are"
                      " not yet precisely known.  Help testing this feature is"
                      " appreciated.", TestingNeeded)

        super().__init__(base_q, parent['v_outflow'][0])

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
        self.parent = {
            'tau_T': parent['tau_T'][0].to(u.s).value,
            'tau_d': parent['tau_d'][0].to(u.s).value,
            'v_outflow': parent['v_outflow'][0].to(u.m / u.s).value,
            'sigma': parent['sigma'][0].to(u.m ** 2).value
        }

        # Same for the fragment info
        self.fragment = {
            'tau_T': fragment['tau_T'][0].to(u.s).value,
            'v_photo': fragment['v_photo'][0].to(u.m / u.s).value
        }

        # Grid settings
        self.radial_points = radial_points
        self.radial_substeps = radial_substeps
        self.angular_points = angular_points
        self.angular_substeps = angular_substeps

        # Helps define cutoff for radial grid at this percentage of parents
        # lost to decay
        self.parent_destruction_level = parent_destruction_level
        # Default here is lower than parents because they are born farther from
        # nucleus, tracking them too long will stretch the radial grid a bit
        # too much
        self.fragment_destruction_level = fragment_destruction_level

        # If a fragment has to travel longer than this many lifetimes to
        # contribute to the density at a point, ignore it
        self.max_fragment_lifetimes = max_fragment_lifetimes

        # Print progress during density calculations?
        self.print_progress = print_progress

        """Initialize data structures to hold our calculations"""
        self.vmodel = {}

        # Calculate up a few things
        self._setup_calculations()

        # Build the radial grid
        self.vmodel['fast_radial_grid'] = self._make_radial_logspace_grid()
        self.vmodel['radial_grid'] = self.vmodel['fast_radial_grid'] * (u.m)

        # Angular grid
        self.vmodel['d_alpha'] = self.vmodel[
            'epsilon_max'] / self.angular_points
        # Make array of angles adjusted up away from zero, to keep from
        # calculating a radial line's contribution to itself
        self.vmodel['angular_grid'] = np.linspace(
            0, self.vmodel['epsilon_max'], num=self.angular_points,
            endpoint=False
        )
        # This maps addition over the whole array automatically
        self.vmodel['angular_grid'] += self.vmodel['d_alpha'] / 2

        # Makes a 2d array full of zero values
        self.vmodel['density_grid'] = np.zeros((self.radial_points,
                                                self.angular_points))

        # Do the main computation
        self._compute_fragment_density()

        # Turn our grid into a function of r
        self._interpolate_radial_density()

        # Get the column density at our grid points and interpolate
        self._compute_column_density()

        # Count up the number of fragments in the grid versus theoretical value
        self.vmodel['num_fragments_theory'] = self.calc_num_fragments_theory()
        self.vmodel['num_fragments_grid'] = self.calc_num_fragments_grid()

    @classmethod
    def binned_production(cls, qs, fragment, parent, ts, **kwargs):
        """ Alternate constructor for vectorial model

            Parameters
            ----------
            qs : `~astropy.units.Quantity`
                List of steady production rates, per time, with length equal to
                that of ``ts``.

            parent: `~sbpy.data.Phys`
                Same as __init__

            fragment: `~sbpy.data.Phys`
                Same as __init__

            ts : `~astropy.units.Quantity`
                List of times corresponding to when the production qs begin,
                with positive times indicating the past.

            kwargs: variable, optional
                Any additional parameters in kwargs are passed on to __init__,
                which are documented above and may be passed in here.

            Returns
            -------
            VectorialModel
                Instance of the VectorialModel class

            Examples
            --------
            This specifies that from 30 days ago to 7 days ago, the production
            was 1.e27, changes to 3.e27 between 7 and 5 days ago, then falls to
            2.e27 from 5 days ago until now
            >>> q_example = [1.e27, 3.e27, 1.e27] * (1/u.s)
            >>> t_example = [30, 7, 5] * u.day

            Notes
            -----
            Preserves Festou's original fortran method of describing time
            dependence in the model - time bins of steady production at
            specified intervals

            The base production of the model is taken from the first element in
            the production array, which assumes the arrays are time-ordered
            from oldest to most recent.  The base production extends backward
            in time to infinity, so take care when using this method for time
            dependence if that is not what is intended.
        """
        return cls(base_q=qs[0],
                   q_t=VectorialModel._make_binned_production(qs, ts),
                   fragment=fragment, parent=parent, **kwargs)

    def _make_binned_production(qs, ts):
        """ Produces a time dependence function out of lists given to
            binned_production constructor

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
            We create a model-compatible function for time dependence out of
            the information specified in the arrays qs and ts
            The resulting function gives a steady specified production within
            the given time windows
        """

        q_invsecs = qs.to(1 / (u.s)).value
        t_at_p_secs = ts.to(u.s).value
        base_q = qs[0]

        # extend the arrays to simplify our comparisons in binned_q_t
        # last bin stops at observation time, t = 0
        t_at_p_secs = np.append(t_at_p_secs, [0])
        # junk value for production because this is in the future
        q_invsecs = np.append(q_invsecs, [0])

        def q_t(t):
            # too long in the past?
            if t > t_at_p_secs[0] or t < 0:
                return base_q
            # in the future?
            if t < 0:
                return 0

            # loop over all elements except the last so we can always look at
            # [index+1] for the comparison
            for i in range(len(q_invsecs) - 1):
                if t < t_at_p_secs[i] and t > t_at_p_secs[i + 1]:
                    return q_invsecs[i]

        return q_t

    def _make_steady_production(self):
        """ Produces a time dependence function that contributes no extra
            parents at any time

            Notes
            -----
            If no q_t is given, we use this as our time dependence as the
            model needs a q_t to run
        """
        # No additional production at any time
        def q_t(t):
            return 0

        return q_t

    def _setup_calculations(self):
        """ Miscellaneus calculations to inform the model later

            Notes
            -----
            Calculates the collision sphere radius, coma radius, time to
            permanent flow regime, the maximum radius our grid could possibly
            need to extend out to, and the maximum angle that a fragment's
            trajectory can deviate from its parent's trajectory (which is
            assumed to be radial)
        """

        """
            Calculate collision sphere radius based on the first production
            rate, Eq. (5) in Festou 1981

            Note that this is only calculated with the base production rate,
            because it is assumed that the base production rate has had roughly
            enough time to reach a steady state before letting production vary
            with time.
        """
        # This v_therm factor comes from molecular flux of ideal gas moving
        # through a surface, in our case the surface of the collision sphere
        # The gas is collisional inside this sphere and free beyond it, so we
        # can model this as effusion
        v_therm = self.parent['v_outflow'] * 0.25
        q = self.base_q
        vp = self.parent['v_outflow']
        vf = self.fragment['v_photo']

        # Eq. 5 of Festou 1981
        self.vmodel['collision_sphere_radius'] = (
            (self.parent['sigma'] * q * v_therm) / (vp * vp)
        ) * u.m

        # Calculates the radius of the coma (parents only), given our input
        # parameters
        # NOTE: Equation (16) of Festou 1981 where alpha is the percent
        # destruction of molecules
        parent_beta_r = -np.log(1.0 - self.parent_destruction_level)
        parent_r = parent_beta_r * vp * self.parent['tau_T']
        self.vmodel['coma_radius'] = parent_r * u.m

        # Calculates the time needed to hit a steady, permanent production of
        # fragments
        fragment_beta_r = -np.log(1.0 - self.fragment_destruction_level)

        perm_flow_radius = (
            self.vmodel['coma_radius'].value +
            ((vp + vf) * fragment_beta_r * self.fragment['tau_T'])
        )

        t_secs = (
            self.vmodel['coma_radius'].value / vp +
            (perm_flow_radius - self.vmodel['coma_radius'].value)
            / (vp + vf)
        )
        self.vmodel['t_perm_flow'] = (t_secs * u.s).to(u.day)

        # This is the total radial size that parents & fragments occupy, beyond
        # which we assume zero density
        self.vmodel['max_grid_radius'] = perm_flow_radius * u.m

        # Two cases for angular range of ejection of fragment based on relative
        # velocities of parent and fragment species
        if(vf < vp):
            self.vmodel['epsilon_max'] = np.arcsin(vf / vp)
        else:
            self.vmodel['epsilon_max'] = np.pi

    def production_at_time(self, t):
        """ Get production rate at time t

            Parameters
            ----------
            t : float
                Time in seconds, with positive values representing the past

            Returns
            -------
            numpy.float64
                Production rate, unitless, at the specified time

        """

        return self.base_q + self.q_t(t)

    def _make_radial_logspace_grid(self):
        """ Create an appropriate radial grid based on the model parameters

            Returns
            -------
            ndarray
                Logarithmically spaced samples of the radial space around the
                coma, out to a maximum distance

            Notes
            -----
            Creates a grid (in meters) with numpy's logspace function that
            covers the expected radial size, stretching from 2 times the
            collision sphere radius (near the nucleus be dragons) out to the
            calculated max.  If we get too close to the nucleus things go very
            badly so don't do it, dear reader
        """
        start_power = np.log10(
            self.vmodel['collision_sphere_radius'].value * 2
        )
        end_power = np.log10(self.vmodel['max_grid_radius'].value)
        return np.logspace(
            start_power, end_power,
            num=self.radial_points, endpoint=True
        )

    def _compute_fragment_density(self):
        """ Computes the density of fragments as a function of radius

            Notes
            -----
            Computes the density at different radii and due to each ejection
            angle, performing the radial integration of eq. (36), Festou 1981
            with only one fragment velocity.  The resulting units will be in
            1/(m^3) as we work in m, s, and m/s.

            The density is first found on a radial grid, then interpolated to
            find density as a function of arbitrary radius.  We use our results
            from the grid to calculate the total number of fragments in the
            coma for comparison to the theoretical number we expect, to provide
            the user with a rough idea of how well the chosen radial and
            angular grid sizes have captured the appropriate amount of
            particles.  Note that some level of disagreement is expected
            because the parent_destruction_level and fragment_destruction_level
            parameters cut the grid off before all particles can dissociate,
            and thus some escape the model and come up missing in the fragment
            count based on the grid.
        """
        vp = self.parent['v_outflow']
        vf = self.fragment['v_photo']

        # Follow fragments until they have been totally destroyed
        time_limit = self.max_fragment_lifetimes * self.fragment['tau_T']
        r_coma = self.vmodel['coma_radius'].value
        r_limit = r_coma

        # temporary radial array for when we loop through 0 to epsilonMax
        ejection_radii = np.zeros(self.radial_substeps)

        p_tau_T = self.parent['tau_T']
        f_tau_T = self.fragment['tau_T']
        p_tau_d = self.parent['tau_d']

        # Compute this once ahead of time
        # More factors to fill out integral similar to eq. (36) Festou 1981
        integration_factor = (
            (1 / (4 * np.pi * p_tau_d)) *
            self.vmodel['d_alpha'] / (4.0 * np.pi)
        )

        # Calculate the density contributions over the volume of the comet
        # atmosphere due to one ejection axis
        for j in range(0, self.angular_points):
            cur_angle = self.vmodel['angular_grid'][j]
            # Loop through the radial points along this axis
            for i in range(0, self.radial_points):

                cur_r = self.vmodel['fast_radial_grid'][i]
                x = cur_r * np.sin(cur_angle)
                y = cur_r * np.cos(cur_angle)

                # Decide how granular our epsilon should be
                d_epsilon_steps = len(ejection_radii)
                d_epsilon = (
                            (self.vmodel['epsilon_max'] - cur_angle)
                    / d_epsilon_steps
                )

                # Maximum radius that contributes to point x,y when there is a
                # a max ejection angle
                if(self.vmodel['epsilon_max'] < np.pi):
                    r_limit = y - (x / np.tan(self.vmodel['epsilon_max']))
                # Set the last element to be r_coma or the above limit
                ejection_radii[d_epsilon_steps - 1] = r_limit

                # We already filled out the very last element in the array
                # above, so it's d_epsilon_steps - 1
                for dE in range(0, d_epsilon_steps - 1):
                    ejection_radii[dE] = (
                        y -
                        x / np.tan((dE + 1) * d_epsilon + cur_angle)
                    )

                ejection_radii_start = 0
                # Number of slices along the contributing axis for each step
                num_slices = self.angular_substeps

                # Loop over radial chunk that contributes to x,y
                for ejection_radii_end in ejection_radii:

                    # We are slicing up this axis into pieces
                    dr = (
                         (ejection_radii_end - ejection_radii_start) /
                        num_slices
                    )

                    # Loop over tiny slices along this chunk
                    for m in range(0, num_slices):

                        # TODO: We could probably eliminate m by making a
                        # linear space from ejection_radii_start to
                        # ejection_radii_end

                        # Current distance along contributing axis
                        cur_r = (m + 0.5) * dr + ejection_radii_start
                        # This is the distance from the NP axis point to the
                        # current point on the ray, squared
                        sep_dist = np.sqrt(x * x + (cur_r - y) * (cur_r - y))

                        cos_eject = (y - cur_r) / sep_dist
                        sin_eject = x / sep_dist

                        # Calculate sqrt(vR^2 - u^2 sin^2 gamma)
                        v_factor = np.sqrt(
                            vf * vf - (vp * vp) * sin_eject ** 2)

                        # The first (and largest) of the two solutions for the
                        # velocity when it arrives
                        v_one = vp * cos_eject + v_factor

                        # Time taken to travel from the dissociation point at
                        # v1, reject if the time is too large (and all
                        # fragments have decayed)
                        t_frag_one = sep_dist / v_one
                        if t_frag_one > time_limit:
                            continue

                        # This is the total time between parent emission from
                        # nucleus and fragment arriving at our point of
                        # interest, which we then use to look up Q at that time
                        # in the past
                        t_total_one = (cur_r / vp) + t_frag_one

                        # Division by parent velocity makes this production per
                        # unit distance for radial integration q(r, epsilon)
                        # given by eq. 32, Festou 1981
                        prod_one = self.production_at_time(t_total_one) / vp
                        q_r_eps_one = (
                            (v_one * v_one * prod_one) /
                            (vf * np.abs(v_one - vp * cos_eject))
                        )

                        # Parent extinction when traveling along to the
                        # dissociation site
                        p_extinction = np.e ** (-cur_r / (p_tau_T * vp))
                        # Fragment extinction when traveling at speed v1
                        f_extinction_one = np.e ** (-t_frag_one / f_tau_T)

                        # First differential addition to the density
                        # integrating along dr, similar to eq. (36) Festou
                        # 1981, due to the first velocity
                        n_r_one = (
                            (p_extinction * f_extinction_one * q_r_eps_one) /
                            (sep_dist ** 2 * v_one)
                        )

                        # Add this contribution to the density grid
                        self.vmodel['density_grid'][i][j] = (
                            self.vmodel['density_grid'][i][j] +
                            n_r_one * dr
                        )

                        # Check if there is a second solution for the velocity
                        if vf > vp:
                            continue

                        # Compute the contribution from the second solution for
                        # v in the same way
                        v_two = vp * cos_eject - v_factor
                        t_frag_two = sep_dist / v_two
                        if t_frag_two > time_limit:
                            continue
                        t_total_two = (cur_r / vp) + t_frag_two
                        prod_two = self.production_at_time(t_total_two) / vp
                        q_r_eps_two = (
                            (v_two * v_two * prod_two) /
                            (vf * np.abs(v_two - vp * cos_eject))
                        )
                        f_extinction_two = np.e ** (-t_frag_two / f_tau_T)
                        n_r_two = (
                            (p_extinction * f_extinction_two * q_r_eps_two) /
                            (v_two * sep_dist ** 2)
                        )
                        self.vmodel['density_grid'][i][j] = (
                            self.vmodel['density_grid'][i][j] +
                            n_r_two * dr
                        )

                    # Next starting radial point is the current end point
                    ejection_radii_start = ejection_radii_end

            if(self.print_progress is True):
                progress_percent = (j + 1) * 100 / self.angular_points
                print(f'Computing: {progress_percent:3.1f} %', end='\r')

        # Loops automatically over the 2d grid
        self.vmodel['density_grid'] *= integration_factor
        # phew

        """
            Performs angular part of the integration to yield density in m^-3
            as a function of radius.  Assumes spherical symmetry of parent
            production.

            Fills vmodel['radial_density'] and vmodel['fast_radial_density']
            with and without units respectively
            Fills vmodel['r_dens_interpolation'] with cubic spline
            interpolation of the radial density, which takes radial coordinate
            in m and outputs the density at that coord in m^-3, without units
            attached
        """

        # Make array to hold our data, no units
        self.vmodel['fast_radial_density'] = np.zeros(self.radial_points)

        # loop through grid array
        for i in range(0, self.radial_points):
            for j in range(0, self.angular_points):
                # Current angle is theta
                theta = self.vmodel['angular_grid'][j]
                # Integration factors from angular part of integral, similar to
                # eq. (36) Festou 1981
                dens_contribution = (
                    2.0 * np.pi * np.sin(theta) *
                    self.vmodel['density_grid'][i][j]
                )
                self.vmodel['fast_radial_density'][i] += dens_contribution

        # Tag with proper units
        self.vmodel['radial_density'] = (
            self.vmodel['fast_radial_density'] / (u.m ** 3)
        )

    def _interpolate_radial_density(self):
        """ Interpolate the radial density.

        Takes our fragment density grid and constructs density as a function of
        arbitrary radius
        """

        if not scipy:
            raise RequiredPackageUnavailable('scipy')
        # Interpolate this radial density grid with a cubic spline for lookup
        # at non-grid radii, input in m, output in 1/m^3
        self.vmodel['r_dens_interpolation'] = (
            CubicSpline(self.vmodel['fast_radial_grid'],
                        self.vmodel['fast_radial_density'],
                        bc_type='natural')
        )

    def _column_density_at_rho(self, rho):
        """ Calculate the column density of fragments at a given impact parameter

            Parameters
            ----------
            rho : float
                Impact parameter of the column density integration, in meters

            Returns
            -------
            float
                Column density at the given impact parameter in m^-2, no
                astropy units attached

            Notes
            -----
            We return zero column density beyond the edge of our grid, so if
            there is still significant column density near the edge of the grid
            this can lead to strange graphing results and sharp cutoffs.
        """

        r_max = self.vmodel['max_grid_radius'].value
        if(rho > r_max):
            return 0
        rhosq = rho ** 2
        z_max = np.sqrt(r_max ** 2 - rhosq)

        def column_density_integrand(z):
            return self.vmodel['r_dens_interpolation'](np.sqrt(z ** 2 + rhosq))

        # Romberg is significantly slower for impact parameters near the
        # nucleus, and becomes much faster at roughly 60 times the collision
        # sphere radius, after a few tests
        # The results of both were the same to within .1% or better, generally
        # TODO: test this for a range of productions to see if this holds

        if rho < (60 * self.vmodel['collision_sphere_radius'].value):
            c_dens = (
                quad(column_density_integrand,
                     -z_max, z_max, limit=3000)
            )[0]
        else:
            c_dens = (
                2 * romberg(column_density_integrand,
                            0, z_max, rtol=0.0001, divmax=20)
            )

        # result is in 1/m^2
        return c_dens

    def _compute_column_density(self):
        """ Compute the column density on the grid and interpolate the results

            Computes the fragment column density on a grid and interpolates for
            fragment column density as a function of arbitrary radius.

            Notes
            -----
            The interpolator returns column density in m^-2, no astropy units
            attached
        """
        if not scipy:
            raise RequiredPackageUnavailable('scipy')

        # make a grid for the column density values
        column_density_grid = self._make_radial_logspace_grid()
        # vectorize so we can hand the whole grid to our column density
        # computing function in one line
        cd_vectorized = np.vectorize(self._column_density_at_rho)
        # array now holds corresponding column density values
        column_densities = cd_vectorized(column_density_grid)

        self.vmodel['fast_column_density_grid'] = column_density_grid
        self.vmodel['column_density_grid'] = column_density_grid * u.m
        self.vmodel['column_densities'] = column_densities / (u.m ** 2)
        # Interpolation gives column density in m^-2
        self.vmodel['column_density_interpolation'] = (
            CubicSpline(
                column_density_grid,
                column_densities,
                bc_type='natural'
            )
        )

    def calc_num_fragments_theory(self):
        """ The total number of fragment species we expect in the coma

            Returns
            -------
            float
                Total number of fragment species we expect in the coma
                theoretically

            Notes
            -----
            This needs to be rewritten to better handle time dependence; the
            original implementation was designed for abrupt but small parent
            production changes.
        """

        # TODO: Re-implement this as differential equations with parent
        # photodissociation as a source of fragments, and fragment
        # photodissociation as a sink
        vp = self.parent['v_outflow']
        vf = self.fragment['v_photo']
        p_tau_T = self.parent['tau_T']
        f_tau_T = self.fragment['tau_T']
        p_tau_d = self.parent['tau_d']
        t_perm = self.vmodel['t_perm_flow'].to(u.s).value

        alpha = f_tau_T * p_tau_T / p_tau_d

        max_r = self.vmodel['max_grid_radius'].value
        last_density_element = len(self.vmodel['fast_radial_density']) - 1
        edge_adjust = (
            (np.pi * max_r * max_r * (vf + vp) *
             self.vmodel['fast_radial_density'][last_density_element])
        )

        num_time_slices = 10000
        # array starting at tperm (farthest back in time we expect something to
        # hang around), stepped down to zero (when the observation takes place)
        time_slices = np.linspace(t_perm, 0, num_time_slices, endpoint=False)

        theory_total = 0
        for i, t in enumerate(time_slices[:-1], start=0):
            if t > t_perm:
                t1 = t_perm
            else:
                t1 = t
            t1 /= p_tau_T

            if i != (num_time_slices - 1):
                t2 = time_slices[i + 1]
            else:
                t2 = 0
            t2 /= p_tau_T
            mult_factor = -np.e ** (-t1) + np.e ** (-t2)
            theory_total += self.production_at_time(t) * mult_factor

        theory_total *= alpha
        theory_total -= edge_adjust

        return theory_total

    def calc_num_fragments_grid(self):
        """ Total number of fragments in the coma.

            Calculates the total number of fragment species by integrating the
            density grid over its volume

            Returns
            -------
            float
                Number of fragments in the coma based on our grid calculations

            Notes
            -----
            Outbursts/time dependent production in general will make this
            result poor due to the grid being sized to capture a certain
            fraction of parents/fragments at the oldest (first) production
            rate.  The farther you get from this base production, the farther
            the model will deviate from capturing the requested percentage of
            particles.
        """
        max_r = self.vmodel['max_grid_radius'].value

        def vol_integrand(r, r_func):
            return (r_func(r) * r ** 2)

        r_int = romberg(
            vol_integrand, 0, max_r,
            args=(self.vmodel['r_dens_interpolation'], ),
            rtol=0.0001, divmax=20)
        return 4 * np.pi * r_int

    def _column_density(self, rho):
        """ Gives fragment column density at arbitrary impact parameter

            Parameters
            ----------
            rho : float
                Impact parameter, in meters, no astropy units attached

            Returns
            -------
            float
                Fragment column density at given impact parameter, in m^-2
        """
        return self.vmodel['column_density_interpolation'](rho)

    def _volume_density(self, r):
        """ Gives fragment volume density at arbitrary radius

            Parameters
            ----------
            r : float
                Distance from nucles, in meters, no astropy units attached

            Returns
            -------
            float
                Fragment volume density at specified radius
        """
        return self.vmodel['r_dens_interpolation'](r)
