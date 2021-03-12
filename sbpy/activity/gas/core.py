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

import copy
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
from ...exceptions import RequiredPackageUnavailable
from .. core import (RectangularAperture, GaussianAperture,
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

        return self._volume_density(r.to_value('m')) / u.m**3

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

        rho = rho.to_value('m', equiv)
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
            `~scipy.integrate.quad`.


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
        points = rho * np.logspace(-4, 4, limit // 2)
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
            `~scipy.integrate.quad` (circular, annular, Gaussian) and
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
                return (rho * np.exp(-rho**2 / sigma**2 / 2)
                        * self._column_density(rho))

            sigma = aper.sigma.to_value('m')
            N, err = quad(f, 0, np.inf, args=(sigma,), epsabs=epsabs)
            N *= 2 * np.pi
            err *= 2 * np.pi

        return N, err


class Haser(GasComa):
    """Haser coma model.

    Some functions require ``scipy``.


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
        n = (self.Q / self.v).to_value('1/m') / r**2 / 4 / np.pi
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
        x = (rho / parent).to_value(u.dimensionless_unscaled)
        N = (self.Q * rho / self.v).to_value(u.dimensionless_unscaled)
        if self.daughter is None or self.daughter == 0:
            N *= 1 / x - self._K1(x) + np.pi / 2 - self._iK0(x)
        else:
            daughter = self.daughter.to(rho.unit)
            y = (rho / daughter).to_value('')
            N *= ((daughter / (parent - daughter)).to_value('')
                  * (self._iK0(y) - self._iK0(x) + x**-1 - y**-1
                     + self._K1(y) - self._K1(x)))

        return N
    total_number.__doc__ = GasComa.total_number.__doc__


class VectorialModel(GasComa):

    """ NOTE: We allow Q to vary with time, so I suppose we pass along Q to base class without using it """
    # TODO: We would have to rewrite the time-varying production stuff to accept Q here and then let the param structure
    # define some variation of this initial production rate
    # TODO: these inputParameters could be moved to a different data structure like an astropy QTable with named fields
    @u.quantity_input(Q=(u.s**-1, u.mol / u.s), v=u.m / u.s)
    # def __init__(self, Q, v, vModelParams):
    def __init__(self, Q, v, inputParameters):
        super().__init__(Q, v)

        # Deep copy so we can make changes without affecting data user passed in
        self.vModelParams = copy.deepcopy(inputParameters)

        self.vModelParams['ParentSpecies']['Velocity'] = v.to(u.m/u.s).value
        # fill a few short names for readability later
        self.par = self.vModelParams['ParentSpecies']
        self.frag = self.vModelParams['FragmentSpecies']
        self.comet = self.vModelParams['Comet']

        self.comet['TimeAtProductions'] = self.comet['TimeAtProductions'].to(u.s).value
        self.numProductionRates = len(self.comet['ProductionRates'])

        self.par['TotalLifetime'] = self.par['TotalLifetime'].to(u.s).value
        self.par['DissociativeLifetime'] = self.par['DissociativeLifetime'].to(u.s).value

        self.frag['Velocity'] = self.frag['Velocity'].to(u.m/u.s).value
        self.frag['TotalLifetime'] = self.frag['TotalLifetime'].to(u.s).value

        self.sigma = self.par['Sigma'].to(u.m**2).value

        """Initialize data structures to hold our calculations"""
        self.vModel = {}

        # Calculate up a few things
        self._setupCalculations()

        # Build the radial grid
        self.vModel['FastRadialGrid'] = self._makeLogspaceGrid()
        self.vModel['RadialGrid'] = self.vModel['FastRadialGrid']*(u.m)
        self.vModel['NumRadialGridpoints'] = len(self.vModel['FastRadialGrid'])

        # Angular grid
        self.vModel['NumAngularGridpoints'] = self.vModelParams['Grid']['AngularGridpoints']
        self.vModel['dAlpha'] = self.vModel['EpsilonMax']/self.vModel['NumAngularGridpoints']
        # Make array of angles adjusted up away from zero, to keep from calculating a radial line's contribution
        # to itself
        self.vModel['AngularAlphaGrid'] = np.linspace(0, self.vModel['EpsilonMax'], num=self.vModel['NumAngularGridpoints'], endpoint=False)
        # This maps addition over the whole array automatically
        self.vModel['AngularAlphaGrid'] += self.vModel['dAlpha']/2

        # makes a 2d array full of zero values
        self.vModel['DensityGrid'] = np.zeros((self.vModel['NumRadialGridpoints'],
                                               self.vModel['NumAngularGridpoints']))

        self._computeDensityGrid()
        self._interpolateColumnDensity()

    def _setupCalculations(self):
        """ A few calculations before the heavy lifting starts """

        """ Calculate collision sphere radius based on oldest production rate, Eq. (5) in Festou 1981

            Note that this is only calculated with the first production rate, because it is assumed that the first
            production rate has had enough time to reach a steady state before letting production vary with time.
            We also use that vtherm = 0.25 * outflow velocity
        """
        # this factor comes from molecular flux of ideal gas moving through a surface, in our case the surface of the
        # collision sphere
        vtherm = self.par['Velocity']*0.25
        q = self.comet['ProductionRates'][0]
        vp = self.par['Velocity']
        vf = self.frag['Velocity']

        self.vModel['MiscOutput'] = {}

        # Eq. 5 of Festou 1981
        self.vModel['MiscOutput']['CollisionSphereRadius'] = ((self.sigma * q * vtherm)/(vp * vp))*u.m

        """ Calculates the radius of the coma given our input parameters """
        # NOTE: Equation (16) of Festou 1981 where alpha is the percent destruction of molecules
        parentBetaR = -np.log(1.0 - self.par['DestructionLevel'])
        parentR = parentBetaR * vp * self.par['TotalLifetime']
        fragR = vp * self.comet['TimeAtProductions'][0]
        self.vModel['MiscOutput']['ComaRadius'] = min(parentR, fragR)*u.m

        """ Calculates the time needed to hit a steady, permanent production """
        fragmentBetaR = -np.log(1.0 - self.frag['DestructionLevel'])
        # Permanent flow regime
        permFlowR = self.vModel['MiscOutput']['ComaRadius'].value + ((vp + vf) * fragmentBetaR * self.frag['TotalLifetime'])

        timeInSecs = self.vModel['MiscOutput']['ComaRadius'].value/vp + (permFlowR - self.vModel['MiscOutput']['ComaRadius'].value)/(vp + vf)
        self.vModel['MiscOutput']['TimeToPermanentFlowRegime'] = (timeInSecs * u.s).to(u.day)

        """ this is the total radial size that parents & fragments occupy, beyond which we assume zero density
        """
        # Calculate the lesser of the radii of two situations
        #       Permanent flow regime
        permFlowR = self.vModel['MiscOutput']['ComaRadius'].value + ((vp + vf) * fragmentBetaR * self.frag['TotalLifetime'])
        #       or Outburst situation
        outburstR = (vp + vf) * self.comet['TimeAtProductions'][0]
        self.vModel['MaxRadiusOfGrid'] = min(permFlowR, outburstR)*u.m

        # Two cases for angular range of ejection of fragment based on relative velocities of parent and fragment species
        if(vf < vp):
            self.vModel['EpsilonMax'] = np.arcsin(vf/vp)
        else:
            self.vModel['EpsilonMax'] = np.pi

    def productionRateAtTime(self, t):
        """ Get production rate at time t in s, with positive values representing the past, so
            productionRateAtTime(5) will provide the production 5 seconds ago

            For times in the past beyond the specified production rates, returns zero
        """
        if t > self.comet['TimeAtProductions'][0]:
            return 0.0

        for i in range(0, self.numProductionRates):
            binStartTime = self.comet['TimeAtProductions'][i]
            if i == (self.numProductionRates - 1):
                # We're at the end of the array, so stop time is zero seconds ago
                binStopTime = 0
            else:
                # Otherwise we go all the way to the start of the next one
                binStopTime = self.comet['TimeAtProductions'][i+1]

            # NOTE: remember that these times are in seconds ago, so the comparison is backward
            if t < binStartTime and t >= binStopTime:
                return self.comet['ProductionRates'][i]

    def _makeLogspaceGrid(self):
        """ Creates a grid (in meters) with numpy's logspace function that covers the expected radial size,
            stretching from 2 times the collision sphere radius (near the nucleus be dragons) out to the calculated max
            If we get too close to the nucleus things go very badly so don't do it, dear reader
        """
        rStartpointPower = np.log10(self.vModel['MiscOutput']['CollisionSphereRadius'].value * 2)
        rEndpointPower = np.log10(self.vModel['MaxRadiusOfGrid'].value)
        return np.logspace(rStartpointPower, rEndpointPower.astype(float), num=self.vModelParams['Grid']['RadialGridpoints'], endpoint=True)

    def _computeDensityGrid(self):
        """ Computes the density at different radii and due to each ejection angle, performing the
           radial integration of eq. (36), Festou 1981 with only one fragment velocity.
           The resulting units will be in 1/(m^3) because we work in m, s, and m/s.
        """
        vp = self.par['Velocity']
        vf = self.frag['Velocity']

        # Follow fragments until they have been totally destroyed
        timeLimit = 8.0 * self.frag['TotalLifetime']
        rComa = self.vModel['MiscOutput']['ComaRadius'].value
        rLimit = rComa

        # temporary radial array for when we loop through 0 to epsilonMax
        ejectionRadii = np.zeros(self.vModelParams['Grid']['SubgridRadialSteps'])

        pTotLifetime = self.par['TotalLifetime']
        fTotLifetime = self.frag['TotalLifetime']
        pDisLifetime = self.par['DissociativeLifetime']

        # Compute this once ahead of time
        # More factors to fill out integral similar to eq. (36) Festou 1981
        IntegrationFactor = (1/(4 * np.pi * pDisLifetime)) * self.vModel['dAlpha']/(4.0 * np.pi)

        # Build the 2d elementary density table
        # Loop through alpha
        for j in range(0, self.vModel['NumAngularGridpoints']):
            curAngle = self.vModel['AngularAlphaGrid'][j]
            # Loop through the radial points along this axis
            for i in range(0, self.vModel['NumRadialGridpoints']):

                curR = self.vModel['FastRadialGrid'][i]
                x = curR * np.sin(curAngle)
                y = curR * np.cos(curAngle)

                # Decide how granular our epsilon should be
                dEpsilonSteps = len(ejectionRadii)
                dEpsilon = (self.vModel['EpsilonMax'] - curAngle)/dEpsilonSteps

                # Maximum radius that contributes to point x,y when there is a a max ejection angle
                if(self.vModel['EpsilonMax'] < np.pi):
                    rLimit = y - (x/np.tan(self.vModel['EpsilonMax']))
                # Set the last element to be rComa or the above limit
                ejectionRadii[dEpsilonSteps-1] = rLimit

                for dE in range(0, dEpsilonSteps-1):  # We already filled out the very last element in the array above
                    ejectionRadii[dE] = y - x/np.tan((dE+1)*dEpsilon + curAngle)

                ejectionRadiiStart = 0
                # Number of slices along the contributing axis for each step
                NumRadialSlices = self.vModelParams['Grid']['SubgridAngularSteps']

                # Loop over radial chunk that contributes to x,y
                for ejectionRadiiEnd in ejectionRadii:

                    # We are slicing up this axis into pieces
                    dr = (ejectionRadiiEnd - ejectionRadiiStart)/NumRadialSlices

                    # Loop over tiny slices along this chunk
                    for m in range(0, NumRadialSlices):

                        # TODO: We could probably eliminate m by making a linear space from ejectionRadiiStart to
                        # ejectionRadiiEnd

                        # Current distance along contributing axis
                        R = (m + 0.5)*dr + ejectionRadiiStart
                        # This is the distance from the NP axis point to the current point on the ray, squared
                        sepDist = np.sqrt(x * x + (R - y)*(R - y))

                        cosEjection = (y - R)/sepDist
                        sinEjection = x/sepDist

                        # Calculate sqrt(vR^2 - u^2 sin^2 gamma)
                        vFactor = np.sqrt(vf * vf - (vp * vp)*sinEjection**2)

                        # The first (and largest) of the two solutions for the velocity when it arrives
                        vOne = vp*cosEjection + vFactor

                        # Time taken to travel from the dissociation point at v1, reject if the time is too large (and all
                        # fragments have decayed)
                        tFragmentOne = sepDist/vOne
                        if tFragmentOne > timeLimit:
                            continue

                        # This is the total time between parent emission from nucleus and fragment arriving at our point of
                        # interest, which we then use to look up Q at that time in the past
                        tTotalOne = (R/vp) + tFragmentOne

                        # Division by parent velocity makes this production per unit distance for radial integration
                        # q(r, epsilon) given by eq. 32, Festou 1981
                        prodOne = self.productionRateAtTime(tTotalOne)/vp
                        qREpsOne = (vOne*vOne*prodOne)/(vf * np.abs(vOne - vp*cosEjection))

                        # Parent extinction when traveling along to the dissociation site
                        pExtinction = np.e**(-R/(pTotLifetime * vp))
                        # Fragment extinction when traveling at speed v1
                        fExtinctionOne = np.e**(-tFragmentOne/fTotLifetime)

                        # First differential addition to the density integrating along dr, similar to eq. (36) Festou 1981,
                        # due to the first velocity
                        densityOne = (pExtinction * fExtinctionOne * qREpsOne)/(sepDist**2 * vOne)

                        # Add this contribution to the density grid
                        self.vModel['DensityGrid'][i][j] = self.vModel['DensityGrid'][i][j] + densityOne*dr

                        # Check if there is a second solution for the velocity
                        if vf > vp:
                            continue

                        # Compute the contribution from the second solution for v in the same way
                        vTwo = vp*cosEjection - vFactor
                        tFragmentTwo = sepDist/vTwo
                        if tFragmentTwo > timeLimit:
                            continue
                        tTotalTwo = (R/vp) + tFragmentTwo
                        prodTwo = self.productionRateAtTime(tTotalTwo)/vp
                        qREpsTwo = (vTwo * vTwo * prodTwo)/(vf * np.abs(vTwo - vp*cosEjection))
                        fExtinctionTwo = np.e**(-tFragmentTwo/fTotLifetime)
                        densityTwo = (pExtinction * fExtinctionTwo * qREpsTwo)/(vTwo * sepDist**2)
                        self.vModel['DensityGrid'][i][j] = self.vModel['DensityGrid'][i][j] + densityTwo*dr

                    # Next starting radial point is the current end point
                    ejectionRadiiStart = ejectionRadiiEnd

            if(self.vModelParams['Misc']['PrintDensityProgress'] is True):
                progressPercent = (j+1)*100/self.vModel['NumAngularGridpoints']
                print(f'Computing: {progressPercent:3.1f} %', end='\r')

        # Loops automatically over the 2d grid
        self.vModel['DensityGrid'] *= IntegrationFactor
        # phew

        """ Performs angular part of the integration to yield density in m^-3 as a function of radius.
            Assumes spherical symmetry of parent production.

            Fills vModel['RadialDensity'] and vModel['FastRadialDensity'] with and without units respectively
            Fills vModel['rDensInterpolator'] with cubic spline interpolation of the radial density,
                which takes radial coordinate in m and outputs the density at that coord in m^-3
        """

        # Make array to hold our data, no units
        self.vModel['FastRadialDensity'] = np.zeros(self.vModel['NumRadialGridpoints'])

        # loop through grid array
        for i in range(0, self.vModel['NumRadialGridpoints']):
            for j in range(0, self.vModel['NumAngularGridpoints']):
                # Current angle is theta
                theta = self.vModel['AngularAlphaGrid'][j]
                # Integration factors from angular part of integral, similar to eq. (36) Festou 1981
                densityToAdd = 2.0 * np.pi * np.sin(theta) * self.vModel['DensityGrid'][i][j]
                self.vModel['FastRadialDensity'][i] += densityToAdd

        # Tag with proper units
        self.vModel['RadialDensity'] = self.vModel['FastRadialDensity']/(u.m**3)

        # Turn our grid into a function of r
        self._interpolateRadialDensity()

        # Count up the number of fragments in the grid versus theoretical value
        self.vModel['MiscOutput']['NumFragmentsTheory'] = self.calcNumFragmentsTheory()
        self.vModel['MiscOutput']['NumFragmentsFromGrid'] = self.calcNumFragmentsFromGrid()

    def _interpolateRadialDensity(self):
        # Interpolate this radial density grid with a cubic spline for lookup at non-grid radii, input in m, out in 1/m^3
        self.vModel['rDensInterpolator'] = CubicSpline(self.vModel['FastRadialGrid'], self.vModel['FastRadialDensity'], bc_type='natural')

    def _calculateColumnDensity(self, rho):
        """ Returns the number of fragment species per m^2 at the given impact parameter, in m """
        rMax = self.vModel['MaxRadiusOfGrid'].value
        if(rho > rMax):
            return 0
        rhosq = rho**2
        zMax = np.sqrt(rMax**2 - rhosq)

        def columnDensityIntegrand(z):
            return self.vModel['rDensInterpolator'](np.sqrt(z**2 + rhosq))

        # Romberg is significantly slower for impact parameters near the nucleus, and becomes much faster at roughly 60 times
        # the collision sphere radius, after a few tests
        # The results of both were the same to within .1% or better, generally

        if rho < (60 * self.vModel['MiscOutput']['CollisionSphereRadius'].value):
            cDens = (quad(columnDensityIntegrand, -zMax, zMax, limit=1000))[0]
        else:
            cDens = 2 * romberg(columnDensityIntegrand, 0, zMax, rtol=0.0001, divmax=20)

        # result is in 1/m^2
        return cDens

    def _interpolateColumnDensity(self):
        """ computes the column density on a grid and produces an interpolation function based on it.
            The interpolator takes input in m from nucleus and return column density in m^-2
        """
        cDensGrid = self._makeLogspaceGrid()
        cdVec = np.vectorize(self._calculateColumnDensity)
        columnDensities = cdVec(cDensGrid)

        self.vModel['ColumnDensity'] = {}
        self.vModel['ColumnDensity']['FastCDGrid'] = cDensGrid
        self.vModel['ColumnDensity']['CDGrid'] = cDensGrid*u.m
        self.vModel['ColumnDensity']['Values'] = columnDensities/(u.m**2)
        # Interpolator gives column density in m^-2
        self.vModel['ColumnDensity']['Interpolator'] = CubicSpline(cDensGrid, columnDensities, bc_type='natural')

    def calcNumFragmentsTheory(self):
        """ Returns the total number of fragment species we expect in the coma theoretically """
        vp = self.par['Velocity']
        vf = self.frag['Velocity']
        pTotLifetime = self.par['TotalLifetime']
        fTotLifetime = self.frag['TotalLifetime']
        pDisLifetime = self.par['DissociativeLifetime']
        pRates = self.comet['ProductionRates']
        pTimes = self.comet['TimeAtProductions']
        tPerm = self.vModel['MiscOutput']['TimeToPermanentFlowRegime'].to(u.s).value

        mR = self.vModel['MaxRadiusOfGrid'].value
        lastDensityElement = len(self.vModel['FastRadialDensity'])-1

        theoryTot = 0
        for i in range(0, len(pRates)):
            if(pTimes[i] > tPerm):
                t1 = tPerm/pTotLifetime
            else:
                t1 = pTimes[i]/pTotLifetime

            if i != (self.numProductionRates - 1):
                t2 = pTimes[i+1]/pTotLifetime
            else:
                t2 = 0
            theoryTot += pRates[i]*(-np.e**(-t1) + np.e**(-t2))

        theoryTot = theoryTot*(fTotLifetime*pTotLifetime/pDisLifetime) \
            - (np.pi * mR * mR * (vf + vp) * self.vModel['FastRadialDensity'][lastDensityElement])

        return theoryTot

    def calcNumFragmentsFromGrid(self):
        """ Returns the total number of fragment species by integrating the density grid over the volume """
        maxR = self.vModel['MaxRadiusOfGrid'].value

        def volIntegrand(r, rFunc):
            return (rFunc(r) * r**2)

        # Maybe stay away from r = 0 here, but spline type of 'natural' seems to handle origin well
        rInt = romberg(volIntegrand, 0, maxR, args=(self.vModel['rDensInterpolator'], ),
                       rtol=0.0001, divmax=20)
        return 4*np.pi*rInt

    def _column_density(self, rho):
        """ rho in meters """
        return self.vModel['ColumnDensity']['Interpolator'](rho)

    def _volume_density(self, rho):
        """ rho in meters """
        return self.vModel['rDensInterpolator'](rho)
