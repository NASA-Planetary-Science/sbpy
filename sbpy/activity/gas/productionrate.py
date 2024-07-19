# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
===========================
sbpy Production Rate Module
===========================

:author: Giannina Guzman (gguzman2@villanova.edu)

created on June 26, 2019
"""

import tempfile
import numpy as np
import astropy.constants as con
import astropy.units as u

try:
    from astroquery.jplspec import JPLSpec
    from astroquery.lamda import Lamda
except ImportError:
    pass

try:
    import pyradex
except ImportError:
    pyradex = None

from ...bib import register
from ...data import Phys
from ...utils.decorators import requires
from ...utils import required_packages

__all__ = ['LTE', 'NonLTE', 'einstein_coeff',
           'intensity_conversion', 'beta_factor', 'total_number',
           'from_Haser']

__doctest_requires__ = {
    "LTE.from_Drahus": ["astroquery>=0.4.7"],
    "NonLTE.from_pyradex": ["astroquery>=0.4.7", "pyradex"],
    "from_Haser": ["astroquery>=0.4.7"],
}


def intensity_conversion(mol_data):
    """
    Returns conversion of the integrated line intensity at 300 K
    (from the JPL molecular spectra catalog) to a chosen temperature

    Parameters
    ----------
    mol_data : `~sbpy.data.Phys`
        `~sbpy.data.Phys` object that contains the following data,
        using `~astropy.units` for the required units:

            * Transition frequency in MHz
            * Temperature in Kelvins
            * Integrated line intensity at 300 K in MHz * nm**2
            * Partition function at 300 K (Dimensionless)
            * Partition function at designated temperature (Dimensionless)
            * Upper state degeneracy (Dimensionless)
            * Upper level energy in Joules
            * Lower level energy in Joules
            * Degrees of freedom (Dimensionless)

        Keywords that can be used for these values are found under
        `~sbpy.data.Conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.phys.from_jplspec` to obtain these values in order to
        maintain consistency and because all calculations can be handled within
        `sbpy` scope if JPLSpec is used. Yet, if you wish to use your own
        molecular data, it is possible. Make sure to inform yourself on the
        values needed for each function, their units, and their interchangeable
        keywords as part of the `~sbpy.data.Phys` data class.

    Returns
    -------
    intl : `~astropy.Quantity`
        Integrated line intensity at designated temperature in MHz * nm**2,
        which can be appended to the original `sbpy.data.Phys` object for
        future calculations

    References
    ----------
    Picket et al 1998, JQSRT 60, 883-890

    """

    if not isinstance(mol_data, Phys):
        raise ValueError('mol_data must be a `sbpy.data.Phys` instance.')

    temp = mol_data['Temperature'][0]
    lgint = mol_data['lgint300'][0]
    part300 = mol_data['partfn300'][0]
    partition = mol_data['partfn'][0]
    eup_J = mol_data['eup_j'][0]
    elo_J = mol_data['elo_J'][0]
    df = mol_data['degfr'][0]

    register(intensity_conversion, {'conversion': '1998JQSRT..60..883P'})

    k = con.k_B.to('J/K')  # Boltzmann constant

    if (eup_J - elo_J) < (k * min(temp, 300 * u.K)):

        if df in (0, 2):
            n = 1
        else:
            n = 3./2
        intl = lgint*(300*u.K/temp)**(n+1)*np.exp(-(1/temp - 1/(300*u.K))
                                                  * elo_J/k)

    else:

        intl = lgint*(part300/partition)*(np.exp(-elo_J/(k*temp)) -
                                          np.exp(-eup_J/(k*temp))) / \
            (np.exp(-elo_J/(k*300 * u.K)) - np.exp(-eup_J/(k*300*u.K)))

    return intl


def einstein_coeff(mol_data):
    """
    Einstein coefficient from molecular data

    Parameters
    ----------
    mol_data : `~sbpy.data.phys`
        `~sbpy.data.phys` object that contains the following data,
        using `astropy.units` for the required units:

            * Transition frequency in MHz
            * Temperature in Kelvins
            * Integrated line intensity at 300 K in MHz * nm**2
            * Integrated line intensity at desired temperature
            * Partition function at 300 K (Dimensionless)
            * Partition function at designated temperature (Dimensionless)
            * Upper state degeneracy (Dimensionless)
            * Upper level energy in Joules
            * Lower level energy in Joules
            * Degrees of freedom (Dimensionless)

        Keywords that can be used for these values are found under
        `~sbpy.data.Conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.phys.from_jplspec` to obtain these values in order to
        maintain consistency. Yet, if you wish to use your own molecular data,
        it is possible. Make sure to inform yourself on the values needed for
        each function, their units, and their interchangeable keywords as part
        of the data class.

    Returns
    -------
    einstein_coeff : `~astropy.Quantity`
        Spontaneous emission coefficient (1/s), which can be appended
        to the original `sbpy.phys` object for future calculations

    """

    if not isinstance(mol_data, Phys):
        raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

    temp = mol_data['Temperature'][0]
    lgint = mol_data['lgint300'][0]
    part300 = mol_data['partfn300'][0]
    partition = mol_data['partfn'][0]
    eup_J = mol_data['eup_j'][0]
    elo_J = mol_data['elo_J'][0]
    df = mol_data['degfr'][0]
    t_freq = mol_data['t_freq'][0]
    gu = mol_data['dgup'][0]

    h = con.h.to('J*s')  # Planck constant

    k = con.k_B.to('J/K')  # Boltzmann constant

    intl = mol_data['lgint'][0]

    if (h*t_freq/(k*temp)).decompose().value and \
            (h*t_freq/(k*300*u.K)).decompose().value < 1:

        au = (lgint*t_freq
              * (part300/gu)*np.exp(eup_J / (k*300*u.K))*(1.748e-9)).value

    else:

        au = (intl*(t_freq)**2 *
              (partition/gu)*(np.exp(-(elo_J/(k*temp)).value) -
                              np.exp(-(eup_J/(k*temp)).value))**(-1)
              * (2.7964e-16)).value

    au = au / u.s

    return au


def beta_factor(mol_data, ephemobj):
    """
    Returns beta factor based on timescales from `~sbpy.activity.gas`
    and distance from the Sun using an `~sbpy.data.ephem` object.
    The calculation is:
    parent photodissociation timescale * (distance from comet to Sun)**2
    If you wish to provide your own beta factor, you can calculate the equation
    expressed in units of AU**2 * s , all that is needed is the timescale
    of the molecule and the distance of the comet from the Sun. Once you
    have the beta factor you can append it to your mol_data phys object
    with the name 'beta' or any of its alternative names.

    Parameters
    ----------

    mol_data : `~sbpy.data.phys`
        `sbpy.data.phys` object that contains AT LEAST the following data:
                | mol_tag: Molecular identifier (`int` or `str`)

        This field can be given by the user directly or found using
        `~sbpy.data.phys.from_jplspec`. If the mol_tag is an integer, the
        program will assume it is the JPL Spectral Molecular Catalog identifier
        of the molecule and will treat it as such. If mol_tag is a string,
        then it will be assumed to be the human-readable name of the molecule.
        The molecule MUST be defined in `sbpy.activity.gas.timescale`,
        otherwise this function cannot be used and the beta factor
        will have to be provided by the user directly for calculations. The
        user can obtain the beta factor from the formula provided above.
        Keywords that can be used for these values are found under
        `~sbpy.data.Conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.Phys.from_jplspec` to obtain
        these values in order to maintain consistency. Yet, if you wish to
        use your own molecular data, it is possible. Make sure to inform
        yourself on the values needed for each function, their units, and
        their interchangeable keywords as part of the Phys data class.

    ephemobj : `~sbpy.data.ephem`
        `sbpy.data.ephem` object holding ephemeride information including
        distance from comet to Sun ['r'] and from comet to observer ['delta']


    Returns
    -------
    q : `~astropy.units.Quantity`
        Beta factor 'beta', which can be appended
        to the original `sbpy.phys` object for future calculations

    """
    # imported here to avoid circular dependency with activity.gas
    from .core import photo_timescale
    from ...data import Ephem

    if not isinstance(ephemobj, Ephem):
        raise ValueError('ephemobj must be a `sbpy.data.ephem` instance.')
    if not isinstance(mol_data, Phys):
        raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

    orb = ephemobj
    delta = (orb['delta']).to('m')
    r = (orb['r'])

    if not isinstance(mol_data['mol_tag'][0], str):
        required_packages(
            "astroquery", message=f"mol_tag = {mol_data['mol_tag'][0]} requires astroquery")

        cat = JPLSpec.get_species_table()
        mol = cat[cat['TAG'] == mol_data['mol_tag'][0]]
        name = mol['NAME'].data[0]

    else:
        name = mol_data['mol_tag'][0]

    timescale = photo_timescale(name)

    if timescale.ndim != 0:
        # array
        timescale = timescale[0]

    beta = (timescale) * r**2

    return beta


def total_number(mol_data, aper, b):
    """
    Equation relating number of molecules with column density, the aperture,
    and geometry given and accounting for photodissociation, derived from data
    provided.  Feel free to use your own total number to calculate production
    rate or use this function with your own molecular data as long as you are
    aware of the needed data.

    Parameters
    ----------
    integrated_line : `~astropy.units.Quantity`
        Integrated flux of emission line.

    mol_data : `sbpy.data.phys`
        `sbpy.data.phys` object that contains AT LEAST the following data:
                | Transition frequency in MHz
                | Einstein Coefficient (1/s)
                | Beta: (Timescale (in s) * r^2 (in au))
                | Column Density in 1/m^2

        The values above can either be given by the user or obtained from the
        functions `~sbpy.activity.gas.productionrate.einstein_coeff` and
        `~sbpy.activity.gas.productionrate.beta_factor`
        Keywords that can be used for these values are found under
        `~sbpy.data.Conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.phys.from_jplspec` to obtain
        these values in order to maintain consistency. Yet, if you wish to
        use your own molecular data, it is possible. Make sure to inform
        yourself on the values needed for each function, their units, and
        their interchangeable keywords as part of the Phys data class.

    Returns
    -------
    total_number: `astropy.units.Quantity`
        An astropy Quantity containing the total number of molecules within the
        aperture (Dimensionless)

    """

    if not isinstance(mol_data, Phys):
        raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

    beta = mol_data['beta'][0]

    sigma = (1./2. * beta * b * con.c / (mol_data['t_freq'][0] * aper)).value

    tnumber = mol_data['cdensity'][0].decompose() * sigma * u.m**2 / \
        np.sqrt(np.log(2))

    return tnumber


def from_Haser(coma, mol_data, aper=25 * u.m):
    """
    Calculate production rate for `GasComa`

    Parameters
    ----------
    coma : `sbpy.activity.gas.GasComa`
        Gas coma model for ratio calculation of production rate, the
        production rate `Q` that the gas coma model expects should be an
        educated first guess. A good way to get this guess would be to
        use the function `from_drahus` found under `sbpy.activity.gas.LTE`
        The molecule name used for the `parent` argument of the coma model
        should be the same name or equivalent JPLSpec identifier used to
        calculate the total number of molecules.

    mol_data: `sbpy.data.phys`
        `sbpy.data.phys` object that contains AT LEAST the following data:

                | Total Number of Molecules (See
                | `~sbpy.activity.gas.total_number` for a calculation
                | of this datum if you don't wish to provide it yourself)

        This field can be given by the user directly or calculated using the
        necessary combinations of the following functions:
        `~sbpy.data.phys.from_jplspec`,
        `~sbpy.activity.gas.productionrate.einstein_coeff`,
        `~sbpy.activity.gas.productionrate.beta_factor`, and
        `~sbpy.activity.gas.productionrate.total_number`.
        Keywords that can be used for these values are found under
        `~sbpy.data.Conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.Phys.from_jplspec` to obtain
        these values in order to maintain consistency. Yet, if you wish to
        use your own molecular data, it is possible. Make sure to inform
        yourself on the values needed for each function, their units, and
        their interchangeable keywords as part of the Phys data class.

    aper : `~astropy.units.Quantity`
        Telescope aperture in meters. Default is 25 m

    Returns
    -------
    Q : `~astropy.units.Quantity`
        production rate

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> from sbpy.data import Ephem, Phys
    >>> from sbpy.activity import (Haser, LTE, photo_timescale, einstein_coeff,
    ...                            from_Haser)
    >>> from sbpy.activity import (intensity_conversion, beta_factor,
    ...                            total_number)

    >>> aper = 10 * u.m
    >>> mol_tag = 28001
    >>> temp_estimate = 25. * u.K
    >>> target = 'C/2016 R2'
    >>> b = 0.74
    >>> vgas = 0.5 * u.km / u.s
    >>> transition_freq = (230.53799 * u.GHz).to('MHz')
    >>> integrated_flux = 0.26 * u.K * u.km / u.s

    >>> time = Time('2017-12-22 05:24:20', format = 'iso')
    >>> ephemobj = Ephem.from_horizons(target,
    ...                                epochs=time) # doctest: +REMOTE_DATA

    >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq,
    ...                              mol_tag) # doctest: +REMOTE_DATA

    >>> intl = intensity_conversion(mol_data) # doctest: +REMOTE_DATA
    >>> mol_data.apply([intl.value] * intl.unit,
    ...                name='intl') # doctest: +REMOTE_DATA

    >>> au = einstein_coeff(mol_data) # doctest: +REMOTE_DATA
    >>> mol_data.apply([au.value] * au.unit,
    ...                name='eincoeff') # doctest: +REMOTE_DATA

    >>> beta = beta_factor(mol_data, ephemobj) # doctest: +REMOTE_DATA
    >>> mol_data.apply([beta.value] * beta.unit,
    ...                name='beta') # doctest: +REMOTE_DATA

    >>> lte = LTE()

    >>> cdensity = lte.cdensity_Bockelee(integrated_flux,
    ...                                  mol_data) # doctest: +REMOTE_DATA
    >>> mol_data.apply([cdensity.value] * cdensity.unit,
    ...                name='cdensity') # doctest: +REMOTE_DATA

    >>> tnum = total_number(mol_data, aper, b) # doctest: +REMOTE_DATA
    >>> mol_data.apply([tnum.value] * tnum.unit,
    ...                name='total_number') # doctest: +REMOTE_DATA

    >>> Q_estimate = 2.8*10**(28) / u.s
    >>> parent = photo_timescale('CO') * vgas
    >>> coma = Haser(Q_estimate, vgas, parent)

    >>> Q = from_Haser(coma, mol_data, aper=aper) # doctest: +REMOTE_DATA

    >>> Q # doctest: +REMOTE_DATA +FLOAT_CMP
        <Quantity [9.35795579e+27] 1 / s>

    References
    ----------
    Haser 1957, Bulletin de la Societe Royale des Sciences de Liege
    43, 740.
    Newburn and Johnson 1978, Icarus 35, 360-368.

    """

    from .core import GasComa

    if not isinstance(coma, GasComa):
        raise ValueError('coma must be a GasComa instance.')
    if not isinstance(mol_data, Phys):
        raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

    register('Spectroscopy', {'Total Number (eq. 15)': '2004come.book..391B'})

    # integrated_line = self.integrated_flux(transition_freq) - not yet implemented

    molecules = mol_data['total_number']

    model_molecules = coma.total_number(aper)

    Q = coma.Q * molecules/model_molecules

    return Q


class LTE():
    """ LTE Methods for calculating production_rate """

    def cdensity_Bockelee(self, integrated_flux, mol_data):
        """
        Basic equation relating column density with observed integrated flux
        without the need for an initial column density to be given.
        This is found in equation 10 in
        https://ui.adsabs.harvard.edu/abs/2004come.book..391B
        and is derived from data from JPLSpec, feel free to use your own column
        density to calculate production rate or use this function with your own
        molecular data as long as you are aware of the needed data.

        Parameters
        ----------
        integrated_flux : `~astropy.units.Quantity`
            Integrated flux of emission line.

        mol_data : `sbpy.data.phys`
            `sbpy.data.phys` object that contains AT LEAST the following data:
                    | Transition frequency in MHz
                    | Einstein Coefficient (1/s)

            This function will calculate the column
            density from Bockelee-Morvan et al. 2004 and append it to the phys
            object as 'Column Density' or any of its alternative field names.
            The values above can either be given by the user or obtained from
            the functions `~sbpy.activity.gas.productionrate.einstein_coeff`
            and `~sbpy.activity.gas.productionrate.beta_factor`
            Keywords that can be used for these values are found under
            `~sbpy.data.Conf.fieldnames` documentation. We recommend the use of
            the JPL Molecular Spectral Catalog and the use of
            `~sbpy.data.phys.from_jplspec` to obtain
            these values in order to maintain consistency. Yet, if you wish to
            use your own molecular data, it is possible. Make sure to inform
            yourself on the values needed for each function, their units, and
            their interchangeable keywords as part of the Phys data class.

        Returns
        -------
        Column Density : `astropy.units.Quantity`
            Column density from Bockelee-Morvan et al. 2004 as astropy Quantity
            (1/m^2)

        """

        if not isinstance(mol_data, Phys):
            raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

        register('Spectroscopy', {
                 'Total Number (eq. 10)': '2004come.book..391B'})

        cdensity = (8*np.pi*con.k_B*mol_data['t_freq'][0]**2 /
                    (con.h*con.c**3 * mol_data['eincoeff'][0])).decompose()
        cdensity *= integrated_flux

        return cdensity

    def from_Drahus(self, integrated_flux, mol_data, ephemobj, vgas=1 * u.km/u.s,
                    aper=25 * u.m, b=1.2):
        """
        Returns production rate based on Drahus 2012 model referenced.
        Does not include photodissociation, good for first guesses for
        more computationally intensive methods or for the Haser model
        under `sbpy.activity.gas.productionrate.from_Haser`

        Parameters
        ----------
        integrated_flux : `~astropy.units.Quantity`
            Line integral derived from spectral data in Kelvins * km/s

        mol_data : `sbpy.data.phys`
            `sbpy.data.phys` object that contains the following data:
                    | Transition frequency in MHz
                    | Temperature in Kelvins
                    | Partition function at designated temperature (unitless)
                    | Upper state degeneracy (unitless)
                    | Upper level energy in Joules
                    | Degrees of freedom (unitless)
                    | Einstein Coefficient (1/s)

            These fields can be given by the user directly or calculated using
            `~sbpy.data.phys.from_jplspec`,
            `~sbpy.activity.gas.productionrate.einstein_coeff`,
            Keywords that can be used for these values are found under
            `~sbpy.data.Conf.fieldnames` documentation. We recommend the use of
            the JPL Molecular Spectral Catalog and the use of
            `~sbpy.data.Phys.from_jplspec` to obtain
            these values in order to maintain consistency. Yet, if you wish to
            use your own molecular data, it is possible. Make sure to inform
            yourself on the values needed for each function, their units, and
            their interchangeable keywords as part of the Phys data class.

        ephemobj : `sbpy.data.ephem`
            `sbpy.data.ephem` object holding ephemeride information including
            distance from comet to Sun ['r'] and from comet to observer
            ['delta']

        vgas : `~astropy.units.Quantity`
            Gas velocity approximation in km / s. Default is 1 km / s

        aper : `~astropy.units.Quantity`
            Telescope aperture in meters. Default is 25 m

        b : int
            Dimensionless factor intrinsic to every antenna. Typical
            value, and the default for this model, is 1.22. See
            references for more information on this parameter.

        Returns
        -------
        q : `~astropy.units.Quantity`
            Production rate, not including photodissociation

        Examples
        --------
        >>> import astropy.units as u
        >>> from astropy.time import Time
        >>> from sbpy.data import Ephem, Phys
        >>> from sbpy.activity import LTE, einstein_coeff, intensity_conversion

        >>> temp_estimate = 47. * u.K
        >>> target = '103P'
        >>> vgas = 0.8 * u.km / u.s
        >>> aper = 30 * u.m
        >>> b = 1.13
        >>> mol_tag = 27001
        >>> transition_freq = (265.886434 * u.GHz).to('MHz')
        >>> integrated_flux = 1.22 * u.K * u.km / u.s

        >>> time = Time('2010-11-3 00:48:06', format='iso')
        >>> ephemobj = Ephem.from_horizons(
        ...     target, epochs=time, closest_apparition=True,
        ...     id_type='designation') # doctest: +REMOTE_DATA

        >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq,
        ...                              mol_tag) # doctest: +REMOTE_DATA

        >>> intl = intensity_conversion(mol_data) # doctest: +REMOTE_DATA
        >>> mol_data.apply([intl.value] * intl.unit,
        ...                name='intl') # doctest: +REMOTE_DATA

        >>> au = einstein_coeff(mol_data) # doctest: +REMOTE_DATA
        >>> mol_data.apply([au.value] * au.unit,
        ...                name='eincoeff') # doctest: +REMOTE_DATA

        >>> lte = LTE()
        >>> q = lte.from_Drahus(integrated_flux, mol_data,
        ...                     ephemobj, vgas, aper, b=b) # doctest: +REMOTE_DATA

        >>> q  # doctest: +REMOTE_DATA +FLOAT_CMP
        <MaskedQuantity 1.09899965e+25 1 / s>


        References
        ----------
        Drahus et al. September 2012. The Sources of HCN and CH3OH and the
        Rotational Temperature in Comet 103P/Hartley 2 from Time-resolved
        Millimeter Spectroscopy. The Astrophysical Journal, Volume 756,
        Issue 1.

        """

        register('Spectroscopy', {'Production Rate (No photodissociation)':
                                  '2012ApJ...756...80D'})

        if not isinstance(mol_data, Phys):
            raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

        t_freq = mol_data['t_freq'][0]
        temp = mol_data['Temperature'][0]
        partition = mol_data['partfn']
        gu = mol_data['dgup'][0]
        eup_J = mol_data['eup_j'][0]
        h = con.h.to('J*s')  # Planck constant
        k = con.k_B.to('J/K')  # Boltzmann constant
        c = con.c.to('m/s')  # speed of light
        vgas = vgas.to('m/s')

        au = mol_data['eincoeff'][0]

        delta = ephemobj["delta"][0]

        calc = ((16*np.pi*k*t_freq.decompose() *
                 partition*vgas) / (np.sqrt(np.pi*np.log(2))
                                    * h * c**2 * au * gu *
                                    np.exp(-eup_J/(k*temp)))).decompose()

        q = integrated_flux*(calc * b * delta / aper)

        q = q.to(u.Hz, equivalencies=u.spectral()).decompose()[0]

        return q


class NonLTE():
    """
    Class for non LTE production rate models.

    """

    @staticmethod
    @requires("pyradex", "astroquery")
    def from_pyradex(integrated_flux, mol_data, line_width=1.0 * u.km / u.s,
                     escapeProbGeom='lvg', iter=100,
                     collider_density={'H2': 900*2.2}):
        """
        Calculate production rate from the Non-LTE iterative code pyradex
        Presently, only the LAMDA catalog is supported by this function.
        In the future a function will be provided by sbpy to build your own
        molecular data file from JPLSpec for use in this function.
        Collider is assumed to be H2O for the default setting since we are only
        considering comet production rates. See documentation for pyradex
        installation tips

        Parameters
        ----------

        integrated_flux : `~astropy.units.Quantity`
            The integrated flux in K * km/s

        mol_data : `sbpy.data.phys`
            `sbpy.data.phys` object that contains AT LEAST the following data:

                    | mol_tag: molecule of interest as string or int JPLSpec identifier
                    | temp: kinetic temperature in gas coma (unit K)
                    | cdensity : cdensity estimate (can be calculated from cdensity_Bockelee) (unit 1/cm^2)
                    | temp_back: (optional) background temperature in K (default is 2.730 K)
                    | lamda_name: (optional) LAMDA molecule identifier to avoid ambiguity. `Lamda.molecule_dict` provides list

            Keywords that can be used for these values are found under
            `~sbpy.data.Conf.fieldnames` documentation. Make sure to inform
            yourself on the values needed for each function, their units, and
            their interchangeable keywords as part of the Phys data class.

        line_width : `~astropy.units.Quantity`
            The FWHM line width (really, the single-zone velocity width to
            scale the column density by: this is most sensibly interpreted as a
            velocity gradient (dv/length)) in km/s (default is 1.0 km/s)

        escapeProbGeom : str
            Which escape probability method to use, available choices are
            'sphere', 'lvg', and 'slab'

        iter : int
            Number of iterations you wish to perform. Default is 100, more
            iterations will take more time to do but will offer a better range
            of results to compare against. The range of guesses is built by
            having the column density guess and subtracting/adding an order of
            magnitude for the start and end values of the loop, respectively.
            i.e.  a guess of 1e15 will create a range between 1e14 and 1e16

        collider_density : dict
            Dictionary of colliders and their densities in cm^-3. Allowed
            entries are any of the following : h2,oh2,ph2,e,He,H,H+
            See `~Pyradex` documentation for more information.
            Default dictionary is {'H2' : 900*2.2} where 900 is the
            collider density of H2 and 2.2 is the value for the
            square root of the ratio of reduced masses of H2O/H2
            as follows:
            (mu_H2O/mu_H2)**0.5 = ((18**2/18*2)/((18*2)/(18+2)))**0.5 = 2.2
            in order to scale the collisional excitation to H2O as the main
            collisional partner. (Walker, et al. 2014; de val Borro, et al.
            2017; & Schoier, et al. 2004)

        Returns
        -------
        column density : `~astropy.units.Quantity`
            column density to use for the Haser model ratio calculation
            Note: it is normal for pyradex/RADEX to output warnings depending
            on the setup the user gives it (such as not having found a
            molecular data file so it searches for it on LAMDA catalog)

        Examples
        --------
        >>> from sbpy.activity import NonLTE
        >>> from sbpy.data import Phys
        >>> import astropy.units as u

        >>> transition_freq = (177.196 * u.GHz).to(u.MHz)
        >>> mol_tag =  29002
        >>> cdensity_guess = (1.89*10.**(14) / (u.cm * u.cm))
        >>> temp_estimate = 20. * u.K

        >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq,
        ...                              mol_tag)  # doctest: +REMOTE_DATA
        >>> mol_data.apply([cdensity_guess.value] * cdensity_guess.unit,
        ...                 name= 'cdensity')  # doctest: +REMOTE_DATA
        >>> mol_data.apply(['HCO+@xpol'],
        ...                 name='lamda_name')  # doctest: +REMOTE_DATA

        >>> nonLTE = NonLTE()
        >>> try:  # doctest: +SKIP
        ...     cdensity = nonLTE.from_pyradex(1.234 * u.K * u.km / u.s,
        ...                                    mol_data, iter=600,
        ...                                    collider_density={'H2': 900})  # doctest: +REMOTE_DATA
        ...     print(cdensity)  # doctest: +REMOTE_DATA
        ... except ImportError:
        ...     pass
        Closest Integrated Flux:[1.24925956] K km / s
        Given Integrated Flux: 1.234 K km / s
        [1.06363773e+14] 1 / cm2


        References
        ----------
        Haser 1957, Bulletin de la Societe Royale des Sciences de Liege
        43, 740.

        Walker, et al., On the Validity of Collider-mass Scaling for Molecular
        Rotational Excitation, APJ, August 2014.

        van der Tak, et al., A computer program for fast non-LTE analysis of
        interstellar line spectra. With diagnostic plots to interpret observed
        line intensity ratios. A&A, February 12 2013.

        de Val Borro, et al., Measuring molecular abundances in comet C/2014 Q2
        (Lovejoy) using the APEX telescope. Monthly Notices of the Royal
        Astronomical Society, October 27 2017.

        Schoier, et al., An atomic and molecular database for analysis of
        submillimetre line observations. A&A, November 4 2004.

        """

        if not isinstance(mol_data, Phys):
            raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

        register('Production Rates', {'Radex': '2007A&A...468..627V'})

        # convert mol_tag JPLSpec identifier to verbose name if needed
        try:
            mol_data['lamda_name']
            name = mol_data['lamda_name'][0]
            name = name.lower()
        except KeyError:
            if not isinstance(mol_data['mol_tag'][0], str):
                cat = JPLSpec.get_species_table()
                mol = cat[cat['TAG'] == mol_data['mol_tag'][0]]
                name = mol['NAME'].data[0]
                name = name.lower()
            else:
                name = mol_data['mol_tag'][0]
                name = name.lower()

        # try various common instances of molecule names and check them against LAMDA before complaining
        try:
            Lamda.molecule_dict[name]
        except KeyError:
            try_name = "{}@xpol".format(name)
            try:
                Lamda.molecule_dict[try_name]
                name = try_name
            except KeyError:
                print('Molecule name {} not found in LAMDA, module tried {} and also\
                       found no molecule with this identifier within LAMDA. Please\
                       enter LAMDA identifiable name using mol_data["lamda_name"]\
                       . Use Lamda.molecule_dict to see all available options.'.format(name, try_name))
                raise

        # define Temperature
        temp = mol_data['temp']

        # check for optional values within mol_data
        if 'temp_back' in mol_data:
            tbackground = mol_data['temp_back']
        else:
            tbackground = 2.730 * u.K

        # define cdensity and iteration parameters
        cdensity = mol_data['cdensity'].to(1 / (u.cm * u.cm))
        cdensity_low = cdensity - (cdensity*0.9)
        cdensity_high = cdensity + (cdensity*9)
        # range for 400 iterations
        cdensity_range = np.linspace(cdensity_low, cdensity_high, iter)
        fluxes = []
        column_density = []

        with tempfile.TemporaryDirectory() as datapath:
            for i in cdensity_range:
                R = pyradex.Radex(column=i, deltav=line_width,
                                  tbackground=tbackground, species=name,
                                  temperature=temp, datapath=datapath,
                                  escapeProbGeom=escapeProbGeom,
                                  collider_densities=collider_density)

                table = R()

                # find closest matching frequency to user defined
                indx = (np.abs(table['frequency']-mol_data['t_freq'])).argmin()
                radexfreq = table['frequency'][indx]
                # get table for that frequency
                values = table[table['frequency'] == radexfreq]
                # use eq in io.f from Pyradex to get integrated flux in K * km/s
                int_flux_pyradex = 1.0645 * values['T_B'] * line_width

                fluxes.append(int_flux_pyradex)
                column_density.append(i)

        # closest matching integrated flux from pyradex

        fluxes = np.array(fluxes)

        index_flux = (
            np.abs(fluxes-integrated_flux.to(u.K * u.km / u.s).value)).argmin()

        # corresponding column density in 1/cm^2
        column_density = column_density[index_flux]
        print('Closest Integrated Flux:{}'.format(
            fluxes[index_flux] * u.K * u.km / u.s))
        print('Given Integrated Flux: {}'.format(integrated_flux))

        return column_density
