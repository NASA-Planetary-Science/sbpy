# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Production Rate Module

created on June 26, 2019
"""

import numpy as np
import astropy.constants as con
import astropy.units as u
from astropy.time import Time
from astroquery.jplhorizons import Horizons, conf
from astroquery.jplspec import JPLSpec
from ...bib import register
from ...data import Phys

conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

__all__ = ['LTE', 'NonLTE', 'einstein_coeff',
           'intensity_conversion', 'beta_factor', 'total_number_nocd']


def intensity_conversion(mol_data):
    """
    Returns conversion of the integrated line intensity at 300 K
    (from the JPL molecular spectra catalog) to a chosen temperature

    Parameters
    ----------
    mol_data : `~sbpy.data.phys`
        `~sbpy.data.phys` object that contains the following data,
        using `~astropy.units` for the required units:

            .. _making-a-list:

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
        `~sbpy.data.conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.phys.from_jplspec` to obtain these values in order to maintain
        consistency and because all calculations can be handled within `sbpy`
        scope if JPLSpec is used. Yet, if you wish to use your own molecular data,
        it is possible. Make sure to inform yourself on the values needed for each
        function, their units, and their interchangeable keywords as part of
        the `~sbpy.data.phys` data class.

    Returns
    -------
    intl : `~astropy.Quantity`
        Integrated line intensity at designated temperature in MHz * nm**2,
        which can be appended to the original `sbpy.phys` object for future
        calculations

    """

    if not isinstance(mol_data, Phys):
        raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

    temp = mol_data['temp'][0]
    lgint = mol_data['lgint300'][0]
    part300 = mol_data['partfn300'][0]
    partition = mol_data['partfn'][0]
    energy_J = mol_data['eup_j'][0]
    elo_J = mol_data['elo_J'][0]
    df = mol_data['degfr'][0]

    k = con.k_B.to('J/K')  # Boltzmann constant

    if (((energy_J - elo_J).value / (k * temp).value) and
            ((energy_J - elo_J).value / (k * 300 * u.K).value) < 1):

        if df == 0 or 2:

            intl = lgint*(300*u.K / temp)**(2)*np.exp(-(1/temp - 1/(300*u.K))
                                                      * elo_J/k)

        else:

            intl = lgint*(300*u.K/temp)**(5/2)*np.exp(-(1/temp - 1/(300*u.K))
                                                      * elo_J/k)

    else:

        intl = lgint*(part300/partition)*(np.exp(-elo_J/(k*temp)) -
                                          np.exp(-energy_J/(k*temp))) / \
            (np.exp(-elo_J/(k*300 * u.K)) - np.exp(-energy_J/(k*300*u.K)))

    return intl


def einstein_coeff(mol_data):
    """
    Einstein coefficient from molecular data

    Parameters
    ----------

    mol_data : `~sbpy.data.phys`
        `~sbpy.data.phys` object that contains the following data,
        using `astropy.units` for the required units:

            .. _making-a-list:

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
        `~sbpy.data.conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.phys.from_jplspec` to obtain these values in order to maintain
        consistency. Yet, if you wish to use your own molecular data, it is
        possible. Make sure to inform yourself on the values needed for each
        function, their units, and their interchangeable keywords as part of
        the data class.

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
    energy_J = mol_data['eup_j'][0]
    elo_J = mol_data['elo_J'][0]
    df = mol_data['degfreedom'][0]
    t_freq = mol_data['t_freq'][0]
    gu = mol_data['dgup'][0]

    h = con.h.to('J*s')  # Planck constant

    k = con.k_B.to('J/K')  # Boltzmann constant

    intl = mol_data['lgint'][0]

    if (h*t_freq/(k*temp)).decompose().value and \
            (h*t_freq/(k*300*u.K)).decompose().value < 1:

        au = (lgint*t_freq
              * (part300/gu)*np.exp(energy_J / (k*300*u.K))*(1.748e-9)).value

    else:

        au = (intl*(t_freq)**2 *
              (partition/gu)*(np.exp(-(elo_J/(k*temp)).value) -
                              np.exp(-(energy_J/(k*temp)).value))**(-1)
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
        The molecule MUST be defined in `sbpy.activity.gas.timescale`, otherwise
        this function cannot be used and the beta factor
        will have to be provided by the user directly for calculations. The
        user can obtain the beta factor from the formula provided above.
        Keywords that can be used for these values are found under
        `~sbpy.data.conf.fieldnames` documentation. We recommend the use of the
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


def total_number_nocd(integrated_flux, mol_data, aper, b):
    """
    Basic equation relating number of molecules with observed integrated flux
    without the need for column density to be given
    This is given by equation 10 in
    https://ui.adsabs.harvard.edu/#abs/2004come.book..391B
    and is derived from data from JPLSpec, feel free to use your own total number
    to calculate production rate or use this function with your own molecular data
    as long as you are aware of the needed data

    Parameters
    ----------
    integrated_line : `~astropy.units.Quantity`
        Integrated flux of emission line.

    mol_data : `sbpy.data.phys`
        `sbpy.data.phys` object that contains AT LEAST the following data:
                | Transition frequency in MHz
                | Einstein Coefficient (1/s)
                | Beta: (Timescale (in s) * r^2 (in au))
                | Optional: user-defined Column Density in 1/m^2

        The values above can either be given by the user or obtained from the
        functions `~sbpy.activity.gas.productionrate.einstein_coeff` and
        `~sbpy.activity.gas.productionrate.beta_factor`
        Keywords that can be used for these values are found under
        `~sbpy.data.conf.fieldnames` documentation. We recommend the use of the
        JPL Molecular Spectral Catalog and the use of
        `~sbpy.data.phys.from_jplspec` to obtain
        these values in order to maintain consistency. Yet, if you wish to
        use your own molecular data, it is possible. Make sure to inform
        yourself on the values needed for each function, their units, and
        their interchangeable keywords as part of the Phys data class.

    Returns
    -------
    total_number : float
        Total number of molecules within the aperture (Dimensionless)

    """

    if not isinstance(mol_data, Phys):
        raise ValueError('mol_data must be a `sbpy.data.phys` instance.')

    register('Spectroscopy', {'Total Number (eq. 10)': '2004come.book..391B'})

    if mol_data['Column Density'][0].value:
        cdensity = mol_data['Column Density'][0]
    else:
        cdensity = integrated_flux
        cdensity *= (8*np.pi*con.k_B*mol_data['t_freq'][0]**2 /
                     (con.h*con.c**3 * mol_data['eincoeff'][0])).decompose()

    beta = mol_data['beta'][0]

    sigma = (1./2. * beta * b * con.c / (mol_data['t_freq'][0] * aper)).value

    total_number = cdensity.decompose() * sigma * u.m * u.m / np.sqrt(np.log(2))

    return total_number


class LTE():
    """ LTE Methods for calculating production_rate """

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
            `~sbpy.data.conf.fieldnames` documentation. We recommend the use of the
            JPL Molecular Spectral Catalog and the use of
            `~sbpy.data.Phys.from_jplspec` to obtain
            these values in order to maintain consistency. Yet, if you wish to
            use your own molecular data, it is possible. Make sure to inform
            yourself on the values needed for each function, their units, and
            their interchangeable keywords as part of the Phys data class.

        ephemobj : `sbpy.data.ephem`
            `sbpy.data.ephem` object holding ephemeride information including
            distance from comet to Sun ['r'] and from comet to observer ['delta']

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
        >>> import astropy.units as u  # doctest: +SKIP
        >>> from astropy.time import Time # doctest: +SKIP
        >>> from sbpy.data import Ephem, Phys # doctest: +SKIP
        >>> from sbpy.activity import (LTE, einstein_coeff,
        ...                            intensity_conversion, beta_factor)  # doctest: +SKIP

        >>> temp_estimate = 47. * u.K  # doctest: +SKIP
        >>> target = '103P'  # doctest: +SKIP
        >>> vgas = 0.8 * u.km / u.s  # doctest: +SKIP
        >>> aper = 30 * u.m  # doctest: +SKIP
        >>> b = 1.13  # doctest: +SKIP
        >>> mol_tag = 27001  # doctest: +SKIP
        >>> transition_freq = (265.886434 * u.GHz).to('MHz')  # doctest: +SKIP
        >>> integrated_flux = 1.22 * u.K * u.km / u.s  # doctest: +SKIP

        >>> time = Time('2010-11-3 00:48:06', format='iso')  # doctest: +SKIP
        >>> ephemobj = Ephem.from_horizons(target, epochs=time.jd, id_type='id') # doctest: +SKIP

        >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag) # doctest: +SKIP

        >>> intl = intensity_conversion(mol_data) # doctest: +SKIP
        >>> mol_data.table.add_column(Column([intl.value] * intl.unit,
        ...                                  name='intl')) # doctest: +SKIP

        >>> au = einstein_coeff(mol_data) # doctest: +SKIP
        >>> mol_data.table.add_column(Column([au.value] * au.unit,
        ...                                  name='eincoeff')) # doctest: +SKIP

        >>> lte = LTE() # doctest: +SKIP
        >>> q = lte.from_Drahus(integrated_flux, mol_data,
        ...                     ephemobj, vgas, aper, b=b)  # doctest: +SKIP

        >>> q  # doctest: +SKIP
        <Quantity 1.05828096e+25 1 / s>


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
        energy_J = mol_data['eup_j'][0]
        h = con.h.to('J*s')  # Planck constant
        k = con.k_B.to('J/K')  # Boltzmann constant
        c = con.c.to('m/s')  # speed of light
        vgas = vgas.to('m/s')

        au = mol_data['eincoeff'][0]

        delta = ephemobj["delta"][0]

        calc = ((16*np.pi*k*t_freq.decompose() *
                 partition*vgas) / (np.sqrt(np.pi*np.log(2))
                                    * h * c**2 * au * gu *
                                    np.exp(-energy_J/(k*temp)))).decompose()

        q = integrated_flux*(calc * b * delta / aper)

        q = q.to(u.Hz, equivalencies=u.spectral()).decompose()[0]

        return q

    def from_Haser(self, coma, mol_data, aper=25 * u.m):
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
                    | `~sbpy.activity.gas.total_number_nocd` for a calculation
                    | of this datum if you don't wish to provide it yourself)

            This field can be given by the user directly or calculated using the
            necessary combinations of the following functions:
            `~sbpy.data.phys.from_jplspec`,
            `~sbpy.activity.gas.productionrate.einstein_coeff`,
            `~sbpy.activity.gas.productionrate.beta_factor`, and
            `~sbpy.activity.gas.productionrate.total_number`.
            Keywords that can be used for these values are found under
            `~sbpy.data.conf.fieldnames` documentation. We recommend the use of the
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
        >>> import astropy.units as u  # doctest: +SKIP
        >>> from astropy.time import Time # doctest: +SKIP
        >>> from sbpy.data import Ephem, Phys # doctest: +SKIP
        >>> from sbpy.activity import Haser, LTE, photo_timescale, einstein_coeff # doctest: +SKIP
        >>> from sbpy.activity import intensity_conversion, beta_factor, total_number_nocd  # doctest: +SKIP

        >>> aper = 10 * u.m # doctest: +SKIP
        >>> mol_tag = 28001 # doctest: +SKIP
        >>> temp_estimate = 25. * u.K # doctest: +SKIP
        >>> target = 'C/2016 R2' # doctest: +SKIP
        >>> b = 0.74 # doctest: +SKIP
        >>> vgas = 0.5 * u.km / u.s # doctest: +SKIP
        >>> transition_freq = (230.53799 * u.GHz).to('MHz') # doctest: +SKIP
        >>> integrated_flux = 0.26 * u.K * u.km / u.s # doctest: +SKIP

        >>> time = Time('2017-12-22 05:24:20', format = 'iso') # doctest: +SKIP
        >>> ephemobj = Ephem.from_horizons(target, epochs=time.jd) # doctest: +SKIP

        >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag) # doctest: +SKIP

        >>> intl = intensity_conversion(mol_data) # doctest: +SKIP
        >>> mol_data.table.add_column(Column([intl.value] * intl.unit,
        ...                                  name='intl')) # doctest: +SKIP

        >>> au = einstein_coeff(mol_data) # doctest: +SKIP
        >>> mol_data.table.add_column(Column([au.value] * au.unit,
        ...                                  name='eincoeff')) # doctest: +SKIP

        >>> beta = beta_factor(mol_data, ephemobj) # doctest: +SKIP
        >>> mol_data.table.add_column(Column([beta.value] * beta.unit,
        ...                                  name='beta')) # doctest: +SKIP

        >>> tnum = total_number_nocd(integrated_flux, mol_data, aper, b) # doctest: +SKIP
        >>> mol_data.table.add_column(Column([tnum],
        ...                                  name='total_number_nocd')) # doctest: +SKIP

        >>> Q_estimate = 2.8*10**(28) / u.s # doctest: +SKIP
        >>> parent = photo_timescale('CO') * vgas # doctest: +SKIP
        >>> coma = Haser(Q_estimate, vgas, parent) # doctest: +SKIP

        >>> lte = LTE() # doctest: +SKIP
        >>> Q = lte.from_Haser(coma, mol_data, aper=aper) # doctest: +SKIP

        >>> Q # doctest: +SKIP
            <Quantity [[9.35795579e+27]] 1 / s>

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

        # integrated_line = self.integrated_flux(transition_freq) - not yet implemented

        molecules = mol_data['total_number_nocd']

        model_molecules = coma.total_number(aper)

        Q = coma.Q * molecules/model_molecules

        return Q


class NonLTE():
    """
    Class method for non LTE production rate models
    Not Yet implemented

    """
