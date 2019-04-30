# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Spectroscopy Module

created on June 23, 2017
"""

import numpy as np
from scipy import interpolate
import astropy.constants as con
import astropy.units as u
from astropy.time import Time
from astroquery.jplspec import JPLSpec
from astroquery.jplhorizons import Horizons, conf
from ..bib import register
from ..data import Phys

conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

__all__ = ['Spectrum', 'SpectralModel', 'molecular_data', 'einstein_coeff',
           'intensity_conversion']


def molecular_data(temp_estimate, transition_freq, mol_tag):
    """
    Returns relevant constants from JPLSpec catalog and energy calculations

    Parameters
    ----------
    temp_estimate : `~astropy.units.Quantity`
        Estimated temperature in Kelvins

    transition_freq : `~astropy.units.Quantity`
        Transition frequency in MHz

    mol_tag : int or str
        Molecule identifier. Make sure it is an exclusive identifier.

    Returns
    -------
    Molecular data : `~sbpy.data.Phys` instance
        Quantities in the following order:
            | Transition frequency
            | Temperature
            | Integrated line intensity at 300 K
            | Partition function at 300 K
            | Partition function at designated temperature
            | Upper state degeneracy
            | Upper level energy in Joules
            | Lower level energy in Joules
            | Degrees of freedom

    """

    query = JPLSpec.query_lines(min_frequency=(transition_freq - (1 * u.MHz)),
                                max_frequency=(transition_freq + (1 * u.MHz)),
                                molecule=mol_tag)

    freq_list = query['FREQ']

    t_freq = min(list(freq_list.quantity),
                 key=lambda x: abs(x-transition_freq))

    data = query[query['FREQ'] == t_freq.value]

    df = int(data['DR'].data)

    lgint = float(data['LGINT'].data)

    lgint = 10**lgint * u.nm * u.nm * u.MHz

    elo = float(data['ELO'].data) / u.cm

    gu = float(data['GUP'].data)

    cat = JPLSpec.get_species_table()

    mol = cat[cat['TAG'] == mol_tag]

    temp_list = cat.meta['Temperature (K)'] * u.K

    part = list(mol['QLOG1', 'QLOG2', 'QLOG3', 'QLOG4', 'QLOG5', 'QLOG6',
                    'QLOG7'][0])

    temp = temp_estimate

    f = interpolate.interp1d(temp_list, part, 'linear')

    partition = 10**(f(temp_estimate.value))

    part300 = 10 ** (float(mol['QLOG1'].data))

    # yields in 1/cm
    energy = elo + (t_freq.to(1/u.cm, equivalencies=u.spectral()))

    energy_J = energy.to(u.J, equivalencies=u.spectral())
    elo_J = elo.to(u.J, equivalencies=u.spectral())

    quantities = [t_freq, temp, lgint, part300, partition, gu, energy_J, elo_J, df]
    names = ('Transition frequency',
             'Temperature',
             'Integrated line intensity at 300 K',
             'Partition function at 300 K',
             'Partition function at designated temperature',
             'Upper state degeneracy',
             'Upper level energy in Joules',
             'Lower level energy in Joules',
             'Degrees of freedom')
    result = Phys.from_array(quantities, names)

    return result


def intensity_conversion(temp_estimate, transition_freq, mol_tag):
    """
    Returns conversion of integrated line intensity at designated temperature.

    Parameters
    ----------
    transition_freq : `~astropy.units.Quantity`
        Transition frequency in MHz

    temp_estimate : `~astropy.units.Quantity`
        Estimated temperature in Kelvins

    mol_tag : int or str
        Molecule identifier. Make sure it is an exclusive identifier.

    Returns
    -------
    intl : `~astropy.units.Quantity`
        Integrated line intensity at designated temperature

    """

    mol_data = molecular_data(temp_estimate, transition_freq, mol_tag)

    temp = mol_data[0][1]
    lgint = mol_data[0][2]
    part300 = mol_data[0][3]
    partition = mol_data[0][4]
    energy_J = mol_data[0][6]
    elo_J = mol_data[0][7]
    df = mol_data[0][8]

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


def einstein_coeff(temp_estimate, transition_freq, mol_tag):
    """
    Einstein coefficient from molecular data

    Parameters
    ----------
    transition_freq : `~astropy.units.Quantity`
        Transition frequency in MHz

    temp_estimate : `~astropy.units.Quantity`
        Estimated temperature in Kelvins

    mol_tag : int or str
        Molecule identifier. Make sure it is an exclusive identifier.

    Returns
    -------
    einstein_coeff : `~astropy.units.Quantity`
        Spontaneous emission coefficient

    """

    mol_data = molecular_data(temp_estimate, transition_freq, mol_tag)

    t_freq = mol_data[0][0]
    temp = mol_data[0][1]
    lgint = mol_data[0][2]
    part300 = mol_data[0][3]
    partition = mol_data[0][4]
    gu = mol_data[0][5]
    energy_J = mol_data[0][6]
    elo_J = mol_data[0][7]

    h = con.h.to('J*s')  # Planck constant

    k = con.k_B.to('J/K')  # Boltzmann constant

    intl = intensity_conversion(temp_estimate, transition_freq, mol_tag)

    if (h*t_freq/(k*temp)).decompose().value and \
            (h*t_freq/(k*300*u.K)).decompose().value < 1:

        au = (lgint*t_freq
              * (part300/gu)*np.exp(energy_J / (k*300*u.K))*(1.748e-9)).value

    else:

        au = (intl*(t_freq)**2 *
              (partition/gu)*(np.exp(-(elo_J/(k*temp)).value) -
                              np.exp(-(energy_J/(k*temp)).value))**(-1)
              * (2.7964e-16)).value

    return au / u.s


def photod_rate(ephemobj, mol_tag):
    # imported here to avoid circular dependency with activity.gas
    from ..activity.gas import photo_timescale
    from ..data import Ephem

    if not isinstance(ephemobj, Ephem):
        raise ValueError('ephemobj must be a `sbpy.data.ephem` instance.')

    orb = ephemobj
    delta = (orb['delta']).to('m')
    r = (orb['r'])
    cat = JPLSpec.get_species_table()
    mol = cat[cat['TAG'] == mol_tag]
    name = mol['NAME'].data[0]

    timescale = photo_timescale(name)

    if timescale.ndim != 0:
        # array
        timescale = timescale[0]

    beta = (timescale) * r**2

    result = Phys.from_array((beta, delta), ('beta', 'delta'))

    return result


def total_number(integrated_line, temp_estimate, transition_freq, mol_tag,
                 aper, b, ephemobj):
    """
    Basic equation relating number of molecules with observed integrated flux.
    This is given by equation 10 in
    https://ui.adsabs.harvard.edu/#abs/2004come.book..391B

    Parameters
    ----------
    integrated_line : `~astropy.units.Quantity`
        Integrated flux of emission line.
    transition_freq : `~astropy.units.Quantity`
        Transition frequency

    Returns
    -------
    total_number : float
        Total number of molecules within the aperture

    """
    register('Spectroscopy', {'Total Number (eq. 10)': '2004come.book..391B'})
    cdensity = integrated_line
    cdensity *= (8*np.pi*con.k_B*transition_freq**2 /
                 (con.h*con.c**3 * einstein_coeff(temp_estimate,
                                                  transition_freq,
                                                  mol_tag))).decompose()

    photod = photod_rate(ephemobj, mol_tag)

    beta = photod["beta"]

    sigma = (1./2. * beta * b * con.c / (transition_freq * aper)).value

    total_number = cdensity.decompose() * sigma * u.m * u.m / np.sqrt(np.log(2))

    return total_number


class SpectralModel():
    """Range of spectral models"""
    def haser():
        """Haser model

        should allow direct creation of a `sbpy.actvity.Haser` instance"""
        pass

    def emission_lines():
        """Emission lines"""
        pass

    def reflectance():
        """Reflectance spectrum (asteroids)"""


class Spectrum():

    def __init__(self, flux, dispersionaxis, unit):
        self.flux = flux
        self.axis = dispersionaxis
        self.unit = unit

    @classmethod
    def read(cls, filename, columns='auto'):
        """Read spectrum from file

        Parameters
        ----------
        filename : str, mandatory
            data file name
        columns : str or list-like, optional, default: 'auto'
            file format, `auto` will try to recognize the format
            automatically

        Returns
        -------
        `Spectrum` instance

        Examples
        --------
        >>> spec = Spectrum.read('2012_XY.dat') # doctest: +SKIP

        not yet implemented

        """

    def write(self, filename, columns='all'):
        """Write spectrum to file

        Parameters
        ----------
        filename : str, mandatory
            data file name
        columns : str or list-like, optional: default: 'all'
            file format; `all` will write all fields to the file

        Examples
        --------
        >>> spec = Spectrum.read('2012_XY.dat') # doctest: +SKIP
        >>> spec.write('2012_XY.dat.bak') # doctest: +SKIP

        not yet implemented

        """

    def convert_units(self, **kwargs):
        """Convert Spectrum units as provided by user

        Examples
        --------
        >>> spec.convert_units(flux_unit=u.K) # doctest: +SKIP
        >>> spec.convert_units(dispersion_unit=u.km/u.s) # doctest: +SKIP

        not yet implemented

        """

    def baseline(self, subtract=False):
        """fit baseline to `Spectrum` instance

        Parameters
        ----------
        subtract : bool, optional, default=False
            if `True`, subtract the baseline

        Returns
        -------
        float

        Examples
        --------
        >>> baseline = spec.baseline() # doctest: +SKIP
        >>> spec.baseline(subtract=True) # doctest: +SKIP

        not yet implemented

        """

    def slope(self, subtract=False):
        """fit slope to `Spectrum` instance

        Parameters
        ----------
        subtract : bool, optional, default=False
            if `True`, subtract the slope

        Returns
        -------
        float

        Examples
        --------
        >>> slope = spec.slope() # doctest: +SKIP
        >>> spec.slope(subtract=True) # doctest: +SKIP

        not yet implemented

        """

    def integrated_flux(self, frequency, interval=1*u.km/u.s):
        """
        Calculate integrated flux of emission line.

        Parameters
        ----------
        frequency : `~astropy.units.Quantity`
            Transition frequency
        interval : `~astropy.units.Quantity`
            line width

        Examples
        --------
        >>> flux = spec.integrated_flux(frequency=556.9*u.GHz, # doctest: +SKIP
        >>>                             interval=1.7*u.km/u.s) # doctest: +SKIP

        not yet implemented

        """

    def fit(self, spec):
        """Fit `SpectralModel` to different model types

        Parameters
        ----------
        spec : str, mandatory
            `SpectralModel` instance to fit

        Examples
        --------
        >>> spec_model = SpectralModel(type='Haser', molecule='H2O')        # doctest: +SKIP

        >>> spec.fit(spec_model) # doctest: +SKIP
        >>> print(spec.fit_info) # doctest: +SKIP

        not yet implemented

        """

    def prodrate_np(self, spectra, temp_estimate, transition_freq,
                    mol_tag, ephemobj, vgas=1 * u.km/u.s,
                    aper=25 * u.m, b=1.2):
        """
        | Returns production rate based on Drahus 2012 model referenced. Includes
        | no photodissociation

        Parameters
        ----------
        spectra : `~astropy.units.Quantity`
            Temperature brightness integral derived from spectra

        transition_freq : `~astropy.units.Quantity` or list of this type
            Transition frequency in MHz. If a list of more than one frequency
            is given the function will find the average constants from all the
            lines and use these to calculate the final production rate.
            According to Drahus himself:
            "One advantage of this approach is that the calculated production
            rate is based on the high-S/N observation of the averaged line
            profile and not the much lower S/N of a single line- as
            long as the theory behind holds. In some cases, individual lines
            may not even be detectable, but if their average is, one will be
            able to calculate the production rate for some assumed rotational
            temperature."

        temp_estimate : `~astropy.units.Quantity`
            Estimated temperature in Kelvins

        ephemobj: `~sbpy.data.Ephem` object
            An `~sbpy.data.Ephem` object that holds ephemerides information

        mol_tag : int or str
            Molecule identifier. Make sure it is an exclusive identifier.

        vgas : `~astropy.units.Quantity`
            Gas velocity approximation in km / s. Default is 1 km / s

        aper : `~astropy.units.Quantity`
            Telescope aperture in meters. Default is 25 m

        b : int
            | Dimensionless factor intrinsic to every antenna. Typical
            | value, and the default for this model, is 1.22. See
            | references for more information on this parameter.

        Returns
        -------
        q : `~astropy.units.Quantity`
            Production rate, not including photodissociation

        Examples
        --------
        >>> import astropy.units as u  # doctest: +SKIP

        >>> from astropy.time import Time # doctest: +SKIP

        >>> from sbpy.data import Ephem # doctest: +SKIP

        >>> from sbpy.spectroscopy import prodrate_np  # doctest: +SKIP

        >>> temp_estimate = 33. * u.K  # doctest: +SKIP

        >>> target = '103P'  # doctest: +SKIP

        >>> vgas = 0.8 * u.km / u.s  # doctest: +SKIP

        >>> aper = 30 * u.m  # doctest: +SKIP

        >>> b = 1.13  # doctest: +SKIP

        >>> mol_tag = 27001  # doctest: +SKIP

        >>> transition_freq = 265.886434 * u.MHz  # doctest: +SKIP

        >>> spectra = 1.22 * u.K * u.km / u.s  # doctest: +SKIP

        >>> time = Time('2010-11-3 00:48:06', format='iso')  # doctest: +SKIP

        >>> ephemobj = Ephem(target, epochs=time.jd, id_type='id') # doctest: +SKIP

        >>> q = prodrate_np(spectra, temp_estimate, transition_freq, # doctest: +SKIP
                            mol_tag, ephemobj, vgas, aper, b=b)

        >>> q  # doctest: +SKIP
        <Quantity 1.0432591198553935e+25 1 / s>


        References
        ----------
        Drahus et al. September 2012. The Sources of HCN and CH3OH and the
        Rotational Temperature in Comet 103P/Hartley 2 from Time-resolved
        Millimeter Spectroscopy. The Astrophysical Journal, Volume 756,
        Issue 1.

        """

        register('Spectroscopy', {'Production Rate (No photodissociation)':
                                  '2012ApJ...756...80D'})

        assert isinstance(temp_estimate, u.Quantity)
        assert isinstance(vgas, u.Quantity)
        assert isinstance(aper, u.Quantity)
        assert isinstance(spectra, u.Quantity)  # K * km / s

        if type(transition_freq) == list:

            t_freq_list = []
            temp_list = []
            partition_list = []
            gu_list = []
            energy_J_list = []
            au_list = []
            h = con.h.to('J*s')  # Planck constant
            k = con.k_B.to('J/K')  # Boltzmann constant
            c = con.c.to('m/s')  # speed of light
            vgas = vgas.to('m/s')

            for i in range(0, len(transition_freq)):

                assert isinstance(transition_freq[i], u.Quantity)

                mol_data = molecular_data(temp_estimate, transition_freq[i],
                                          mol_tag)

                t_freq = mol_data[0][0]
                temp = mol_data[0][1]
                partition = mol_data[0][4]
                gu = mol_data[0][5]
                energy_J = mol_data[0][6]

                au = einstein_coeff(temp_estimate, transition_freq[i], mol_tag)

                au_list.append(au)
                t_freq_list.append(t_freq)
                temp_list.append(temp)
                partition_list.append(partition)
                gu_list.append(gu)
                energy_J_list.append(energy_J)

            t_freq = t_freq_list[0]
            temp = temp_list[0]
            au = sum(au_list) / float(len(au_list))
            partition = sum(partition_list) / float(len(partition_list))
            gu = sum(gu_list) / float(len(gu_list))
            energy_J = sum(energy_J_list) / float(len(energy_J_list))

            photod = photod_rate(ephemobj, mol_tag)

            delta = photod["delta"]

            calc = ((16*np.pi*k*t_freq.decompose() *
                     partition*vgas) / (np.sqrt(np.pi*np.log(2))
                                        * h * c**2 * au * gu *
                                        np.exp(-energy_J/(k*temp)))).decompose()

            q = spectra*(calc * b * delta / aper)

            q = q.to(u.Hz, equivalencies=u.spectral()).decompose()[0]

        else:

            assert isinstance(transition_freq, u.Quantity)

            mol_data = molecular_data(temp_estimate, transition_freq, mol_tag)

            t_freq = mol_data[0][0]
            temp = mol_data[0][1]
            partition = mol_data[0][4]
            gu = mol_data[0][5]
            energy_J = mol_data[0][6]
            h = con.h.to('J*s')  # Planck constant
            k = con.k_B.to('J/K')  # Boltzmann constant
            c = con.c.to('m/s')  # speed of light
            vgas = vgas.to('m/s')

            au = einstein_coeff(temp_estimate, transition_freq, mol_tag)

            photod = photod_rate(ephemobj, mol_tag)

            delta = photod["delta"]

            calc = ((16*np.pi*k*t_freq.decompose() *
                     partition*vgas) / (np.sqrt(np.pi*np.log(2))
                                        * h * c**2 * au * gu *
                                        np.exp(-energy_J/(k*temp)))).decompose()

            q = spectra*(calc * b * delta / aper)

            q = q.to(u.Hz, equivalencies=u.spectral()).decompose()[0]

        return q

    def production_rate(self, coma, integrated_line, temp_estimate,
                        transition_freq, mol_tag, ephemobj,
                        aper=25 * u.m, b=1.2):
        """
        Calculate production rate for `GasComa`

        Parameters
        ----------
        coma : `sbpy.activity.gas.GasComa`
            Gas coma model

        integrated_line : `~astropy.units.Quantity`
            Temperature brightness integral derived from spectra

        temp_estimate : `~astropy.units.Quantity`
            Estimated temperature in Kelvins

        transition_freq : `~astropy.units.Quantity`
            Transition frequency being analyzed

        mol_tag : int or str
            Molecule identifier. Make sure it is an exclusive identifier.

        ephemobj: `~sbpy.data.Ephem` object
            An `~sbpy.data.Ephem` object that holds ephemerides information

        aper : `~astropy.units.Quantity`
            Telescope aperture in meters. Default is 25 m

        b : int
            | Dimensionless factor intrinsic to every antenna. Typical
            | value, and the default for this model, is 1.22. See
            | references for more information on this parameter.

        Returns
        -------
        Q : `~astropy.units.Quantity`
            production rate

        Examples
        --------
        >>> from sbpy.activity.gas import Haser # doctest: +SKIP
        >>> from sbpy.data import Ephem # doctest: +SKIP
        >>> from astropy.time import Time # doctest: +SKIP

        >>> coma = Haser(Q, v, parent) # doctest: +SKIP
        >>> Q = spec.production_rate(coma, molecule='H2O') # doctest: +SKIP

        >>> Q_estimate = 2.8*10**(28) / u.s # doctest: +SKIP
        >>> transition_freq = (230.53799 * u.GHz).to('MHz') # doctest: +SKIP
        >>> aper = 10 * u.m # doctest: +SKIP
        >>> mol_tag = 28001 # doctest: +SKIP
        >>> temp_estimate = 25. * u.K # doctest: +SKIP
        >>> target = 'C/2016 R2' # doctest: +SKIP
        >>> b = 0.74 # doctest: +SKIP
        >>> vgas = 0.5 * u.km / u.s # doctest: +SKIP

        >>> time = Time('2017-12-22 05:24:20', format = 'iso') # doctest: +SKIP
        >>> ephemobj = Ephem(target, epochs=time.jd) # doctest: +SKIP
        >>> spectra = 0.26 * u.K * u.km / u.s # doctest: +SKIP

        >>> parent = photo_timescale('CO') * vgas # doctest: +SKIP

        >>> coma = Haser(Q_estimate, vgas, parent) # doctest: +SKIP

        >>> Q = spec.production_rate(coma, spectra, temp_estimate, # doctest: +SKIP
                                     transition_freq, mol_tag, ephemobj,
                                     aper=aper, b=b) # doctest: +SKIP

        >>> print(Q) # doctest: +SKIP
            <Quantity [1.64403219e+28] 1 / s>

        References
        ----------
        Haser 1957, Bulletin de la Societe Royale des Sciences de Liege
        43, 740.
        Newburn and Johnson 1978, Icarus 35, 360-368.

        """

        from ..activity.gas import GasComa

        if not isinstance(coma, GasComa):
            raise ValueError('coma must be a GasComa instance.')

        # integrated_line = self.integrated_flux(transition_freq) - not yet implemented

        molecules = total_number(integrated_line, temp_estimate,
                                 transition_freq, mol_tag, aper, b, ephemobj)

        model_molecules = coma.total_number(aper)

        Q = coma.Q * molecules/model_molecules

        return Q

    def plot(self):
        """Plot spectrum

        Returns
        -------
        `matplotlib.pyplot` instance

        not yet implemented
        """
