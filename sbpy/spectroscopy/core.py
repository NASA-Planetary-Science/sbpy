# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Spectroscopy Module

created on June 23, 2017
"""

from abc import ABC
import numpy as np
from scipy import interpolate
import astropy.constants as con
import astropy.units as u
from astropy.time import Time
from astroquery.jplspec import JPLSpec
from astroquery.jplhorizons import Horizons, conf
from ..bib import register
from ..exceptions import SinglePointSpectrumError

conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

__all__ = ['Spectrum', 'SpectralModel', 'molecular_data', 'einstein_coeff',
           'intensity_conversion']


def molecular_data(temp_estimate, transition_freq, mol_tag):
    """
    Returns relevant constants from JPLSpec catalog and energy calculations

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
    Molecular data : list
        List of constants in the following order:
            | Transtion frequency
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

    result = []

    result.append(t_freq)

    result.extend((temp, lgint, part300, partition, gu, energy_J,
                   elo_J, df))

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

    temp = mol_data[1]
    lgint = mol_data[2]
    part300 = mol_data[3]
    partition = mol_data[4]
    energy_J = mol_data[6]
    elo_J = mol_data[7]
    df = mol_data[8]

    k = con.k_B.to('J/K')  # boltzman constant

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

    t_freq = mol_data[0]
    temp = mol_data[1]
    lgint = mol_data[2]
    part300 = mol_data[3]
    partition = mol_data[4]
    gu = mol_data[5]
    energy_J = mol_data[6]
    elo_J = mol_data[7]

    h = con.h.to('J*s')  # planck constant

    k = con.k_B.to('J/K')  # boltzman constant

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


def photod_rate(time, time_scale, target, id_type, observatory, time_format,
                mol_tag):
    # imported here to avoid circular dependency with activity.gas
    from ..activity.gas import photo_timescale

    epoch = Time(time, scale=time_scale, format=time_format)
    obj = Horizons(id=target, epochs=epoch.jd, location=observatory,
                   id_type=id_type)

    try:
        orb = obj.ephemerides()

    except ValueError:

        raise

    delta = (orb['delta'].data * u.au).to('m')
    r = (orb['r'].data)
    cat = JPLSpec.get_species_table()
    mol = cat[cat['TAG'] == mol_tag]
    name = mol['NAME'].data[0]

    timescale = photo_timescale(name)

    if timescale.ndim != 0:
        # array
        timescale = timescale[0]

    beta = (timescale) * r**2

    return beta, delta


def total_number(integrated_line, temp_estimate, transition_freq, mol_tag,
                 aper, b, time, target, time_scale, id_type, observatory,
                 time_format):
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

    delta, beta = photod_rate(time, time_scale, target, id_type, observatory,
                              time_format, mol_tag)

    sigma = (1./2. * delta * b * con.c / (transition_freq * aper)).value

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
                    mol_tag, time, target, vgas=1 * u.km/u.s,
                    aper=25 * u.m, observatory='500', b=1.2,
                    time_format='iso', time_scale='utc',
                    id_type='designation'):
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

        time : str
            Time of observation of any format supported by `~astropy.time`

        target : str
            | Target designation, if there is more than one aparition you
            | will be prompted to pick a more specific identifier from a
            | displayed table and change the parameter id_type to 'id'.
            | Look at `~astroquery.jplhorizons` for more information.

        mol_tag : int or str
            Molecule identifier. Make sure it is an exclusive identifier.

        vgas : `~astropy.units.Quantity`
            Gas velocity approximation in km / s. Default is 1 km / s

        aper : `~astropy.units.Quantity`
            Telescope aperture in meters. Default is 25 m

        observatory : str
            | Observatory identifier as per `~astroquery.jplhorizons`
            | Default is geocentric ('500')

        b : int
            | Dimensionless factor intrinsic to every antenna. Typical
            | value, and the default for this model, is 1.22. See
            | references for more information on this parameter.

        time_format : str
            | Time format, see `~astropy.time` for more information
            | Default is 'iso' which corresponds to 'YYYY-MM-DD HH:MM:SS'

        time_scale : str
            | Time scale, see `~astropy.time` for mor information.
            | Default is 'utc'

        id_type : str
            | ID type for target. See `~astroquery.jplhorizons` for more.
            | Default is 'designation'

        Returns
        -------
        q : list
            Production rate, not including photodissociation

        Examples
        --------
        >>> import astropy.units as u  # doctest: +SKIP

        >>> from sbpy.spectroscopy import prodrate_np  # doctest: +SKIP

        >>> temp_estimate = 33. * u.K  # doctest: +SKIP

        >>> target = '900918'  # doctest: +SKIP

        >>> vgas = 0.8 * u.km / u.s  # doctest: +SKIP

        >>> aper = 30 * u.m  # doctest: +SKIP

        >>> b = 1.13  # doctest: +SKIP

        >>> mol_tag = 27001  # doctest: +SKIP

        >>> transition_freq = 265.886434 * u.MHz  # doctest: +SKIP

        >>> spectra = 1.22 * u.K * u.km / u.s  # doctest: +SKIP

        >>> time = '2010-11-3 00:48:06'  # doctest: +SKIP

        >>> q = prodrate_np(spectra, temp_estimate, transition_freq, # doctest: +SKIP
                            mol_tag, time, target, vgas, aper,
                            b=b, id_type='id')

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
        assert isinstance(time, str), "Input time as string, i.e. '2018-05-14'"
        assert isinstance(time_scale, str), "Input time scale as string, i.e.\
                                             'utc' see astropy.time.Time for \
                                              more info"
        assert isinstance(target, str), "Object name should be a string and\
                                         should be identifiable through \
                                         astroquery.jplhorizons"
        assert isinstance(observatory, str), "Observatory should be a string\
                                              and identified through \
                                              astroquery.jplhorizons,\
                                              i.e. '500' (geocentric)"

        if type(transition_freq) == list:

            t_freq_list = []
            temp_list = []
            partition_list = []
            gu_list = []
            energy_J_list = []
            au_list = []
            h = con.h.to('J*s')  # planck constant
            k = con.k_B.to('J/K')  # boltzman constant
            c = con.c.to('m/s')  # speed of light
            vgas = vgas.to('m/s')

            for i in range(0, len(transition_freq)):

                assert isinstance(transition_freq[i], u.Quantity)

                mol_data = molecular_data(temp_estimate, transition_freq[i],
                                          mol_tag)

                t_freq = mol_data[0]
                temp = mol_data[1]
                partition = mol_data[4]
                gu = mol_data[5]
                energy_J = mol_data[6]

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

            beta, delta = photod_rate(time, time_scale, target, id_type,
                                      observatory, time_format, mol_tag)

            calc = ((16*np.pi*k*t_freq.decompose() *
                     partition*vgas) / (np.sqrt(np.pi*np.log(2))
                                        * h * c**2 * au * gu *
                                        np.exp(-energy_J/(k*temp)))).decompose()

            q = spectra*(calc * b * delta / aper)

            q = q.to(u.Hz, equivalencies=u.spectral()).decompose()[0]

        else:

            assert isinstance(transition_freq, u.Quantity)

            mol_data = molecular_data(temp_estimate, transition_freq, mol_tag)

            t_freq = mol_data[0]
            temp = mol_data[1]
            partition = mol_data[4]
            gu = mol_data[5]
            energy_J = mol_data[6]
            h = con.h.to('J*s')  # planck constant
            k = con.k_B.to('J/K')  # boltzman constant
            c = con.c.to('m/s')  # speed of light
            vgas = vgas.to('m/s')

            au = einstein_coeff(temp_estimate, transition_freq, mol_tag)

            beta, delta = photod_rate(time, time_scale, target, id_type,
                                      observatory, time_format, mol_tag)

            calc = ((16*np.pi*k*t_freq.decompose() *
                     partition*vgas) / (np.sqrt(np.pi*np.log(2))
                                        * h * c**2 * au * gu *
                                        np.exp(-energy_J/(k*temp)))).decompose()

            q = spectra*(calc * b * delta / aper)

            q = q.to(u.Hz, equivalencies=u.spectral()).decompose()[0]

        return q

    def production_rate(self, coma, integrated_line, temp_estimate,
                        transition_freq, mol_tag, time, target,
                        aper=25 * u.m, observatory='500', b=1.2,
                        time_format='iso', time_scale='utc',
                        id_type='designation'):
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

        time : str
            Time of observation of any format supported by `~astropy.time`

        target : str
            | Target designation, if there is more than one aparition you
            | will be prompted to pick a more specific identifier from a
            | displayed table and change the parameter id_type to 'id'.
            | Look at `~astroquery.jplhorizons` for more information.

        aper : `~astropy.units.Quantity`
            Telescope aperture in meters. Default is 25 m

        observatory : str
            | Observatory identifier as per `~astroquery.jplhorizons`
            | Default is geocentric ('500')

        b : int
            | Dimensionless factor intrinsic to every antenna. Typical
            | value, and the default for this model, is 1.22. See
            | references for more information on this parameter.

        time_format : str
            | Time format, see `~astropy.time` for more information
            | Default is 'iso' which corresponds to 'YYYY-MM-DD HH:MM:SS'

        time_scale : str
            | Time scale, see `~astropy.time` for mor information.
            | Default is 'utc'

        id_type : str
            | ID type for target. See `~astroquery.jplhorizons` for more.
            | Default is 'designation'

        Returns
        -------
        Q : `~astropy.units.Quantity`
            production rate

        Examples
        --------
        >>> from sbpy.activity.gas import Haser # doctest: +SKIP
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

        >>> time = '2017-12-22 05:24:20' # doctest: +SKIP
        >>> spectra = 0.26 * u.K * u.km / u.s # doctest: +SKIP

        >>> parent = photo_timescale('CO') * vgas # doctest: +SKIP

        >>> coma = Haser(Q_estimate, vgas, parent) # doctest: +SKIP

        >>> Q = spec.production_rate(coma, spectra, temp_estimate, # doctest: +SKIP
                                     transition_freq, mol_tag, time, target,
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
                                 transition_freq, mol_tag, aper, b, time, target,
                                 time_scale, id_type, observatory, time_format)

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


class SpectralSource(ABC):
    """Abstract base class for SBPy spectral sources.

    Requires `~synphot`.


    Parameters
    ----------
    source : `~synphot.SourceSpectrum`
        The source spectrum.

    description : string, optional
        A brief description of the source spectrum.


    Attributes
    ----------
    wave        - Wavelengths of the source spectrum.
    fluxd       - Source spectrum.
    description - Brief description of the source spectrum.
    meta        - Meta data from ``source``, if any.

    """

    def __init__(self, source, description=None):
        try:
            import synphot
        except ImportError:
            raise ImportError(
                'synphot required for {}.'.format(self.__class__.__name__))

        self._source = source
        self._description = description

    @classmethod
    def from_array(cls, wave, fluxd, meta=None, **kwargs):
        """Create standard from arrays.


        Parameters
        ----------
        wave : `~astropy.units.Quantity`
            The spectral wavelengths.

        fluxd : `~astropy.units.Quantity`
            The spectral flux densities.

        meta : dict, optional
            Meta data.

        **kwargs
            Passed to object initialization.

        """

        try:
            import synphot
        except ImportError:
            raise ImportError(
                'synphot required for {}.'.format(cls.__name__))

        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=wave, lookup_table=fluxd,
            meta=meta)

        return cls(source, **kwargs)

    @classmethod
    def from_file(cls, filename, wave_unit=None, flux_unit=None,
                  cache=True, **kwargs):
        """Load the source spectrum from a file.

        NaNs are dropped.


        Parameters
        ----------
        filename : string
            The name of the file.  See
            `~synphot.SourceSpectrum.from_file` for details.

        wave_unit, flux_unit : str or `~astropy.units.Unit`, optional
            Wavelength and flux units in the file.

        cache : bool, optional
            If ``True``, cache the contents of URLs.

        **kwargs
            Passed to object initialization.

        """

        try:
            import synphot
        except ImportError:
            raise ImportError(
                'synphot required for {}.'.format(cls.__name__))

        from astropy.utils.data import download_file, _is_url

        if filename.lower().endswith(('.fits', '.fit', '.fz')):
            read_spec = synphot.specio.read_fits_spec
        else:
            read_spec = synphot.specio.read_ascii_spec

        # URL cache because synphot.SourceSpectrum.from_file does not
        if _is_url(filename):
            fn = download_file(filename, cache=True)
        else:
            fn = filename

        spec = read_spec(fn, wave_unit=wave_unit, flux_unit=flux_unit)
        i = np.isfinite(spec[1] * spec[2])
        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=spec[1][i], lookup_table=spec[2][i],
            meta={'header': spec[0]})

        return cls(source, **kwargs)

    @property
    def description(self):
        """Description of the source spectrum."""
        return self._description

    @property
    def wave(self):
        """Wavelengths of the source spectrum."""
        return self.source.waveset

    @property
    def fluxd(self):
        """The source spectrum."""
        return self.source(self._source.waveset, flux_unit='W / (m2 um)')

    @property
    def source(self):
        return self._source

    @property
    def meta(self):
        self._source.meta

    def __call__(self, wave_or_freq, unit=None):
        """Evaluate or interpolate the source spectrum.


        Parameters
        ----------
        wave_or_freq : `~astropy.units.Quantity`
            Requested wavelengths or frequencies of the resulting
            spectrum.

        unit : string, `~astropy.units.Unit`, optional
            Spectral units of the output (flux density).  If ``None``,
            the default depends on ``wave_or_freq``: W/(m2 μm) for
            wavelengths, Jy for frequencies.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`
            The spectrum evaluated at or interpolated to the requested
            wavelengths or frequencies.

        """

        from .. import units as sbu  # avoid circular dependency

        if unit is None:
            if wave_or_freq.unit.is_equivalent('m'):
                unit = u.Unit('W/(m2 um)')
            else:
                unit = u.Jy
        else:
            unit = u.Unit(unit)

        if unit.is_equivalent(sbu.VEGA):
            fluxd = self.source(wave_or_freq, 'W/(m2 um)').to(
                unit, sbu.spectral_density_vega(wave_or_freq))
        else:
            fluxd = self.source(wave_or_freq, unit)

        return fluxd

    def observe(self, wfb, unit=None, **kwargs):
        """Observe source as through filters or spectrometer.


        Parameters
        ----------
        wfb : `~astropy.units.Quantity`, `~synphot.SpectralElement`
            Wavelengths, frequencies, or bandpasses.  May also be a
            list of ``SpectralElement``s.

        unit : string, `~astropy.units.Unit`, optional
            Spectral units of the output (flux density).  If ``None``,
            the default depends on ``wfb``: W/(m2 μm) for wavelengths
            or bandpasses, Jy for frequencies.

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        fluxd : `~astropy.units.Quantity`
            The spectrum rebinned.


        Raises
        ------
        SinglePointSpectrumError - If requested wavelengths or
            frequencies has only one value.


        Notes
        -----
        Method adapted from AstroBetter post by Jessica Lu:
        http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/

        """

        import synphot
        from .. import units as sbu  # avoid circular dependency

        # promote single bandpasses to a list, but preserve number of
        # dimensions
        if isinstance(wfb, (synphot.SpectralElement, str)):
            ndim = 0
            wfb = [wfb]
        else:
            ndim = np.ndim(wfb)

        if isinstance(wfb, (tuple, list)):
            if unit is None:
                unit = u.Unit('W/(m2 um)')
            else:
                unit = u.Unit(unit)

            fluxd = np.ones(len(wfb)) * unit
            for i in range(len(wfb)):
                fluxd[i] = self.filt(wfb[i], unit=unit, **kwargs)[1]
        else:
            if np.size(wfb) == 1:
                raise SinglePointSpectrumError(
                    'Multiple wavelengths or frequencies required for '
                    'observe.  Consider interpolation with {}() instead.'
                    .format(self.__class__.__name__))

            if unit is None:
                if wfb.unit.is_equivalent('m'):
                    unit = u.Unit('W/(m2 um)')
                else:
                    unit = u.Jy
            else:
                unit = u.Unit(unit)

            specele = synphot.SpectralElement(synphot.ConstFlux1D(1))

            # Use force='taper' to prevent PartialOverlap execption.
            # Specele is defined over all wavelengths, but most
            # spectral standards are not.
            kwargs['force'] = kwargs.get('force', 'taper')

            obs = synphot.Observation(
                self.source, specele, binset=wfb, **kwargs)

            if unit.is_equivalent(sbu.VEGAmag):
                fluxd = obs.sample_binned(flux_unit='W/(m2 um)').to(
                    unit, sbu.spectral_density_vega(wfb))
            else:
                fluxd = obs.sample_binned(flux_unit=unit)

        if np.ndim(fluxd) != ndim:
            # likely need a squeeze
            fluxd = fluxd.squeeze()

        return fluxd

    def filt(self, bp, unit='W / (m2 um)', **kwargs):
        """Observe source through a single filter.


        Parameters
        ----------
        bp : string or `~synphot.SpectralElement`
            The name of a filter, or a transmission spectrum as a
            `~synphot.SpectralElement`.  See Notes for built-in filter
            names.

        unit : string, `~astropy.units.Unit`, optional
            Spectral flux density units of the output.

        **kwargs
            Additional keyword arguments for
            `~synphot.observation.Observation`, e.g., ``force``.


        Returns
        -------
        wave : `~astropy.units.Quantity`
            Effective wavelength.

        fluxd : `~astropy.units.Quantity` or float
            Spectral flux density.


        Notes
        -----

        Filter reference data is from STScI's Calibration Reference
        Data System.

        * ``'bessel_j'`` (Bessel *J*)
        * ``'bessel_h'`` (Bessel *H*)
        * ``'bessel_k'`` (Bessel *K*)
        * ``'cousins_r'`` (Cousins *R*)
        * ``'cousins_i'`` (Cousins *I*)
        * ``'johnson_u'`` (Johnson *U*)
        * ``'johnson_b'`` (Johnson *B*)
        * ``'johnson_v'`` (Johnson *V*)
        * ``'johnson_r'`` (Johnson *R*)
        * ``'johnson_i'`` (Johnson *I*)
        * ``'johnson_j'`` (Johnson *J*)
        * ``'johnson_k'`` (Johnson *K*)

        """

        import synphot
        from .. import units as sbu  # avoid circular dependency

        if isinstance(bp, str):
            bp = synphot.SpectralElement.from_filter(bp)

        obs = synphot.Observation(self.source, bp, **kwargs)
        wave = obs.effective_wavelength()
        _unit = u.Unit(unit)

        if _unit.is_equivalent(sbu.VEGAmag):
            fluxd = obs.effstim('W/(m2 um)').to(
                _unit, sbu.spectral_density_vega(bp))
        else:
            fluxd = obs.effstim(flux_unit=_unit)

        return wave, fluxd

    def color_index(self, wfb, unit, equivalencies=[], equiv_func=None):
        """Color index (magnitudes) and effective wavelengths.


        Parameters
        ----------
        wfb : two-element `~astropy.units.Quantity` or tuple
            Wavelengths, frequencies, or bandpasses of the
            measurement.  See :func:`~observe`.

        unit : string or `~astropy.units.MagUnit`
            Units for the output, e.g., ``astropy.units.ABmag`` or
            ``sbpy.units.VEGAmag``.

        equivalencies : list, optional
            List of unit equivalencies for magnitude-flux-density
            conversion.

        equiv_func : callable, optional
            A function that takes the effective wavelength of the
            observation and generates an equivalency list, e.g., use
            :func:`~sbpy.units.spectral_density_vega` for Vega-based
            magnitude conversions.
            :func:`~sbpy.units.spectral_density` is automatically
            included.


        Returns
        -------
        eff_wave : `~astropy.units.Quantity`
            Effective wavelengths for each ``wfb``.

        ci : `~astropy.units.Quantity`
            Color index, ``m_0 - m_1``, where 0 and 1 are element
            indexes for ``wfb``.

        """

        eff_wave = np.zeros(2) * u.um
        m = np.zeros(2) * u.Unit(unit)
        for i in range(2):
            if isinstance(wfb[i], u.Quantity):
                eff_wave[i] = wfb[i].to(u.um, u.spectral())
                f = self(eff_wave[i])
            else:
                eff_wave[i], f = self.filt(wfb[i])

            eqv = u.spectral_density(eff_wave[i])
            if equiv_func is not None:
                eqv.extend(equiv_func(eff_wave[i]))

            m[i] = f.to(m.unit, equivalencies + eqv)
        ci = m[0] - m[1]

        return eff_wave, ci


class SpectralStandard(SpectralSource, ABC):
    """Abstract base class for SBPy spectral standards.

    Parameters
    ----------
    source : `~synphot.SourceSpectrum`
        The source spectrum.

    description : string, optional
        A brief description of the source spectrum.

    bibcode : string, optional
        Bibliography code for `sbpy.bib.register`.


    Attributes
    ----------
    wave        - Wavelengths of the source spectrum.
    fluxd       - Source spectrum.
    description - Brief description of the source spectrum.
    meta        - Meta data from `source`, if any.

    """

    def __init__(self, source, description=None, bibcode=None):
        super().__init__(source, description=description)
        self._bibcode = bibcode

    @property
    def source(self):
        if self._bibcode is not None:
            register('spectroscopy', {self._description: self._bibcode})
        return self._source
