# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Spectroscopy Module

created on June 23, 2017
"""

from abc import ABC
import numpy as np
import astropy.constants as con
import astropy.units as u
from astroquery.jplspec import JPLSpec
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from ..activity.gas import photo_timescale
from astroquery.jplhorizons import conf

conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

__all__ = ['Spectrum', 'SpectralModel']


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
            | Planck Constant
            | Speed of light
            | Boltzmann constant

    """

    query = JPLSpec.query_lines(min_frequency=(transition_freq - (1 * u.MHz)),
                                max_frequency=(transition_freq + (1 * u.MHz)),
                                molecule=mol_tag)

    freq_list = query['FREQ']

    t_freq = min(list(freq_list.quantity), key=lambda x: abs(x-transition_freq))

    data = query[query['FREQ'] == t_freq.value]

    df = int(data['DR'].data)

    lgint = float(data['LGINT'].data)

    lgint = 10**lgint * u.nm * u.nm * u.MHz

    elo = float(data['ELO'].data) / u.cm

    gu = float(data['GUP'].data)

    cat = JPLSpec.get_species_table()

    mol = cat[cat['TAG'] == mol_tag]

    temp_list = cat.meta['Temperature (K)'] * u.K

    temp = min(temp_list, key=lambda x: abs(x-temp_estimate))

    for i in range(1, 8):
        if mol.columns['QLOG{}'.format(i)].meta !=\
           {'Temperature (K)': temp.value}:

            pass

        else:
            partition = 10 ** (float(mol.columns['QLOG{}'.format(i)].data))

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

    if ((energy_J-elo_J).value/(k*temp).value) and ((energy_J-elo_J).value /
                                                    (k*300 * u.K).value) < 1:

        if df == 0 or 2:

            intl = lgint*(300*u.K / temp)**(2)*np.exp(-(1/temp - 1/(300*u.K))
                                                      * elo_J/k)

        else:

            intl = lgint*(300*u.K/temp)**(5/2)*np.exp(-(1/temp - 1/(300*u.K))
                                                      * elo_J/k)

    else:

        intl = lgint*(part300/partition)*(np.exp(-elo_J/(k*temp)) -
                                          np.exp(-energy_J/k*temp)) / \
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
              (partition/gu)*(np.exp(-elo_J/(k*temp)) -
                              np.exp(-energy_J/(k*temp)))**(-1)
              * (2.7964e-16)).value

    return au / u.s


def photod_rate(time, time_scale, target, id_type, observatory, format,
                mol_tag):

    epoch = Time(time, scale=time_scale, format=format)
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

    if type(timescale) == np.array:
        timescale = timescale[0]

    beta = (timescale) * r**2

    return beta, delta


def total_number(integrated_flux, frequency):
    """
    Basic equation relating number of molecules with observed integrated flux.
    This is given by equation 10 in
    https://ui.adsabs.harvard.edu/#abs/2004come.book..523C

    Parameters
    ----------
    integrated_flux : `~astropy.units.Quantity`
        Integrated flux of emission line.
    frequency : `~astropy.units.Quantity`
        Transition frequency

    Returns
    -------
    total_number : float
        Total number of molecules within the aperture

    not implemented
    """
    total_number = integrated_flux
    total_number *= 8*np.pi*u.k_B*frequency**2/(con.h*con.c**3 *
                                                einstein_coeff(frequency))
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
                    diameter=25 * u.m, observatory='500', b=1.2,
                    format='iso', time_scale='utc',
                    id_type='designation'):

        """
        | Returns production rate based on Drahus 2012 model referenced. Includes
        | no photodissociation

        Parameters
        ----------
        spectra : `~astropy.units.Quantity`
            Temperature brightness integral derived from spectra

        transition_freq : `~astropy.units.Quantity`
            Transition frequency in MHz

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

        diameter : `~astropy.units.Quantity`
            Telescope diameter in meters. Default is 25 m

        observatory : str
            | Observatory identifier as per `~astroquery.jplhorizons`
            | Default is geocentric ('500')

        b : int
            | Dimensionless factor intrinsic to every antenna. Typical
            | value, and the default for this model, is 1.22. See
            | references for more information on this parameter.

        format : str
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

        >>> diameter = 30 * u.m  # doctest: +SKIP

        >>> b = 1.13  # doctest: +SKIP

        >>> mol_tag = 27001  # doctest: +SKIP

        >>> transition_freq = 265.886434 * u.MHz  # doctest: +SKIP

        >>> spectra = 1.22 * u.K * u.km / u.s  # doctest: +SKIP

        >>> time = '2010-11-3 00:48:06'  # doctest: +SKIP

        >>> q = prodrate_np(spectra, temp_estimate, transition_freq, # doctest: +SKIP
                            mol_tag, time, target, vgas, diameter,
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

        assert isinstance(temp_estimate, u.Quantity)
        assert isinstance(vgas, u.Quantity)
        assert isinstance(transition_freq, u.Quantity)
        assert isinstance(diameter, u.Quantity)
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
                                  observatory, format, mol_tag)

        calc = ((16*np.pi*k*t_freq.decompose() *
                 partition*vgas) / ((np.sqrt(np.pi*np.log(2)))
                * h * c**2 * au * gu * np.exp(-energy_J/(k*temp)))).decompose()

        q = spectra*(calc * b * delta / diameter)

        q = q.decompose().to(u.Hz, equivalencies=u.spectral()).decompose()[0]

        return q

    def production_rate(self, coma, molecule, frequency, aper):
        """
        Calculate production rate for `GasComa`

        Parameters
        ----------
        coma : `sbpy.activity.gas.GasComa`
            Gas coma model

        Returns
        -------
        Q : `~astropy.units.Quantity`
            production rate

        Examples
        --------
        >>> from sbpy.activity.gas import Haser
        >>> coma = Haser(Q, v, parent) # doctest: +SKIP
        >>> Q = spec.production_rate(coma, molecule='H2O') # doctest: +SKIP

        not yet implemented

        """

        from ..activity.gas import GasComa

        assert isinstance(coma, GasComa)

        integrated_line = self.integrated_flux(frequency)

        molecules = total_number(integrated_line, frequency)

        model_molecules = coma.total_number(aper)

        Q = coma.q * molecules/model_molecules

        return Q

    def plot(self):
        """Plot spectrum

        Returns
        -------
        `matplotlib.pyplot` instance

        not yet implemented
        """


class SpectralStandard(ABC):
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
    meta        - Meta data from `source`.

    """

    def __init__(self, source, description=None, bibcode=None):
        self._source = source
        self._description = description
        self._bibcode = bibcode

    @classmethod
    def from_array(cls, wave, fluxd, meta=None, **kwargs):
        """Create standard from arrays.

        Parameters
        ----------
        wave : `~astropy.units.Quantity`
          The spectral wavelengths.

        fluxd : `~astropy.units.Quantity`
          The solar flux densities, at 1 au.

        meta : dict, optional
          Meta data.

        **kwargs
          Passed to object initialization.

        """

        import synphot

        source = synphot.SourceSpectrum(
            synphot.Empirical1D, points=wave, lookup_table=fluxd,
            meta=meta)

        return cls(source, **kwargs)

    @classmethod
    def from_file(cls, filename=None, wave_unit=None, flux_unit=None,
                  cache=True, **kwargs):
        """Load the source spectrum from a file.

        Parameters
        ----------
        filename : string
          The name of the file.  See
          `~synphot.SourceSpectrum.from_file` for details.

        wave_unit, flux_unit : str or `~astropy.units.core.Unit`, optional
          Wavelength and flux units.

        cache : bool, optional
          If `True`, cache the contents of URLs.

        **kwargs
          Passed to `Sun` initialization.

        """

        from astropy.utils.data import download_file
        from astropy.utils.data import _is_url
        import synphot
        from synphot.specio import read_fits_spec, read_ascii_spec

        assert filename is not None, "File name required."

        # URL cache because synphot.SourceSpectrum.from_file does not
        if _is_url(filename):
            if filename.lower().endswith(('.fits', '.fit')):
                read_spec = read_fits_spec
            else:
                read_spec = read_ascii_spec

            fn = download_file(filename, cache=True)
            spec = read_spec(fn, wave_unit=wave_unit, flux_unit=flux_unit)
            source = synphot.SourceSpectrum(
                synphot.Empirical1D, points=spec[1], lookup_table=spec[2],
                meta={'header': spec[0]})
        else:
            source = synphot.SourceSpectrum.from_file(
                filename, wave_unit=wave_unit, flux_unit=flux_unit)

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
        from .. import bib
        if self._bibcode is not None:
            bib.register('spectroscopy', {self._description, self._bibcode})
        return self._source

    @property
    def meta(self):
        self._source.meta

    def __call__(self, wave_or_freq, unit=None):
        """Evaluate the source spectrum.

        Parameters
        ----------
        wave_or_freq : `~astropy.units.Quantity`
          Requested wavelength or frequencies of the resulting
          spectrum.  If an array, `wave_or_freq` specifies bin
          centers.  If a single value, fluxd will be interpolated to
          this wavelength/frequency.

        unit : string, `~astropy.units.Unit`, optional
          Spectral units of the output: flux density, 'vegamag',
          'ABmag', or 'STmag'.  If ``None``, return units are W/(m2
          μm) for ``wave_or_freq`` as wavelength, otherwise return Jy.

        Returns
        -------
        fluxd : `~astropy.units.Quantity`
          The spectrum binned to match `wave_or_freq`.  If a single
          point is requested, the original spectrum is interpolated to
          it.

        """

        import numpy as np
        import synphot

        if unit is None:
            if wave_or_freq.unit.is_equivalent('m'):
                unit = u.Unit('W/(m2 um)')
            else:
                unit = u.Jy

        if np.size(wave_or_freq) > 1:
            # Method adapted from http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
            specele = synphot.SpectralElement(synphot.ConstFlux1D(1))
            obs = synphot.Observation(self.source, specele, binset=wave_or_freq,
                                      force='taper')

            # The following is the same as obs.binflux, except sample_binned
            # will do the unit coversions.
            fluxd = obs.sample_binned(flux_unit=unit)
        else:
            fluxd = self.source(wave_or_freq, unit)

        return fluxd

    def filt(self, bp, unit='W / (m2 um)', **kwargs):
        """Spectrum observed through a filter.

        Parameters
        ----------
        bp : string or `~synphot.SpectralElement`
          The name of a filter, or a transmission spectrum as a
          `~synphot.SpectralElement`.  See notes for built-in filter
          names.

        unit : string, `~astropy.units.Unit`, optional
          Spectral units of the output: flux density, 'vegamag',
          'ABmag', or 'STmag'.

        **kwargs
          Additional keyword arguments for
          `~synphot.observation.Observation`.

        Returns
        -------
        wave : `~astropy.units.Quantity`
          Effective wavelength.

        fluxd : `~astropy.units.Quantity` or float
          Flux density or magnitude.


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
        from synphot.units import VEGAMAG
        from .vega import Vega

        assert isinstance(bp, (str, synphot.SpectralElement))

        if isinstance(bp, str):
            bp = synphot.SpectralElement.from_filter(bp)

        obs = synphot.Observation(self.source, bp, **kwargs)
        wave = obs.effective_wavelength()

        if str(unit).lower().strip() == 'vegamag':
            f = obs.effstim('W/(m2 um)')
            f0 = Vega.from_default().filt(bp, unit='W/(m2 um)')[1]
            fluxd = -2.5 * np.log10((f / f0).value) * VEGAMAG
        else:
            fluxd = obs.effstim(unit)

        return wave, fluxd
