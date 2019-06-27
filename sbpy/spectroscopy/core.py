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
from astropy import log
from ..bib import register
from ..data import Phys
from .sun import Sun


conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

__all__ = ['Spectrum', 'SpectralModel', 'SpectralGradient', 'molecular_data',
           'einstein_coeff', 'intensity_conversion']


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

    quantities = [t_freq, temp, lgint, part300,
                  partition, gu, energy_J, elo_J, df]
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


class SpectralGradient(u.SpecificTypeQuantity):
    r"""Convert between magnitude and spectral gradient.

    ``SpectralGradient`` is a `~astropy.units.Quantity` object with
    units of percent change per wavelength.

    Spectral gradient is with respect to the flux density at a
    particular wavelength.  By convention this is the mean flux
    between the two filters, assuming the reflectance spectrum is
    linear with wavelength.  Eq. 1 of A'Hearn et al. [ADT84]_:

    .. math::

        S = \frac{R(λ1) - R(λ0)}{R(λ1) + R(λ0)} \frac{2}{Δλ}

    Δλ is typically measured in units of 100 nm.


    Parameters
    ----------
    value : number, `~astropy.units.Quantity`
        The value(s).

    unit : string, `~astropy.units.Unit`, optional
        The unit of the input value.  Strings must be parseable by
        :mod:`~astropy.units` package.

    wave : `~astropy.units.Quantity`, optional
        Effective wavelengths of the measurement for a solar spectral
        energy distribution.  Required for conversion to other
        wavelengths.

    wave0 : `~astropy.units.Quantity`, optional
        Normalization point.  If ``None``, the mean of ``wave``
        will be assumed.

    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`.

    copy : bool, optional
        See `~astropy.units.Quantity`.


    Examples
    --------
    >>> import astropy.units as u
    >>> from sbpy.units import hundred_nm
    >>> S = SpectralGradient(10 * u.percent / hundred_nm, wave0=5500 * u.AA)
    >>> print(S)
    10.0 % / 100 nm

    >>> from sbpy.units import VEGAmag
    >>> bp = ('johnson_v', 'cousins_r')
    >>> VmR = 15.8 * VEGAmag - 15.3 * VEGAmag
    >>> VmR_sun = 0.37 * u.mag
    >>> S = SpectralGradient.from_color(bp, VmR - VmR_sun)
    ...                              # doctest: +REMOTE_DATA +IGNORE_OUTPUT
    >>> print(S)                     # doctest: +REMOTE_DATA +FLOAT_CMP
    12.29185986266534 % / 100 nm


    References
    ----------
    .. [ADT84] A'Hearn, Dwek & Tokunaga 1984. Infrared Photometry of
        Comet Bowell and Other Comets. ApJ 282, 803-806.

    """

    _equivalent_unit = u.meter**-1
    _include_easy_conversion_members = False

    def __new__(cls, value, unit=None, wave=None, wave0=None,
                dtype=None, copy=None):
        S = super().__new__(cls, value, unit=unit, dtype=dtype,
                            copy=copy)

        if wave is not None:
            if np.size(wave) != 2:
                raise ValueError(
                    'Two wavelengths must be provided, got {}'
                    .format(np.size(wave)))
            S.wave = S._eff_wave(wave)

        if wave0 is None and wave is not None:
            S.wave0 = S.wave.mean()
        else:
            S.wave0 = wave0

        return S

    @classmethod
    def _eff_wave(cls, wfb):
        """Wavelength/frequency/bandpass to wavelength.

        Bandpass is converted to effective wavelength using a solar
        spectrum.

        """

        eff_wave = (0, 0) * u.um
        sun = Sun.from_default()
        for i in range(2):
            if isinstance(wfb, u.Quantity):
                eff_wave[i] = wfb[i].to(u.um)
            else:
                eff_wave[i] = sun.filt(wfb[i])[0]
                log.info('Using λ_eff = {}'.format(eff_wave[i]))

        return eff_wave

    @classmethod
    def from_color(cls, wfb, color):
        r"""Initialize from observed color.


        Parameters
        ----------
        wfb : two-element `~astropy.units.Quantity` or tuple
            Wavelengths, frequencies, or bandpasses of the
            measurement.  If a bandpass, the effective wavelength of a
            solar spectrum will be used.  Bandpasses may be a string
            (name) or `~synphot.SpectralElement` (see
            :func:`~sbpy.spectroscopy.sun.Sun.filt`).

        color : `~astropy.units.Quantity`, optional
            Observed color, ``blue - red`` for magnitudes, ``red
            / blue`` for linear units.  Must be dimensionless and have
            the solar color removed.


        Notes
        -----

        Computes spectral gradient from ``color_index``.
        ``wfb[0]`` is the blue-ward of the two measurements

        .. math::

           S &= \frac{R(λ1) - R(λ0)}{R(λ1) + R(λ0)} \frac{2}{Δλ} \\
           &= \frac{α - 1}{α + 1} \frac{2}{Δλ}

        where R(λ) is the reflectivity, and:

        .. math::

            α = R(λ1) / R(λ0) = 10^{0.4 color_index}

            color_index = Δm - C_{sun}

        Δλ is typically expressed in units of 100 nm.


        Examples
        --------
        >>> import astropy.units as u
        >>> w = [0.4719, 0.6185] * u.um
        >>> S = SpectralGradient.from_color(w, 0.10 * u.mag)
        >>> print(S)                            # doctest: +FLOAT_CMP
        6.27819572 % / 100 nm

        """
        from ..units import hundred_nm

        eff_wave = SpectralGradient._eff_wave(wfb)

        try:
            # works for u.Magnitudes and dimensionless u.Quantity
            alpha = u.Quantity(color, u.dimensionless_unscaled)
        except u.UnitConversionError:
            # works for u.mag
            alpha = color.to(u.dimensionless_unscaled, u.logarithmic())

        dw = eff_wave[0] - eff_wave[1]
        S = ((2 / dw * (alpha - 1) / (alpha + 1))
             .to(u.percent / hundred_nm))

        return SpectralGradient(S, wave=eff_wave)

    def to_color(self, wfb):
        r"""Express as a color index.


        Parameters
        ----------
        wfb : two-element `~astropy.units.Quantity` or tuple
            Wavelengths, frequencies, or bandpasses of the
            measurement.  If a bandpass, the effective wavelength of a
            solar spectrum will be used.  Bandpasses may be a string
            (name) or `~synphot.SpectralElement` (see
            :func:`~sbpy.spectroscopy.sun.Sun.filt`).

        Notes
        -----
        Color index is computed from:

        .. math::

            α = \frac{1 + S Δλ / 2}{1 - S * Δλ / 2}

        where S is the spectral gradient at the mean of λ0 and λ1, and:

        .. math::

            α = R(λ1) / R(λ0) = 10^{0.4 color_index}

            color_index = Δm - C_{sun}

        Δλ is typically expressed in units of 100 nm.


        Returns
        -------
        color : `~astropy.units.Quantity`
            ``blue - red`` color in magnitudes, dimensionless and
            excludes the solar color.


        Examples
        --------
        >>> import astropy.units as u
        >>> from sbpy.units import hundred_nm
        >>> S = SpectralGradient(10 * u.percent / hundred_nm,
        ...                      wave0=0.55 * u.um)
        >>> C = S.to_color((525, 575) * u.nm)
        >>> print(C)    # doctest: +FLOAT_CMP
        0.05429812423309064 mag

        """

        eff_wave = self._eff_wave(wfb)

        S = self.renormalize(eff_wave.mean())
        dw = eff_wave[0] - eff_wave[1]
        beta = (S * dw / 2).decompose()  # dimensionless
        color = ((1 + beta) / (1 - beta)).to(u.mag, u.logarithmic())

        return color

    def renormalize(self, wave0):
        """Re-normalize to another wavelength.

        The slope is linearly extrapolated to the new normalization
        point.  Requires the `wave0` attribute to be defined, see
        `~SpectralGradient`.


        Parameters
        ----------
        wave0 : `~astropy.units.Quantity`
            Wavelength.


        Returns
        -------
        S : ``SpectralGradient``
            Renormalized gradient.


        Examples
        --------
        >>> import astropy.units as u
        >>> from sbpy.spectroscopy import SpectralGradient
        >>> from sbpy.units import hundred_nm
        >>> S1 = SpectralGradient(10 * u.percent / hundred_nm,
        ...                      wave0=0.55 * u.um)
        >>> S2 = S1.renormalize(3.6 * u.um)
        >>> print(S2.renormalize(3.6 * u.um))    # doctest: +FLOAT_CMP
        2.469135802469136 % / 100 nm
        >>> print(S2.wave0)
        3.6 um
        """

        if self.wave0 is None:
            raise ValueError('wave0 attribute must be defined.')

        delta = wave0 - self.wave0
        S0 = 1 + self.to(delta.unit**-1) * delta
        S = self / S0
        S.wave0 = wave0
        return S


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
