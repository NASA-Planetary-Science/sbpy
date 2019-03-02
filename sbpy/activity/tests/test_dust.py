# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
from astropy.utils.data import get_pkg_data_filename
import synphot
from ..dust import *


def test_phase_HalleyMarcus():
    assert np.isclose(phase_HalleyMarcus(0 * u.deg), 1.0)
    assert np.isclose(phase_HalleyMarcus(15 * u.deg), 5.8720e-01)
    assert np.isclose(phase_HalleyMarcus(14.5 * u.deg), 0.5959274462322928)


class TestAfrho:
    def test_init(self):
        afrho = Afrho(1000 * u.cm)
        assert afrho.value == 1000
        assert afrho.unit == u.cm

    def test_scaler_ops(self):
        afrho = Afrho(1000 * u.cm)
        afrho = afrho / 2
        assert afrho == 500 * u.cm

    def test_quantity_ops(self):
        afrho = Afrho(1000 * u.cm)
        v = afrho * 2 * u.cm
        assert v == 2000 * u.cm**2
        assert not isinstance(v, Afrho)

    def test_from_fluxd(self):
        """HST/WFC3 photometry of C/2013 A1 (Siding Spring) (Li et al. 2014).

        Li et al. (2013) quotes the Sun in F606W as 1730 W/m2/um, but
        confirmed with J.-Y. Li that 1707 W/m2/um is a better value

        Li et al. (2014) quotes 1660 cm, which is 1680 cm (3
        significant figures) after solar flux density
        revision.

        """
        fluxd = 4.68e-15 * u.W / u.m**2 / u.um
        aper = 5000 * u.km
        eph = {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}
        S = 1707 * u.W / u.m**2 / u.um
        afrho = Afrho.from_fluxd(None, fluxd, aper, eph, S=S)
        assert np.isclose(afrho.cm, 1680, atol=4)

    def test_from_flam_with_synphot(self):
        wave = 0.55 * u.um
        fluxd = 6.764172537310662e-14 * u.W / u.m**2 / u.um
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        afrho = Afrho.from_fluxd(wave, fluxd, aper, eph)
        assert np.isclose(afrho.cm, 1000)

    def test_from_fnu(self):
        fluxd = 6.161081515869728 * u.mJy
        nu = 2.998e14 / 11.7 * u.Hz
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        S = 1.711e14 * u.Jy
        afrho = Afrho.from_fluxd(nu, fluxd, aper, eph, S=S)
        assert np.isclose(afrho.cm, 1000.0)

    def test_from_filt(self):
        """HST/WFC3 photometry of C/2013 A1 (Siding Spring) (Li et al. 2014).

        Li et al. quotes 1660 cm, but we have revised it to 1683 (see
        test_from_fluxd).

        """

        fluxd = 4.68e-15 * u.W / u.m**2 / u.um
        aper = 5000 * u.km
        eph = {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}
        fn = get_pkg_data_filename(os.path.join(
            '..', '..', 'photometry', 'data',
            'wfc3_uvis_f606w_004_syn.fits'))
        bandpass = synphot.SpectralElement.from_file(fn)
        afrho = Afrho.from_filt(bandpass, fluxd, aper, eph)
        assert np.isclose(afrho.cm, 1680, atol=4)

    @pytest.mark.parametrize('filename, mag, unit, rho, afrho0, eph, unc', (
        ('wfc3_uvis_f606w_004_syn.fits', 16.98, 'vegamag', '5000 km', 1680,
         {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}, 0.05),
        ('wfc3_uvis_f438w_004_syn.fits', 17.91, 'vegamag', '5000 km', 1550,
         {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}, 0.05),
        ('cousins_i_004_syn.fits', 8.49 - 0.53, 'vegamag', '10000 km', 3188,
         {'rh': 1.45 * u.au, 'delta': 0.49 * u.au}, 0.06),
        ('sdss-r.fits', 11.97, 'ABmag', '19.2 arcsec', 34.9,
         {'rh': 1.098 * u.au, 'delta': 0.164 * u.au}, 0.03),
        ('sdss-r.fits', 12.23, 'STmag', '19.2 arcsec', 34.9,
         {'rh': 1.098 * u.au, 'delta': 0.164 * u.au}, 0.03),
    ))
    def test_from_mag_bandpass(self, filename, mag, unit, rho, afrho0,
                               eph, unc):
        """Magnitude to afrho conversions.

        HST/WFC3 photometry of C/2013 A1 (Siding Spring) (Li et
        al. 2014).  Uncertainty is 5%.

        Woodward et al. photometry of C/2007 N3 (Lulin) in I-band.
        Their magnitude has been modifιed from 8.49 to 7.97 according
        to their phase correction (0.03 mag/deg, phase angle 17.77
        deg).  Uncertainty is 0.06 mag.

        Li et al. (2017) photometry of 252P/LINEAR in r'.  They use
        ABmag.  An additional test is performed with this observation
        after converting ABmag to STmag::

            synphot.units.convert_flux(6182, 11.97 * u.ABmag, u.STmag)

        """
        rho = u.Quantity(rho)
        fn = get_pkg_data_filename(os.path.join(
            '..', '..', 'photometry', 'data', filename))
        bandpass = synphot.SpectralElement.from_file(fn)
        afrho = Afrho.from_mag(mag, unit, rho, eph, bandpass=bandpass)
        assert np.isclose(afrho.cm, afrho0, rtol=unc)

    def test_from_mag_m_sun(self):
        """Verify Li et al. (2017) Afρ of 252P/LINEAR."""
        eph = {'rh': 1.098 * u.au, 'delta': 0.164 * u.au}
        afrho = Afrho.from_mag(11.97, 'ignored', 19.2 * u.arcsec, eph,
                               m_sun=-26.93)
        assert np.isclose(afrho.cm, 34.9, rtol=0.005)

    def test_fluxd(self):
        afrho = Afrho(1000, 'cm')
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        S = 1869 * u.W / u.m**2 / u.um
        fluxd = afrho.fluxd(None, aper, eph, S=S)
        assert fluxd.unit.is_equivalent(S.unit)
        assert np.isclose(fluxd.value, 6.730018324465526e-14)

    def test_fluxd_with_wave(self):
        afrho = Afrho(1000, 'cm')
        wave = 1 * u.um
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        unit = u.W / u.m**2 / u.um
        fluxd = afrho.fluxd(wave, aper, eph, unit=unit)
        assert np.isclose(fluxd.value, 2.6930875895493665e-14)

    def test_fluxd_with_freq(self):
        afrho = Afrho(1000, 'cm')
        freq = 299792.458 * u.GHz
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        unit = u.W / u.m**2 / u.um
        fluxd = afrho.fluxd(freq, aper, eph, unit=unit)
        assert np.isclose(fluxd.value, 2.6930875895493665e-14)

    @pytest.mark.parametrize('filename, mag0, unit, rho, afrho, eph, unc', (
        ('wfc3_uvis_f606w_004_syn.fits', 16.98, 'vegamag', '5000 km', 1680,
         {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}, 0.05),
        ('wfc3_uvis_f438w_004_syn.fits', 17.91, 'vegamag', '5000 km', 1550,
         {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}, 0.05),
        ('cousins_i_004_syn.fits', 8.49 - 0.53, 'vegamag', '10000 km', 3188,
         {'rh': 1.45 * u.au, 'delta': 0.49 * u.au}, 0.06),
        ('sdss-r.fits', 11.97, 'ABmag', '19.2 arcsec', 34.9,
         {'rh': 1.098 * u.au, 'delta': 0.164 * u.au}, 0.03),
        ('sdss-r.fits', 12.23, 'STmag', '19.2 arcsec', 34.9,
         {'rh': 1.098 * u.au, 'delta': 0.164 * u.au}, 0.03),
    ))
    def test_mag_bandpass(self, filename, mag0, unit, rho, afrho, eph, unc):
        """Inverse of test_from_mag_bandpass."""
        rho = u.Quantity(rho)
        fn = get_pkg_data_filename(os.path.join(
            '..', '..', 'photometry', 'data', filename))
        bandpass = synphot.SpectralElement.from_file(fn)
        mag = Afrho(afrho * u.cm).mag(unit, rho, eph, bandpass=bandpass)
        assert np.isclose(mag, mag0, rtol=unc)

    def test_mag_m_sun(self):
        """Inverse of test_from_mag_m_sun."""
        eph = {'rh': 1.098 * u.au, 'delta': 0.164 * u.au}
        mag = Afrho(34.9 * u.cm).mag('ignored', 19.2 * u.arcsec, eph,
                                     m_sun=-26.93)
        assert np.isclose(mag, 11.97, rtol=0.0005)

    def test_to_phase(self):
        afrho = Afrho(10 * u.cm).to_phase(15 * u.deg, 0 * u.deg)
        assert np.isclose(afrho.cm, 5.8720)


class TestEfrho:
    def test_init(self):
        efrho = Efrho(1000 * u.cm)
        assert efrho.value == 1000
        assert efrho.unit == u.cm

    def test_scaler_ops(self):
        efrho = Efrho(1000 * u.cm)
        efrho = efrho / 2
        assert efrho == 500 * u.cm

    def test_quantity_ops(self):
        efrho = Efrho(1000 * u.cm)
        v = efrho * 2 * u.cm
        assert v == 2000 * u.cm**2
        assert not isinstance(v, Efrho)

    def test_from_flam(self):
        fluxd = 3.824064455850402e-15 * u.W / u.m**2 / u.um
        wave = 10 * u.um
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        efrho = Efrho.from_fluxd(wave, fluxd, aper, eph)
        assert np.isclose(efrho.cm, 1000)

    def test_from_fnu(self):
        fluxd = 6.0961896974549115 * u.mJy
        nu = 2.998e14 / 11.7 * u.Hz
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        S = 1.711e14 * u.Jy
        efrho = Efrho.from_fluxd(nu, fluxd, aper, eph)
        assert np.isclose(efrho.cm, 33.0)

    def test_fluxd(self):
        efrho = Efrho(260.76955995377915, 'cm')
        wave = 5 * u.um
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        fluxd = efrho.fluxd(wave, aper, eph)
        assert np.isclose(fluxd.value, 1e-16)

    def test_fluxd_unit(self):
        efrho = Efrho(100, 'cm')
        wave = 5 * u.um
        aper = 1 * u.arcsec
        eph = dict(rh=1.5 * u.au, delta=1.0 * u.au)
        Tscale = 1.1
        fluxd = efrho.fluxd(wave, aper, eph, Tscale=Tscale, unit='mJy')
        assert np.isclose(fluxd.value, 0.3197891693353106)

    def test_from_mag_bandpass_fluxd0_error(self):
        with pytest.raises(TypeError):
            Efrho.from_mag(5, 'vegamag', 1 * u.arcsec,
                           dict(rh=1.5 * u.au, delta=1.0 * u.au))

    def test_from_mag_unit_error(self):
        bp = synphot.SpectralElement(
            synphot.Box1D, x_0=11.7 * u.um, width=0.1 * u.um)
        with pytest.raises(ValueError):
            Efrho.from_mag(5, 'asdf', 1 * u.arcsec, bp,
                           dict(rh=1.5 * u.au, delta=1.0 * u.au))

    @pytest.mark.parametrize('unit, efrho0', (
        ('vegamag', 616.1),
        ('abmag', 78750),  # compare with test_from_mag_fluxd0_B
        ('stmag', 3.596e7),
    ))
    def test_from_mag_bandpass(self, unit, efrho0):
        aper = 1 * u.arcsec
        # width = 0.1 um for speed
        bp = synphot.SpectralElement(
            synphot.Box1D, x_0=11.7 * u.um, width=0.1 * u.um)
        eph = dict(rh=1.0 * u.au, delta=1.0 * u.au)
        Tscale = 1.1
        efrho = Efrho.from_mag(5, unit, aper, eph, bandpass=bp,
                               Tscale=Tscale)
        assert np.isclose(efrho.cm, efrho0, rtol=0.001)

    def test_from_mag_fluxd0_bandpass(self):
        # comapre with test_from_mag_bandpass
        aper = 1 * u.arcsec
        # width = 0.1 um for speed
        bp = synphot.SpectralElement(
            synphot.Box1D, x_0=11.7 * u.um, width=0.1 * u.um)
        eph = dict(rh=1.0 * u.au, delta=1.0 * u.au)
        Tscale = 1.1
        fluxd0 = u.Quantity(3631, 'Jy')
        efrho = Efrho.from_mag(5, None, aper, eph, bandpass=bp,
                               fluxd0=fluxd0, Tscale=Tscale)
        assert np.isclose(efrho.cm, 78750, rtol=0.001)

    def test_from_mag_fluxd0_B(self):
        # comapre with test_from_mag_bandpass
        from astropy.modeling.blackbody import blackbody_nu
        aper = 1 * u.arcsec
        eph = dict(rh=1.0 * u.au, delta=1.0 * u.au)
        fluxd0 = u.Quantity(3631, 'Jy')
        Tscale = 1.1
        B = blackbody_nu(11.7 * u.um, 278 * Tscale * u.K)
        efrho = Efrho.from_mag(5, None, aper, eph, B=B, fluxd0=fluxd0)
        assert np.isclose(efrho.cm, 78750, rtol=0.001)

    @pytest.mark.parametrize('unit, efrho0', (
        ('vegamag', 616.1),
        ('abmag', 78750),  # compare with test_from_mag_fluxd0_B
        ('stmag', 3.596e7),
    ))
    def test_mag_bandpass(self, unit, efrho0):
        aper = 1 * u.arcsec
        # width = 0.1 um for speed
        bp = synphot.SpectralElement(
            synphot.Box1D, x_0=11.7 * u.um, width=0.1 * u.um)
        eph = dict(rh=1.0 * u.au, delta=1.0 * u.au)
        Tscale = 1.1
        efrho = Efrho(efrho0, 'cm')
        mag = efrho.mag(unit, aper, eph, bandpass=bp, Tscale=Tscale)
        assert np.isclose(mag, 5, rtol=0.001)

    def test_mag_fluxd0_bandpass(self):
        # comapre with test_mag_bandpass
        aper = 1 * u.arcsec
        # width = 0.1 um for speed
        bp = synphot.SpectralElement(
            synphot.Box1D, x_0=11.7 * u.um, width=0.1 * u.um)
        eph = dict(rh=1.0 * u.au, delta=1.0 * u.au)
        Tscale = 1.1
        fluxd0 = u.Quantity(3631, 'Jy')
        efrho = Efrho(78750, 'cm')
        mag = efrho.mag(None, aper, eph, bandpass=bp,
                        fluxd0=fluxd0, Tscale=Tscale)
        assert np.isclose(mag, 5, rtol=0.001)

    def test_mag_fluxd0_B(self):
        # comapre with test_mag_bandpass
        from astropy.modeling.blackbody import blackbody_nu
        aper = 1 * u.arcsec
        eph = dict(rh=1.0 * u.au, delta=1.0 * u.au)
        Tscale = 1.1
        fluxd0 = u.Quantity(3631, 'Jy')
        B = blackbody_nu(11.7 * u.um, 278 * Tscale * u.K)
        efrho = Efrho(78750, 'cm')
        mag = efrho.mag(None, aper, eph, B=B, fluxd0=fluxd0)
        assert np.isclose(mag, 5, rtol=0.001)
