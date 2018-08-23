# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
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

        Throughput file for this test at:
        http://www.stsci.edu/hst/wfc3/ins_performance/throughputs/Throughput_Tables

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

<<<<<<< HEAD
=======
    def test_from_filt(self):
        """HST/WFC3 photometry of C/2013 A1 (Siding Spring) (Li et al. 2014).

        Li et al. quotes 1660 cm, but we have revised it to 1683 (see
        test_from_fluxd).

        Throughput file for this test at:
        http://www.stsci.edu/hst/wfc3/ins_performance/throughputs/Throughput_Tables

        """

        fluxd = 4.68e-15 * u.W / u.m**2 / u.um
        aper = 5000 * u.km
        eph = {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}
        fn = os.path.join(os.path.dirname(__file__), 'data',
                          'wfc3_uvis_f606w_004_syn.fits')
        bandpass = synphot.SpectralElement.from_file(fn)
        afrho = Afrho.from_filt(bandpass, fluxd, aper, eph)
        assert np.isclose(afrho.cm, 1680, atol=4)

    def test_from_mag(self):
        """HST/WFC3 photometry of C/2013 A1 (Siding Spring) (Li et al. 2014).

        Throughput files for this test:
        http://www.stsci.edu/hst/wfc3/ins_performance/throughputs/Throughput_Tables

        """
        mag = 16.98
        rho = 5000 * u.km
        eph = {'rh': 4.582 * u.au, 'delta': 4.042 * u.au}
        fn = os.path.join(os.path.dirname(__file__), 'data',
                          'wfc3_uvis_f606w_004_syn.fits')
        bandpass = synphot.SpectralElement.from_file(fn)
        afrho = Afrho.from_mag(bandpass, mag, 'vegamag', rho, eph)
        assert np.isclose(afrho.cm, 1680, atol=4)

>>>>>>> Revise dust tests.
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
