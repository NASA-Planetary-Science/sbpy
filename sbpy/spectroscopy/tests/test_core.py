# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import mock
import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
import synphot
from ..core import SpectralStandard
from ... import bib, units, utils


class Star(SpectralStandard):
    pass


class TestSpectralStandard:
    def test_init_importerror(self):
        with mock.patch.dict(sys.modules, {'synphot': None}):
            with pytest.raises(ImportError):
                s = Star(None)

    def test_from_array(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)), 'W/(m2 um)')
        s = Star.from_array(w, f)
        assert np.allclose(f.value, s(w, f.unit).value)

    def test_from_array_importerror(self):

        with mock.patch.dict(sys.modules, {'synphot': None}):
            with pytest.raises(ImportError):
                s = Star.from_array(None, None)

    # from_file tested by Sun and Vega

    def test_from_file_importerror(self):
        with mock.patch.dict(sys.modules, {'synphot': None}):
            with pytest.raises(ImportError):
                s = Star.from_file(None)

    def test_description(self):
        s = Star(None, description='asdf')
        assert s.description == 'asdf'

    def test_wave(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)), 'W/(m2 um)')
        source = synphot.SourceSpectrum(synphot.Empirical1D, points=w,
                                        lookup_table=f)
        s = Star(source)
        assert np.allclose(s.wave.to(w.unit).value, w.value)

    def test_fluxd(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)), 'W/(m2 um)')
        source = synphot.SourceSpectrum(synphot.Empirical1D, points=w,
                                        lookup_table=f)
        s = Star(source)
        assert np.allclose(s.fluxd.to(f.unit).value, f.value)

    # meta tested by Sun

    def test_call_wavelength(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(0.5 * w.value + 0.1, 'W/(m2 um)')
        s = Star.from_array(w, f)
        w = np.linspace(0.31, 0.99) * u.um
        assert np.allclose(s(w).value, 0.5 * w.value + 0.1, rtol=0.0005)

    def test_call_frequency(self):
        nu = u.Quantity(np.linspace(300, 1000), 'THz')
        f = u.Quantity(0.5 * nu.value + 0.1, 'Jy')
        s = Star.from_array(nu, f)
        nu = np.linspace(310, 999) * u.THz
        assert np.allclose(s(nu).value, 0.5 * nu.value + 0.1, rtol=0.002)

    def test_observe_wavelength(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)) * 0.35 / w.value, 'W/(m2 um)')
        s = Star.from_array(w, f)
        w = [0.3, 0.35, 0.4] * u.um
        assert np.isclose(s.observe(w).value[1], 1.0, rtol=0.001)

    def test_observe_frequency(self):
        nu = u.Quantity(np.linspace(300, 1000), 'THz')
        f = u.Quantity(np.ones(len(nu)) * nu.value / 350, 'Jy')
        s = Star.from_array(nu, f)
        nu = [325, 350, 375] * u.THz
        assert np.isclose(s.observe(nu).value[1], 1.0, rtol=0.004)

    def test_observe_bandpass(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)), 'W/(m2 um)')
        s = Star.from_array(w, f)

        bp = synphot.SpectralElement(synphot.Box1D, x_0=0.55 * u.um,
                                     width=0.1 * u.um)
        fluxd = s.observe(bp)
        assert np.allclose(fluxd.value, [1])

    def test_bibcode(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)) * 0.35 / w.value, 'W/(m2 um)')
        s = Star.from_array(w, f, bibcode='asdf', description='fdsa')

        with bib.Tracking():
            s.source

        assert 'asdf' in bib.to_text()

    @remote_data
    def test_filt_str(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)), 'W/(m2 um)')
        s = Star.from_array(w, f)
        fluxd = s.filt('johnson_v', unit='W/(m2 um)')[1]
        assert np.isclose(fluxd.value, 1)

    def test_filt_vegamag(self):
        w = u.Quantity(np.linspace(0.3, 1.0), 'um')
        f = u.Quantity(np.ones(len(w)), 'W/(m2 um)')
        s = Star.from_array(w, f)
        V = utils.get_bandpass('johnson v')
        mag = s.filt(V, unit=units.VEGAmag)[1]
        # -18.60 is -2.5 * log10(3636e-11)
        assert np.isclose(mag.value, -18.60, atol=0.02)
