# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import pytest
import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
from astropy.modeling.blackbody import blackbody_nu, blackbody_lambda
import synphot
from .. import sources
from ..sources import (BlackbodySource, SinglePointSpectrumError,
                       SpectralSource, SynphotRequired)
from ... import bib, units
from ...photometry import bandpass

V = bandpass('johnson v')
I = bandpass('cousins i')


class Star(SpectralSource):
    def __init__(self):
        super().__init__(synphot.SourceSpectrum(
            synphot.ConstFlux1D, amplitude=1 * u.W / u.m**2 / u.um))


class TestSpectralSource:
    def test_init_error(self, monkeypatch):
        with pytest.raises(SynphotRequired):
            monkeypatch.setattr(sources, 'synphot', None)
            Star()

    @pytest.mark.parametrize('wfb, interpolate', (
        ([V], False),
        ([1] * u.um, True),
        ([1, 2, 3] * u.um, False)
    ))
    def test_observe(self, wfb, interpolate):
        s = Star()
        fluxd = s.observe(wfb, unit='W/(m2 um)', interpolate=interpolate)
        assert np.isclose(fluxd.value, 1.0).all()

    def test_observe_bad_wfb(self):
        with pytest.raises(TypeError):
            s = Star()
            s.observe(np.arange(5))

    @pytest.mark.parametrize('wfb, unit, atol', (
        ((V, I), u.ABmag, 0.006),
        ((600 * u.nm, 750 * u.nm), u.ABmag, 1e-6),
        ((750 * u.GHz, 600 * u.GHz), u.ABmag, 1e-6),
    ))
    def test_color_index(self, wfb, unit, atol):
        s = Star()
        eff_wave, ci = s.color_index(wfb, unit)
        test = -5 * np.log10((eff_wave.min() / eff_wave.max()).value)
        assert np.isclose(ci.value, test, atol=atol)

    def test_color_index_typeerror(self):
        s = Star()
        with pytest.raises(TypeError):
            s.color_index(np.arange(2), unit=u.ABmag)


class TestBlackbodySource:
    @pytest.mark.parametrize('T', (
        300, 300 * u.K
    ))
    def test_init_temperature(self, T):
        BB = BlackbodySource(T)
        assert BB.T.value == 300

    def test_init_temperature_error(self):
        with pytest.raises(TypeError):
            BlackbodySource()

    def test_repr(self):
        BB = BlackbodySource(278)
        assert repr(BB) == '<BlackbodySource: T=278.0 K>'

    @pytest.mark.parametrize('B', (blackbody_nu, blackbody_lambda))
    def test_call(self, B):
        w = np.logspace(-0.5, 3) * u.um
        f = B(w, 300 * u.K) * np.pi * u.sr
        BB = BlackbodySource(300 * u.K)
        test = BB(w, unit=f.unit).value
        assert np.allclose(test, f.value)
