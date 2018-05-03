import pytest
import numpy as np
import astropy.units as u
from .. import *

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def test_default_sun_string():
    with default_sun.set('E490_2014'):
        assert default_sun.get().description == sources.E490_2014[
            'description']

    with default_sun.set('E490_2014LR'):
        assert default_sun.get().description == sources.E490_2014LR[
            'description']


@pytest.mark.remote_data
def test_default_sun_string_remote():
    with default_sun.set('Kurucz1993'):
        assert default_sun.get().description == sources.Kurucz1993[
            'description']

    with default_sun.set('Castelli1996'):
        assert default_sun.get().description == sources.Castelli1996[
            'description']


def test_sun_call_single():
    with default_sun.set('E490_2014'):
        sun = default_sun.get()
        f = sun(0.5555 * u.um, unit='W/(m2 um)')
        assert np.isclose(f.value, 1897)

    with default_sun.set('E490_2014LR'):
        sun = default_sun.get()
        f = sun(0.55 * u.um, unit='W/(m2 um)')
        assert np.isclose(f.value, 1878.485)


@pytest.mark.skipif('not HAS_SCIPY')
def test_sun_binning():
    from scipy.integrate import trapz

    # compare Sun's rebinning with an integration over the spectrum
    sun = Sun.from_builtin('E490_2014')

    wave0 = sun.wave.to('um').value
    fluxd0 = sun.fluxd.to('W/(m2 um)').value

    wave = np.linspace(0.35, 0.55, 6)

    d = np.diff(wave)[0] / 2
    left_bins = wave - d
    right_bins = wave + d

    fluxd1 = np.zeros(len(wave))
    for i in range(len(wave)):
        j = (wave0 >= left_bins[i]) * (wave0 <= right_bins[i])
        fluxd1[i] = trapz(fluxd0[j] * wave0[j], wave0[j]) / trapz(
            wave0[j], wave0[j])

    fluxd2 = sun(wave * u.um).value

    assert np.allclose(fluxd1, fluxd2, 0.005)
