import pytest
import numpy as np
import astropy.units as u
from ..core import *
from .. import sources

def test_default_sun_string():
    with default_sun.set('E490_2014'):
        assert default_sun.get().description == sources.E490_2014['description']
    
    with default_sun.set('E490_2014LR'):
        assert default_sun.get().description == sources.E490_2014LR['description']

@pytest.mark.remote_data
def test_default_sun_string_remote():
    with default_sun.set('Kurucz1993'):
        assert default_sun.get().description == sources.Kurucz1993['description']
    
    with default_sun.set('Castelli1996'):
        assert default_sun.get().description == sources.Castelli1996['description']

def test_sun_rebin_single():
    with default_sun.set('E490_2014'):
        f = default_sun.get().rebin(0.5555 * u.um, unit='W/(m2 um)')
        assert np.isclose(f.value, 1897)

    with default_sun.set('E490_2014LR'):
        f = default_sun.get().rebin(0.55 * u.um, unit='W/(m2 um)')
        assert np.isclose(f.value, 1878.485)

