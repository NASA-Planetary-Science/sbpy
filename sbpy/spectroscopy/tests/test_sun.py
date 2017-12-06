import numpy as np
import astropy.units as u
from .. import sun

def test_default_sun_string():
    with sun.default_sun.set('E490_2014'):
        assert sun.default_sun.get() == sun.E490_2014
    
    with sun.default_sun.set('E490_2014LR'):
        assert sun.default_sun.get() == sun.E490_2014LR
    
def test_sun_rebin_single():
    with sun.default_sun.set('E490_2014'):
        f = sun.default_sun.get().rebin(0.5555 * u.um, unit='W/(m2 um)')
        assert np.isclose(f.value, 1897)

    with sun.default_sun.set('E490_2014LR'):
        f = sun.default_sun.get().rebin(0.55 * u.um, unit='W/(m2 um)')
        assert np.isclose(f.value, 1867)

