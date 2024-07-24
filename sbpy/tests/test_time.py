import numpy as np
from astropy.time import Time
from ..time import SpiceEphemerisTime


def test_spice_ephemeris_time():
    """Compare to SPICE result.

    spiceypy.utc2et("2022-08-01")
    712584069.1832777

    """

    t = Time("2022-08-01", scale="utc")
    assert np.isclose(t.et, 712584069.1832777, atol=1e-4, rtol=1e-14)
