# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..gas import *

class TestHaser:
    def test_volume_density(self):
        # reproduce Newburn and Johnson 1978
        # Species, v (km/s), Q, n(1e4 km) cm**-3
        # CN, 0.71, 4.1e23, 0.5
        # C3, 0.59, 5.3e23, 1.1
        # C2, 0.73, 4.3e24, 6.1

        Q = [4.1e23, 5.3e23, 4.3e24] / u.s
        v = [0.71, 0.59, 0.73] * u.km / u.s
        parent = [1.4e4, 0, 1.0e4] * u.km / 1.07**2
        daughter = [1.7e5, 4.6e4, 7.6e4] * u.km / 1.07**2

        n = np.zeros(3)
        for i in range(3):
            coma = Haser(Q[i], v[i], parent[i], daughter[i])
            n[i] = coma.volume_density(1e4 * u.km).to(u.cm**-3).value

        assert np.allclose(n, [0.5, 1.1, 6.1])

    def test_column_density_small_aperture(self):
        # Test column desnity for aperture << lengthscale
        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        aper = 10 * u.km
        parent = 1e4 * u.km
        sigma = Haser(Q, v, parent).column_density(aper)
        ideal = Q / v / 2 / aper
        assert np.isclose(sigma.decompose().value, ideal.decompose().value)

    def test_column_density_large_aperture(self):
        # Test column desnity for aperture >> lengthscale
        Q = 1 / u.s
        v = 1 * u.km / u.s
        aper = 1000 * u.km
        parent = 10 * u.km
        sigma = Haser(Q, v, parent).column_density(aper)
        ideal = 4 * Q * parent / v / np.pi / 4 / aper**2
        assert np.isclose(sigma.decompose().value, ideal.decompose().value)

    def test_total_number(self):
        # reproduce Newburn and Johnson 1978
        # Species, N observed, Q/v (km**-1)
        # CN, 6.4e26, 5.8e23
        # C3, 8.3e28, 9.0e23
        # C2, 7.8e27, 5.9e24

        Nobs = np.array([6.41756750e26, 8.63191842e+28, 7.81278300e27])
        parent = [1.4e4, 0, 1.0e4] * u.km
        daughter = [1.7e5, 4.6e4, 7.6e4] * u.km
        Q = [5.8e23, 9.0e23, 5.9e24] / u.s
        rho = 3300 * u.km

        N = np.zeros(3)
        for i in range(3):
            coma = Haser(Q[i], 1 * u.km / u.s, parent[i], daughter[i])
            N[i] = coma.total_number(rho)

        assert np.allclose(N, Nobs)

