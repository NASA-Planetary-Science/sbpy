# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import astropy.units as u
from ..gas import *

def test_photo_lengthscale():
    gamma = photo_lengthscale('OH', 'CS93')
    assert gamma == 1.6e5 * u.km

def test_photo_timescale():
    tau = photo_timescale('CO2', 'CE83')
    assert tau == 5.0e5 * u.s

class TestHaser:
    def test_volume_density(self):
        """Reproduce Newburn and Johnson 1978.

        Species, v (km/s), Q, n(1e4 km) cm**-3
        CN, 0.71, 4.1e23, 0.5
        C3, 0.59, 5.3e23, 1.1
        C2, 0.73, 4.3e24, 6.1

        """

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
        """Test column desnity for aperture << lengthscale."""
        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        aper = 10 * u.km
        parent = 1e4 * u.km
        sigma = Haser(Q, v, parent).column_density(aper)
        ideal = Q / v / 2 / aper
        assert np.isclose(sigma.decompose().value, ideal.decompose().value)

    def test_column_density_large_aperture(self):
        """Test column desnity for aperture >> lengthscale."""
        Q = 1 / u.s
        v = 1 * u.km / u.s
        aper = 1000 * u.km
        parent = 10 * u.km
        sigma = Haser(Q, v, parent).column_density(aper)
        ideal = 4 * Q * parent / v / np.pi / 4 / aper**2
        assert np.isclose(sigma.decompose().value, ideal.decompose().value)

    def test_total_number_rho(self):
        """Reproduce Newburn and Johnson 1978.

        Species, N observed, Q/v (km**-1)
        CN, 6.4e26, 5.8e23
        C3, 8.3e28, 9.0e23
        C2, 7.8e27, 5.9e24

        """

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

    def test_total_number_circular_ap(self):
        """Reproduce Newburn and Johnson 1978.

        Species, N observed, Q/v (km**-1)
        CN, 6.4e26, 5.8e23
        C3, 8.3e28, 9.0e23
        C2, 7.8e27, 5.9e24

        """
        
        from ..core import CircularAperture
        
        Nobs = np.array([6.41756750e26, 8.63191842e+28, 7.81278300e27])
        parent = [1.4e4, 0, 1.0e4] * u.km
        daughter = [1.7e5, 4.6e4, 7.6e4] * u.km
        Q = [5.8e23, 9.0e23, 5.9e24] / u.s
        aper = CircularAperture(3300 * u.km)

        N = np.zeros(3)
        for i in range(3):
            coma = Haser(Q[i], 1 * u.km / u.s, parent[i], daughter[i])
            N[i] = coma.total_number(aper)

        assert np.allclose(N, Nobs)

    def test_total_number_rectangular_ap(self):
        """

        compare with:
        
        import astropy.units as u
        from sbpy.imageanalysis.utils import rarray
        from sbpy.activity import Haser

        r = rarray((5000, 3300)) * u.km
        coma = Haser(5.8e23 / u.s, 1 * u.km / u.s, 1.4e4 * u.km, 1.7e5 * u.km)
        sigma = coma.column_density(r)
        print(sigma.value.sum())
        --> 3.58626973105e+25
        
        This differs from the test value below by 1.3%.

        """
        
        from ..core import RectangularAperture
        
        parent = 1.4e4 * u.km
        daughter = 1.7e5 * u.km
        Q = 5.8e23 / u.s
        aper = RectangularAperture([5000, 3300] * u.km)

        coma = Haser(Q, 1 * u.km / u.s, parent, daughter)
        N = coma.total_number(aper)

        assert np.isclose(N, 3.6329898239171177e25)

    def test_total_number_gaussian_ap(self):
        """

        Compare with:
        
        import astropy.units as u
        from sbpy.imageanalysis.utils import rarray
        from sbpy.activity import Haser, GaussianAperture

        coma = Haser(5.8e23 / u.s, 1 * u.km / u.s, 1.4e4 * u.km)
        aper = GaussianAperture(1e4 * u.km)

        r = rarray((1000, 1000)) * 100 * u.km
        x = coma.column_density(r) * aper(r)
        print(x.value.sum() * 100**2)
        
        --> 5.146824269306973e+27

        which is within 0.5% of the test value below

        """
        
        from ..core import GaussianAperture
        
        parent = 1.4e4 * u.km
        Q = 5.8e23 / u.s
        aper = GaussianAperture(1e4 * u.km)

        coma = Haser(Q, 1 * u.km / u.s, parent)
        N = coma.total_number(aper)

        assert np.isclose(N, 5.17022685108891e+27)
