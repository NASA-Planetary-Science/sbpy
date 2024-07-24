# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from ... import units as sbu
from ..core import *


class TestCircularAperture:
    def test_init_unit(self):
        with pytest.raises(u.UnitTypeError):
            CircularAperture(1 * u.s)

    def test_str(self):
        assert str(CircularAperture(1 * u.arcsec)
                   ) == 'Circular aperture, radius 1.0 arcsec'

    def test_coma_equivalent_radius(self):
        r = 1 * u.arcsec
        aper = CircularAperture(r)
        assert aper.coma_equivalent_radius() == 1 * u.arcsec

    def test_as_length(self):
        r = 1 * u.arcsec
        aper = CircularAperture(r)
        eph = {'delta': 1 * u.au}
        length = r.to('km', sbu.projected_size(eph))
        assert np.isclose(aper.as_length(eph).dim.value, length.value)

    def test_as_angle(self):
        r = 100 * u.km
        aper = CircularAperture(r)
        eph = {'delta': 1 * u.au}
        angle = r.to('arcsec', sbu.projected_size(eph))
        assert np.isclose(aper.as_angle(eph).dim.value, angle.value)

    def test_from_coma_equivalent(self):
        # test initialization from another aperture
        shape = [1, 2] * u.arcsec
        an_aper = AnnularAperture(shape)
        circ_aper = CircularAperture.from_coma_equivalent(an_aper)
        assert circ_aper.dim == 1 * u.arcsec

        # test initialization from a radius
        circ_aper = CircularAperture.from_coma_equivalent(1 * u.arcsec)
        assert circ_aper.dim == 1 * u.arcsec


class TestAnnularAperture:
    def test_str(self):
        assert str(AnnularAperture([1, 2] * u.arcsec)
                   ) == 'Annular aperture, radii 1.0–2.0 arcsec'

    def test_shape(self):
        assert np.allclose(
            AnnularAperture([1, 2] * u.arcsec).shape.value,
            (1, 2))

    def test_coma_equivalent_radius(self):
        shape = [1, 2] * u.arcsec
        aper = AnnularAperture(shape)
        assert aper.coma_equivalent_radius() == 1 * u.arcsec

    def test_shape_error(self):
        with pytest.raises(ValueError):
            AnnularAperture([1, 2, 3] * u.km)


class TestRectangularAperture:
    def test_str(self):
        assert str(RectangularAperture(
            [1, 2] * u.arcsec)) == 'Rectangular aperture, dimensions 1.0×2.0 arcsec'

    def test_coma_equivalent_radius(self):
        shape = (0.8, 2) * u.arcsec
        aper = RectangularAperture(shape)
        r = aper.coma_equivalent_radius()
        assert np.isclose(r.value, 0.66776816346357259)

    def test_shape_error(self):
        with pytest.raises(ValueError):
            RectangularAperture([1, 2, 3] * u.km)


class TestGaussianAperture:
    def test_init_fwhm(self):
        aper = GaussianAperture(fwhm=1 * u.arcsec)
        assert np.isclose(aper.sigma.value, 1 / 2.3548200450309493)
        assert np.isclose(aper.fwhm.value, 1)

    def test_init_error(self):
        with pytest.raises(ValueError):
            GaussianAperture()

    def test_str(self):
        assert str(GaussianAperture(1 * u.arcsec)
                   ) == 'Gaussian aperture, 1-σ width 1.0 arcsec'

    def test_coma_equivalent_radius(self):
        sigma = 1 * u.arcsec
        aper = GaussianAperture(sigma)
        assert aper.coma_equivalent_radius() == np.sqrt(np.pi / 2) * u.arcsec

    def test_call_angle(self):
        aper = GaussianAperture(1 * u.arcsec)
        assert np.isclose(aper(1 * u.arcsec), 0.6065306597126334)

    def test_call_mixed(self):
        aper = GaussianAperture(1 * u.arcsec)
        assert np.isclose(aper(725.24 * u.km, {'delta': 1 * u.au}),
                          0.60653, rtol=0.0001)

        aper = GaussianAperture(725.24 * u.km)
        assert np.isclose(aper(1 * u.arcsec, {'delta': 1 * u.au}),
                          0.60653, rtol=0.0001)
