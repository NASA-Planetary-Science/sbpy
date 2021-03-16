# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from .. import core
from .. import *
from .... import exceptions as sbe


def test_photo_lengthscale():
    gamma = photo_lengthscale('OH', 'CS93')
    assert gamma == 1.6e5 * u.km


def test_photo_lengthscale_error():
    with pytest.raises(ValueError):
        photo_lengthscale('asdf')

    with pytest.raises(ValueError):
        photo_lengthscale('OH', source='asdf')


def test_photo_timescale():
    tau = photo_timescale('CO2', 'CE83')
    assert tau == 5.0e5 * u.s


def test_photo_timescale_error():
    with pytest.raises(ValueError):
        photo_timescale('asdf')

    with pytest.raises(ValueError):
        photo_timescale('OH', source='asdf')


@pytest.mark.parametrize('band, test', (
    ('OH 0-0', 1.54e-15 * u.erg / u.s),
    ('OH 1-0', 1.79e-16 * u.erg / u.s),
    ('OH 1-1', 2.83e-16 * u.erg / u.s),
    ('OH 2-2', 1.46e-18 * u.erg / u.s),
    ('OH 0-1', 0.00356 * 1.54e-15 * u.erg / u.s),
    ('OH 0-2', 0.00021 * 1.54e-15 * u.erg / u.s),
    ('OH 1-2', 0.00610 * 2.83e-16 * u.erg / u.s),
    ('OH 2-0', 0.274 * 1.46e-18 * u.erg / u.s),
    ('OH 2-1', 1.921 * 1.46e-18 * u.erg / u.s),
))
def test_fluorescence_band_strength_OH_SA88(band, test):
    "Tests values for -1 km/s at 1 au"
    eph = {
        'rh': [1, 2] * u.au,
        'rdot': [-1, -1] * u.km / u.s
    }
    LN = fluorescence_band_strength(band, eph, 'SA88').to(test.unit)
    assert np.allclose(LN.value, test.value / np.r_[1, 2]**2)


def test_fluorescence_band_strength_error():
    with pytest.raises(ValueError):
        fluorescence_band_strength('asdf')

    with pytest.raises(ValueError):
        fluorescence_band_strength('OH 0-0', source='asdf')


def test_gascoma_scipy_error(monkeypatch):
    monkeypatch.setattr(core, 'scipy', None)
    test = Haser(1 / u.s, 1 * u.km / u.s, 1e6 * u.km)
    with pytest.raises(sbe.RequiredPackageUnavailable):
        test._integrate_volume_density(1e5)

    with pytest.raises(sbe.RequiredPackageUnavailable):
        aper = core.CircularAperture(1000 * u.km)
        test._integrate_column_density(aper)


class TestHaser:
    def test_volume_density(self):
        """Test a set of dummy values."""
        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        parent = 1e4 * u.km
        daughter = 1e5 * u.km
        r = np.logspace(1, 7) * u.km
        n = Haser(Q, v, parent, daughter).volume_density(r)
        rel = (daughter / (parent - daughter)
               * (np.exp(-r / parent) - np.exp(-r / daughter)))
        # test radial profile
        assert np.allclose((n / n[0]).value,
                           (rel / rel[0] * (r[0] / r)**2).value)

        # test parent-only coma near nucleus against that expected for
        # a long-lived species; will be close, but not exact
        n = Haser(Q, v, parent).volume_density(10 * u.km)
        assert np.isclose(
            n.decompose().value,
            (Q / v / 4 / np.pi / (10 * u.km)**2).decompose().value,
            rtol=0.001)

    def test_column_density_small_aperture(self):
        """Test column density for aperture << lengthscale.

        Should be within 1% of ideal value.

        """
        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        rho = 1 * u.km
        parent = 1e4 * u.km
        N_avg = 2 * Haser(Q, v, parent).column_density(rho)
        ideal = Q / v / 2 / rho
        assert np.isclose(N_avg.decompose().value, ideal.decompose().value,
                          rtol=0.001)

    def test_column_density_small_angular_aperture(self):
        """Test column density for angular aperture << lengthscale.

        Regression test for PR#243.

        Should be within 1% of ideal value.

        """
        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        rho = 0.001 * u.arcsec
        eph = dict(delta=1 * u.au)
        parent = 1e4 * u.km
        N_avg = 2 * Haser(Q, v, parent).column_density(rho, eph)
        rho_km = (rho * eph['delta'] * 725.24 * u.km / u.arcsec / u.au).to('km')
        ideal = Q / v / 2 / rho_km
        assert np.isclose(N_avg.to_value('1/m2'), ideal.to_value('1/m2'), rtol=0.001)

    def test_column_density(self):
        """
        Test column density for aperture = lengthscale.

        """
        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        rho = 1000 * u.km
        parent = 1000 * u.km
        coma = Haser(Q, v, parent)
        N_avg = coma.column_density(rho)
        integral = coma._integrate_volume_density(rho.to('m').value)[0]
        assert np.isclose(N_avg.decompose().value, integral)

    def test_total_number_large_aperture(self):
        """Test column density for aperture >> lengthscale."""
        Q = 1 / u.s
        v = 1 * u.km / u.s
        rho = 1000 * u.km
        parent = 10 * u.km
        N = Haser(Q, v, parent).total_number(rho)
        ideal = Q * parent / v
        assert np.isclose(N, ideal.decompose().value)

    def test_total_number_large_aperture(self):
        """Test column density for aperture >> lengthscale."""
        Q = 1 / u.s
        v = 1 * u.km / u.s
        rho = 1000 * u.km
        parent = 10 * u.km
        N = Haser(Q, v, parent).total_number(rho)
        ideal = Q * parent / v
        assert np.isclose(N, ideal.decompose().value)

    def test_total_number_circular_aperture_angular(self):
        """Regression test for issue #239.

        https://github.com/NASA-Planetary-Science/sbpy/issues/239

        Code initially from Quanzhi Ye.

        """

        Q = 1e25 / u.s
        v = 1 * u.km / u.s
        parent = 24000 * u.km
        daughter = 160000 * u.km
        ap = core.CircularAperture(1 * u.arcsec).as_length(1 * u.au)
        coma = Haser(Q, v, parent, daughter)
        N = coma.total_number(ap)
        assert np.isclose(N, 5.238964562688742e+26)

    def test_total_number_rho_AC75(self):
        """Reproduce A'Hearn and Cowan 1975

        Assumed 1 km/s.

        N(C2) = 10**(12.9300 + log10(L) + 2 * log10(rh))
        N(CN) = 10**(12.3718 + log10(L) + 2 * log10(rh))
        N(C3) = 10**(13.6    + log10(L) + 2 * log10(rh))

        Species, parent, daughter  (km)
        C2,      1.0e4, 6.61e4
        CN,      1.3e4, 1.48e5
        C3,          0, 4.0e4

        # au, km, erg/s/cm2
        tab = ascii.read('''
        rh,    delta, log rho,  logC2,  logC3,  logCN,   QC2,   QCN,  QC3
        1.773, 2.465, 4.605,  -10.981,      0,      0, 26.16,     0,    0
        1.053, 1.535, 4.399,   -9.624,      0, -9.555, 26.98, 26.52, 26.5
        1.053, 1.535, 4.638,   -9.289, -10.44, -9.042, 26.98, 26.52, 26.5
        0.893, 1.383, 4.592,   -9.061,  -9.85, -8.948, 27.12, 26.54, 27.0
        ''')

        NC2 = 10**(12.9300 + log10(4 * pi * (tab['delta'] * 1.49e13)**2 *
                   10**tab['logC2']) + 2 * log10(tab['rh']))
        NCN = 10**(12.3718 + log10(4 * pi * (tab['delta'] * 1.49e13)**2 *
                   10**tab['logCN']) + 2 * log10(tab['rh']))
        NC3 = 10**(13.6 + log10(4 * pi * (tab['delta'] * 1.49e13)**2 *
                   10**tab['logC3']) + 2 * log10(tab['rh']))
        NC2.name = 'NC2'
        NCN.name = 'NCN'
        NC3.name = 'NC3'

        tab2 = Table((tab['rh'], 10**tab['log rho'], NC2, NCN, NC3,
                      tab['QC2'], tab['QCN'], tab['QC3']))

        """

        # A'Hearn and Cowan 1975:
        # rh      rho     NC2       NCN        NC3      QC2     QCN     QC3
        tabAC = [
            [1.773, 40272, 4.738e30, 0.000000, 0.000000, 26.16, 0.000, 0.00],
            [1.053, 43451, 3.189e31, 1.558e31, 1.054e31, 26.98, 26.52, 26.5],
            [0.893, 39084, 3.147e31, 1.129e31, 2.393e31, 27.12, 26.54, 27.0]
        ]

        # Computed by sbpy.  0.893 C2 and C3 matches are not great,
        # the rest are OK:
        tab = [
            [1.773, 40272, 4.603e30, 0.000000, 0.000000, 26.16, 0.000, 0.00],
            [1.053, 43451, 3.229e31, 1.354e31, 9.527e30, 26.98, 26.52, 26.5],
            [0.893, 39084, 4.097e31, 1.273e31, 2.875e31, 27.12, 26.54, 27.0]
        ]

        for rh, rho, NC2, NCN, NC3, QC2, QCN, QC3 in tab:
            if NC2 > 0:
                parent = 1.0e4 * u.km
                daughter = 6.61e4 * u.km
                Q = 10**QC2 / u.s
                coma = Haser(Q, 1 * u.km / u.s, parent, daughter)
                N = coma.total_number(rho * u.km)
                assert np.isclose(NC2, N, rtol=0.01)

            if NCN > 0:
                parent = 1.3e4 * u.km
                daughter = 1.48e5 * u.km
                Q = 10**QCN / u.s
                coma = Haser(Q, 1 * u.km / u.s, parent, daughter)
                N = coma.total_number(rho * u.km)
                assert np.isclose(NCN, N, rtol=0.01)

            if NC3 > 0:
                parent = 0 * u.km
                daughter = 4.0e4 * u.km
                Q = 10**QC3 / u.s
                coma = Haser(Q, 1 * u.km / u.s, parent, daughter)
                N = coma.total_number(rho * u.km)
                assert np.isclose(NC3, N, rtol=0.01)

    def test_circular_integration_0step(self):
        """Compare total number and integrated column density for circle.

        parent only

        """

        from ..core import CircularAperture

        Nobs = 2.314348613550494e+27
        parent = 1.4e4 * u.km
        Q = 5.8e23 / u.s
        v = 1 * u.km / u.s
        aper = CircularAperture(3300 * u.km)

        coma = Haser(Q, v, parent)
        N1 = coma.total_number(aper)
        N2 = coma._integrate_column_density(aper)[0]
        assert np.allclose(N1, N2)

    def test_circular_integration_1step(self):
        """Compare total number and integrated column density for circle.

        parent + daughter

        """

        from ..core import CircularAperture

        Nobs = 6.41756750e26
        parent = 1.4e4 * u.km
        daughter = 1.7e5 * u.km
        Q = 5.8e23 / u.s
        v = 1 * u.km / u.s
        aper = CircularAperture(3300 * u.km)

        coma = Haser(Q, v, parent, daughter)
        N1 = coma.total_number(aper)
        N2 = coma._integrate_column_density(aper)[0]
        assert np.isclose(N1, N2)

    def test_total_number_annulus(self):
        """Test column density for annular aperture."""

        from ..core import CircularAperture, AnnularAperture

        Q = 1 / u.s
        v = 1 * u.km / u.s
        aper = AnnularAperture((1000, 2000) * u.km)
        parent = 10 * u.km
        N = Haser(Q, v, parent).total_number(aper)

        N1 = Haser(Q, v, parent).total_number(CircularAperture(aper.dim[0]))
        N2 = Haser(Q, v, parent).total_number(CircularAperture(aper.dim[1]))

        assert np.allclose(N, N2 - N1)

    def test_total_number_rectangular_ap(self):
        """

        compare with:

        import astropy.units as u
        from sbpy.imageanalysis.utils import rarray
        from sbpy.activity import Haser

        r = rarray((5000, 3300), subsample=10) * u.km
        parent = 1.4e4 * u.km
        daughter = 1.7e5 * u.km
        Q = 5.8e23 / u.s
        v = 1 * u.km / u.s
        coma = Haser(Q, v, parent, daughter)
        sigma = coma.column_density(r)
        print((sigma * 1 * u.km**2).decompose().sum())
        --> <Quantity 3.449607967230623e+26>

        """

        from ..core import RectangularAperture

        parent = 1.4e4 * u.km
        daughter = 1.7e5 * u.km
        Q = 5.8e23 / u.s
        v = 1 * u.km / u.s
        aper = RectangularAperture([5000, 3300] * u.km)

        coma = Haser(Q, v, parent, daughter)
        N = coma.total_number(aper)

        assert np.isclose(N, 3.449607967230623e+26)

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

        assert np.isclose(N, 5.146824269306973e+27, rtol=0.005)

    def test_missing_scipy(self, monkeypatch):
        monkeypatch.setattr(core, 'scipy', None)
        test = Haser(1 / u.s, 1 * u.km / u.s, 1e6 * u.km)
        with pytest.raises(sbe.RequiredPackageUnavailable):
            test._iK0(1)
        with pytest.raises(sbe.RequiredPackageUnavailable):
            test._K1(1)

class TestVectorialModel:
    def test_fragment_count(self):
        """ Compute theoretical number of fragments vs. integrated value from grid """
        """ Dependent on the default size of the grid - will fail if that is too small """
        vmInput = {}

        # 45 days ago until now, production has been 1e28 molecules/sec
        vmInput['TimeAtProductions'] = [45] * u.day
        vmInput['ProductionRates'] = [1.e28]

        # Parent molecule is H2O
        vmInput['Parent'] = {}
        # A few values used in Festou's original fortran
        vmInput['Parent']['TotalLifetime'] = 86430 * u.s
        vmInput['Parent']['DissociativeLifetime'] = 101730 * u.s
        vmInput['Parent']['Velocity'] = 0.85 * (u.km/u.s)

        # Fragment molecule is OH
        vmInput['Fragment'] = {}
        vmInput['Fragment']['Velocity'] = 1.05 * u.km/u.s
        vmInput['Fragment']['TotalLifetime'] = 129000 * u.s

        coma = VectorialModel(0*(1/u.s), 0 * u.m/u.s, vmInput)

        fragTheory = coma.vModel['NumFragmentsTheory']
        fragGrid = coma.vModel['NumFragmentsFromGrid']
        assert np.isclose(fragTheory, fragGrid, rtol=0.02)

    def test_total_number_large_aperture(self):
        """ Compare theoretical number of fragments vs. integration of column density over a large aperture """
        """ Dependent on the default size of the grid - will fail if that is too small """
        vmInput = {}

        # 45 days ago until now, production has been 1e28 molecules/sec
        vmInput['TimeAtProductions'] = [45] * u.day
        vmInput['ProductionRates'] = [1.e28]

        # Parent molecule is H2O
        vmInput['Parent'] = {}
        # A few values used in Festou's original fortran
        vmInput['Parent']['TotalLifetime'] = 86430 * u.s
        vmInput['Parent']['DissociativeLifetime'] = 101730 * u.s
        vmInput['Parent']['Velocity'] = 0.85 * (u.km/u.s)

        # Fragment molecule is OH
        vmInput['Fragment'] = {}
        vmInput['Fragment']['Velocity'] = 1.05 * u.km/u.s
        vmInput['Fragment']['TotalLifetime'] = 129000 * u.s

        coma = VectorialModel(0*(1/u.s), 0 * u.m/u.s, vmInput)

        fragTheory = coma.vModel['NumFragmentsTheory']
        ap = core.CircularAperture((coma.vModel['MaxRadiusOfGrid'].value) * u.m)
        assert np.isclose(fragTheory, coma.total_number(ap), rtol=0.02)
