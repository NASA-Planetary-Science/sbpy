# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
try:
    import scipy
except ImportError:
    scipy = None
import astropy.units as u
import astropy.constants as const
from .. import core
from .. import (
    photo_lengthscale,
    photo_timescale,
    fluorescence_band_strength,
    VectorialModel,
    Haser,
)
from ....data import Phys


def test_photo_lengthscale():
    gamma = photo_lengthscale("OH", "CS93")
    assert gamma == 1.6e5 * u.km


def test_photo_lengthscale_error():
    with pytest.raises(ValueError):
        photo_lengthscale("asdf")

    with pytest.raises(ValueError):
        photo_lengthscale("OH", source="asdf")


def test_photo_timescale():
    tau = photo_timescale("CO2", "CE83")
    assert tau == 5.0e5 * u.s


def test_photo_timescale_error():
    with pytest.raises(ValueError):
        photo_timescale("asdf")

    with pytest.raises(ValueError):
        photo_timescale("OH", source="asdf")


@pytest.mark.parametrize(
    "band, test",
    (
        ("OH 0-0", 1.54e-15 * u.erg / u.s),
        ("OH 1-0", 1.79e-16 * u.erg / u.s),
        ("OH 1-1", 2.83e-16 * u.erg / u.s),
        ("OH 2-2", 1.46e-18 * u.erg / u.s),
        ("OH 0-1", 0.00356 * 1.54e-15 * u.erg / u.s),
        ("OH 0-2", 0.00021 * 1.54e-15 * u.erg / u.s),
        ("OH 1-2", 0.00610 * 2.83e-16 * u.erg / u.s),
        ("OH 2-0", 0.274 * 1.46e-18 * u.erg / u.s),
        ("OH 2-1", 1.921 * 1.46e-18 * u.erg / u.s),
    ),
)
def test_fluorescence_band_strength_OH_SA88(band, test):
    "Tests values for -1 km/s at 1 au"
    eph = {"rh": [1, 2] * u.au, "rdot": [-1, -1] * u.km / u.s}
    LN = fluorescence_band_strength(band, eph, "SA88").to(test.unit)
    assert np.allclose(LN.value, test.value / np.r_[1, 2] ** 2)


def test_fluorescence_band_strength_error():
    with pytest.raises(ValueError):
        fluorescence_band_strength("asdf")

    with pytest.raises(ValueError):
        fluorescence_band_strength("OH 0-0", source="asdf")


class TestHaser:
    def test_column_density_small_aperture(self):
        """Test column density for aperture << lengthscale.

        Should be within 1% of ideal value.

        """

        pytest.importorskip("scipy")

        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        rho = 1 * u.km
        parent = 1e4 * u.km
        N_avg = 2 * Haser(Q, v, parent).column_density(rho)
        ideal = Q / v / 2 / rho
        assert np.isclose(N_avg.decompose().value,
                          ideal.decompose().value, rtol=0.001)

    def test_column_density_small_angular_aperture(self):
        """Test column density for angular aperture << lengthscale.

        Regression test for PR#243.

        Should be within 1% of ideal value.

        """

        pytest.importorskip("scipy")

        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        rho = 0.001 * u.arcsec
        eph = dict(delta=1 * u.au)
        parent = 1e4 * u.km
        N_avg = 2 * Haser(Q, v, parent).column_density(rho, eph)
        rho_km = (rho * eph["delta"] * 725.24 *
                  u.km / u.arcsec / u.au).to("km")
        ideal = Q / v / 2 / rho_km
        assert np.isclose(N_avg.to_value("1/m2"),
                          ideal.to_value("1/m2"), rtol=0.001)

    def test_column_density(self):
        """
        Test column density for aperture = lengthscale.

        """

        pytest.importorskip("scipy")

        Q = 1e28 / u.s
        v = 1 * u.km / u.s
        rho = 1000 * u.km
        parent = 1000 * u.km
        coma = Haser(Q, v, parent)
        N_avg = coma.column_density(rho)
        integral = coma._integrate_volume_density(rho.to("m").value)[0]
        assert np.isclose(N_avg.decompose().value, integral)

    def test_total_number_large_aperture(self):
        """Test column density for aperture >> lengthscale."""

        pytest.importorskip("scipy")

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

        pytest.importorskip("scipy")

        Q = 1e25 / u.s
        v = 1 * u.km / u.s
        parent = 24000 * u.km
        daughter = 160000 * u.km
        ap = core.CircularAperture(1 * u.arcsec).as_length(1 * u.au)
        coma = Haser(Q, v, parent, daughter)
        N = coma.total_number(ap)
        assert np.isclose(N, 5.238964562688742e26)

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
        # tabAC = [
        #     [1.773, 40272, 4.738e30, 0.000000, 0.000000, 26.16, 0.000, 0.00],
        #     [1.053, 43451, 3.189e31, 1.558e31, 1.054e31, 26.98, 26.52, 26.5],
        #     [0.893, 39084, 3.147e31, 1.129e31, 2.393e31, 27.12, 26.54, 27.0]
        # ]

        pytest.importorskip("scipy")

        # Computed by sbpy.  0.893 C2 and C3 matches are not great,
        # the rest are OK:
        tab = [
            [1.773, 40272, 4.603e30, 0.000000, 0.000000, 26.16, 0.000, 0.00],
            [1.053, 43451, 3.229e31, 1.354e31, 9.527e30, 26.98, 26.52, 26.5],
            [0.893, 39084, 4.097e31, 1.273e31, 2.875e31, 27.12, 26.54, 27.0],
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

        pytest.importorskip("scipy")

        # Nobs = 2.314348613550494e+27
        parent = 1.4e4 * u.km
        Q = 5.8e23 / u.s
        v = 1 * u.km / u.s
        aper = core.CircularAperture(3300 * u.km)

        coma = Haser(Q, v, parent)
        N1 = coma.total_number(aper)
        N2 = coma._integrate_column_density(aper)[0]
        assert np.allclose(N1, N2)

    def test_circular_integration_1step(self):
        """Compare total number and integrated column density for circle.

        parent + daughter

        """

        pytest.importorskip("scipy")

        # Nobs = 6.41756750e26
        parent = 1.4e4 * u.km
        daughter = 1.7e5 * u.km
        Q = 5.8e23 / u.s
        v = 1 * u.km / u.s
        aper = core.CircularAperture(3300 * u.km)

        coma = Haser(Q, v, parent, daughter)
        N1 = coma.total_number(aper)
        N2 = coma._integrate_column_density(aper)[0]
        assert np.isclose(N1, N2)

    def test_total_number_annulus(self):
        """Test column density for annular aperture."""

        pytest.importorskip("scipy")

        Q = 1 / u.s
        v = 1 * u.km / u.s
        aper = core.AnnularAperture((1000, 2000) * u.km)
        parent = 10 * u.km
        N = Haser(Q, v, parent).total_number(aper)

        N1 = Haser(Q, v, parent).total_number(
            core.CircularAperture(aper.dim[0]))
        N2 = Haser(Q, v, parent).total_number(
            core.CircularAperture(aper.dim[1]))

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

        pytest.importorskip("scipy")

        parent = 1.4e4 * u.km
        daughter = 1.7e5 * u.km
        Q = 5.8e23 / u.s
        v = 1 * u.km / u.s
        aper = core.RectangularAperture([5000, 3300] * u.km)

        coma = Haser(Q, v, parent, daughter)
        N = coma.total_number(aper)

        assert np.isclose(N, 3.449607967230623e26)

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

        pytest.importorskip("scipy")

        parent = 1.4e4 * u.km
        Q = 5.8e23 / u.s
        aper = core.GaussianAperture(1e4 * u.km)

        coma = Haser(Q, 1 * u.km / u.s, parent)
        N = coma.total_number(aper)

        assert np.isclose(N, 5.146824269306973e27, rtol=0.005)


@pytest.mark.skipif(scipy is None, reason="requires scipy")
class TestVectorialModel:
    def test_small_vphoto(self):
        """
        The other test using water as parent and hydroxyl as fragment have a
        v_photo > v_outflow, but the model has a slightly different case for
        v_photo < v_outflow
        """

        # Across python, fortran, and rust models, we get very very close to
        # this number of fragments
        num_fragments_grid = 1.2162140e33

        base_q = 1.0e28 * 1 / u.s

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * u.s,
                "tau_d": 101730 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH, but v_photo is modified to be smaller than
        # v_outflow
        fragment = Phys.from_dict(
            {"tau_T": photo_timescale("OH") * 0.93,
             "v_photo": 0.5 * u.km / u.s}
        )

        coma = VectorialModel(
            base_q=base_q, parent=parent, fragment=fragment
        )

        assert np.isclose(
            coma.vmr.num_fragments_grid, num_fragments_grid, rtol=0.02
        )

    def test_time_dependent_function(self):
        """
        Test handing off a time dependence to the model with zero additional
        time-dependent production from q_t: results should match a model with
        steady production specified by base_q
        Also uses a model with print_progress=true to avoid the code coverage
        tests being polluted with trivial branches about the model
        conditionally printing
        """
        def q_t(t):
            # for all times, return zero additional production
            return t * 0

        base_q = 1.0e28 * 1 / u.s

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * u.s,
                "tau_d": 101730 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": photo_timescale("OH") * 0.93,
             "v_photo": 1.05 * u.km / u.s}
        )

        coma_steady = VectorialModel(
            base_q=base_q, parent=parent, fragment=fragment
        )

        coma_q_t = VectorialModel(
            base_q=base_q, parent=parent, fragment=fragment, q_t=q_t,
            print_progress=True
        )

        assert np.isclose(
            coma_steady.vmr.num_fragments_grid, coma_q_t.vmr.num_fragments_grid, rtol=0.02
        )

    def test_binned_production_one_element_list(self):
        """
        Initialize a comet with the fortran-version style of specifying time
        dependent production and make it essentially steady production, and
        test against a comet with the same steady production but initialized
        differently
        This also tests the model dealing with times in the past farther back
        than is specified in the variation list 'ts': it extends the oldest
        production value (here, 1.0e28) back infinitely into the past
        """
        # production from 1 day ago until the present,
        ts = [1] * u.day
        # is a steady value of 1e28
        qs = [1.0e28] / u.s

        base_q = 1.0e28 * 1 / u.s

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * u.s,
                "tau_d": 101730 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": photo_timescale("OH") * 0.93,
             "v_photo": 1.05 * u.km / u.s}
        )

        coma_binned = VectorialModel.binned_production(
            qs=qs, ts=ts, parent=parent, fragment=fragment
        )

        coma_steady = VectorialModel(
            base_q=base_q, parent=parent, fragment=fragment
        )

        assert np.isclose(
            coma_steady.vmr.num_fragments_grid, coma_binned.vmr.num_fragments_grid, rtol=0.02
        )

    def test_binned_production_multi_element_list(self):
        """
        Initialize a comet with the fortran-version style of specifying time
        dependent production and make it steady production, then test against a
        comet with the same steady production but initialized differently
        """
        # Specify that at multiple points in time,
        ts = [60, 50, 40, 30, 20, 10] * u.day
        # production is a steady value of 1e28
        qs = [1.0e28, 1.0e28, 1.0e28, 1.0e28, 1.0e28, 1.0e28] / u.s

        base_q = 1.0e28 * 1 / u.s

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * u.s,
                "tau_d": 101730 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": photo_timescale("OH") * 0.93,
             "v_photo": 1.05 * u.km / u.s}
        )

        coma_binned = VectorialModel.binned_production(
            qs=qs, ts=ts, parent=parent, fragment=fragment
        )

        coma_steady = VectorialModel(
            base_q=base_q, parent=parent, fragment=fragment
        )

        assert np.isclose(
            coma_steady.vmr.num_fragments_grid, coma_binned.vmr.num_fragments_grid, rtol=0.02
        )

    def test_grid_count(self):
        """
        Compute theoretical number of fragments vs. integrated value from
        grid. This is currently only a good estimate for steady production
        due to our method for determining the theoretical count of
        fragments
        """

        base_q = 1.0e28 * 1 / u.s

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * u.s,
                "tau_d": 101730 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": photo_timescale("OH") * 0.93,
             "v_photo": 1.05 * u.km / u.s}
        )

        coma = VectorialModel(base_q=base_q, parent=parent, fragment=fragment)

        assert np.isclose(
            coma.vmr.num_fragments_theory, coma.vmr.num_fragments_grid, rtol=0.02
        )

    def test_total_number_large_aperture(self):
        """
        Compare theoretical number of fragments vs. integration of column
        density over a large aperture
        """

        base_q = 1e28 * 1 / u.s

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * u.s,
                "tau_d": 101730 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": photo_timescale("OH") * 0.93,
             "v_photo": 1.05 * u.km / u.s}
        )

        coma = VectorialModel(base_q=base_q, parent=parent, fragment=fragment)

        ap = core.CircularAperture(coma.vmr.max_grid_radius)
        assert np.isclose(
            coma.vmr.num_fragments_theory, coma.total_number(ap), rtol=0.02
        )

    def test_model_symmetry(self):
        """
        The symmetry of the model allows the parent production to be
        treated as an overall scaling factor, and does not affect the
        features of the calculated densities.

        If we assume an arbitrary fixed count inside an aperture, we can
        use this to calculate what the production would need to be to
        produce this count after running the model with a 'dummy' value for
        the production.

        If we then run another model at this calculated production and use
        the same aperture, we should recover the number of counts assumed
        in the paragraph above.  If we do not, we have broken our model
        somehow and lost the symmetry during our calculations.
        """

        base_production = 1e28

        # # 100,000 x 100,000 km aperture
        ap = core.CircularAperture((1.0e6) * u.km)
        # assumed fragment count inside this aperture
        assumed_count = 1e32

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * u.s,
                "tau_d": 101730 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": photo_timescale("OH") * 0.93,
             "v_photo": 1.05 * u.km / u.s}
        )

        coma = VectorialModel(
            base_q=base_production * (1 / u.s), parent=parent, fragment=fragment
        )

        model_count = coma.total_number(ap)
        calculated_q = (assumed_count / model_count) * base_production

        # # Parent molecule is H2O
        # parent_check = Phys.from_dict({
        #     'tau_T': 86430 * u.s,
        #     'tau_d': 101730 * u.s,
        #     'v_outflow': 1 * u.km / u.s,
        #     'sigma': 3e-16 * u.cm ** 2
        # })
        # # Fragment molecule is OH
        # fragment_check = Phys.from_dict({
        #     'tau_T': photo_timescale('OH') * 0.93,
        #     'v_photo': 1.05 * u.km / u.s
        # })
        # coma_check = VectorialModel(base_q=calculated_q * (1 / u.s),
        #                             parent=parent_check,
        #                             fragment=fragment_check)
        coma_check = VectorialModel(
            base_q=calculated_q * (1 / u.s), parent=parent, fragment=fragment
        )
        count_check = coma_check.total_number(ap)

        assert np.isclose(count_check, assumed_count, rtol=0.001)

    @pytest.mark.parametrize(
        "rh,delta,flux,g,Q",
        (
            [1.2912, 0.7410, 337e-14, 2.33e-4, 1.451e28],
            [1.2949, 0.7651, 280e-14, 2.60e-4, 1.228e28],
            [1.3089, 0.8083, 480e-14, 3.36e-4, 1.967e28],
            [1.3200, 0.8353, 522e-14, 3.73e-4, 2.025e28],
            [1.3366, 0.8720, 560e-14, 4.03e-4, 2.035e28],
        ),
    )
    def test_festou92(self, rh, delta, flux, g, Q):
        """Compare to Festou et al. 1992 production rates of comet 6P/d'Arrest.

        Festou et al. 1992. The Gas Production Rate of Periodic Comet d'Arrest.
        Asteroids, Comets, Meteors 1991, 177.

        https://ui.adsabs.harvard.edu/abs/1992acm..proc..177F/abstract

        IUE observations of OH(0-0) band.

        The table appears to have a typo in the Q column for the 1.32 au
        observation: using 2.025 instead of 3.025.

        """

        # add units
        rh = rh * u.au
        delta = delta * u.au
        flux = flux * u.erg / u.s / u.cm**2
        g = g / u.s
        Q = Q / u.s

        # OH (0-0) luminosity per molecule
        L_N = g / (rh / u.au) ** 2 * const.h * const.c / (3086 * u.AA)

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 65000 * (rh / u.au) ** 2 * u.s,
                "tau_d": 72500 * (rh / u.au) ** 2 * u.s,
                "v_outflow": 0.85 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": 160000 * (rh / u.au) ** 2 * u.s,
             "v_photo": 1.05 * u.km / u.s}
        )

        # https://pds.nasa.gov/ds-view/pds/viewInstrumentProfile.jsp?INSTRUMENT_ID=LWP&INSTRUMENT_HOST_ID=IUE
        # Large-Aperture Length(arcsec)   22.51+/-0.40
        # Large-Aperture Width(arcsec)     9.91+/-0.17
        #
        # 10x20 quoted by Festou et al.
        # effective circle is 8.0 radius
        # half geometric mean is 7.07
        #
        # However, vectorial fortran code has 9.1x15.3 hard coded into it.
        lwp = core.RectangularAperture((9.1, 15.3) * u.arcsec)

        Q0 = 1e28 / u.s
        coma = VectorialModel(base_q=Q0, parent=parent, fragment=fragment)
        N0 = coma.total_number(lwp, eph=delta)

        Q_model = (Q0 * flux / (L_N * N0) * 4 * np.pi * delta**2).to(Q.unit)

        # absolute tolerance: Table 2 has 3 significant figures
        #
        # relative tolerance: actual agreement is 14%, future updates should
        # improve this or else explain the difference
        atol = 1.01 * 10 ** (np.floor(np.log10(Q0.value)) - 2) * Q.unit
        assert u.allclose(Q, Q_model, atol=atol, rtol=0.14)

    @pytest.mark.parametrize(
        "rh,delta,N,Q",
        (
            [1.8662, 0.9683, 0.2424e32, 1.048e29],
            [0.8855, 0.9906, 3.819e32, 5.548e29],
            [0.9787, 0.8337, 1.63e32, 3.69e29],
            [1.0467, 0.7219, 0.8703e32, 2.813e29],
            [1.9059, 1.4031, 1.07e32, 2.76e29],
            [2.0715, 1.7930, 1.01e32, 1.88e29],
        ),
    )
    def test_combi93(self, rh, delta, N, Q):
        """Compare to results of Combi et al. 1993.

        Combi et al. 1993 compared a Monte Carlo approach to the Vectorial
        model for OH.  They find best agreement between the two models near 1
        au, likely due to the assumption that the water outflow speed is
        constant in the VM runs, but the MC model only has 1 km/s speeds near 1
        au.

        """

        # assign units
        rh = rh * u.au
        delta = delta * u.au
        Q = Q / u.s  # Vectorial model run by Roettger
        aper = core.RectangularAperture((10, 15) * u.arcsec)

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 8.2e4 * 0.88 * (rh / u.au) ** 2 * u.s,
                "tau_d": 8.2e4 * (rh / u.au) ** 2 * u.s,
                "v_outflow": 1 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )

        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {"tau_T": 2.0e5 * (rh / u.au) ** 2 * u.s,
             "v_photo": 1.05 * u.km / u.s}
        )

        Q0 = 2e29 / u.s
        coma = VectorialModel(base_q=Q0, parent=parent, fragment=fragment)
        N0 = coma.total_number(aper, eph=delta)

        Q_model = (Q0 * N / N0).to(Q.unit)
        assert u.allclose(Q, Q_model, rtol=0.13)

    def test_vm_fortran(self):
        """Compare to results from vm.f.

        vm.f shared with MSK by Joel Parker.  Comments:

            WRITTEN BY M. C. FESTOU
            (Minor modifications were made by MF on 8 Nov. 1988 on the Tempe
            version) Some unused variables removed on 6 June 1995 in Baltimore
            when adapting the code to run on a unix machine. A little bit of
            additional trimming on 27 June 1995. New printing formats.  Few
            printing format changes made on 27 Sept. 96.  Formula on the
            collision radius added on 1 April 1997.  Modification in sub GAUSS
            (NP and coeff.) and in sub APP on 22-24/4/1997.

        Input fparam.dat:

            Wirtanen (HST-FOS slit)
            2.46900010       1.51400006
                    1
            9.99999988E+26   100.000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.00000000       0.00000000
            0.541000009
            86430.0000
            101730.000
            99.0000000
            OH
            4.15999995E-04
            1.04999995
            129000.000
            95.0000000
            1.28999996       3.66000009

        Results: see wm_fortran_test_output.txt in data directory.

        """

        eph = {"rh": 2.469 * u.au, "delta": 1.514 * u.au}

        # Parent molecule is H2O
        parent = Phys.from_dict(
            {
                "tau_T": 86430 * (eph["rh"] / u.au) ** 2 * u.s,
                "tau_d": 101730 * (eph["rh"] / u.au) ** 2 * u.s,
                "v_outflow": 0.514 * u.km / u.s,
                "sigma": 3e-16 * u.cm**2,
            }
        )
        # Fragment molecule is OH
        fragment = Phys.from_dict(
            {
                "tau_T": 129000 * (eph["rh"] / u.au) ** 2 * u.s,
                "v_photo": 1.05 * u.km / u.s,
            }
        )

        # test values are copy-pasted from wm.f output
        collision_sphere_radius = 0.14e07 * u.cm
        collision_sphere_radius_atol = 0.011e7 * u.cm
        N_fragments_theory = 0.668e33
        # N_fragments_theory_atol = 0.0011e33
        N_fragments_theory_rtol = 0.001
        N_fragments = 0.657e33
        # N_fragments_atol = 0.0011e33
        N_fragments_rtol = 0.003

        fragment_volume_density = """
0.97E+03  0.43E+03    0.68E+04  0.59E+02    0.13E+05  0.31E+02    0.18E+05  0.21E+02
0.24E+05  0.15E+02    0.31E+05  0.11E+02    0.43E+05  0.79E+01    0.54E+05  0.59E+01
0.66E+05  0.46E+01    0.78E+05  0.38E+01    0.90E+05  0.31E+01    0.11E+06  0.24E+01
0.13E+06  0.19E+01    0.14E+06  0.16E+01    0.16E+06  0.13E+01    0.18E+06  0.11E+01
0.21E+06  0.87E+00    0.24E+06  0.69E+00    0.27E+06  0.56E+00    0.30E+06  0.46E+00
0.33E+06  0.38E+00    0.37E+06  0.29E+00    0.42E+06  0.23E+00    0.47E+06  0.18E+00
0.51E+06  0.14E+00    0.56E+06  0.11E+00    0.63E+06  0.85E-01    0.70E+06  0.65E-01
0.77E+06  0.50E-01    0.84E+06  0.39E-01    0.92E+06  0.30E-01    0.10E+07  0.21E-01
0.11E+07  0.16E-01    0.12E+07  0.12E-01    0.13E+07  0.86E-02    0.14E+07  0.64E-02
0.16E+07  0.44E-02    0.17E+07  0.31E-02    0.19E+07  0.23E-02    0.20E+07  0.17E-02
0.22E+07  0.12E-02    0.24E+07  0.78E-03    0.27E+07  0.51E-03    0.29E+07  0.34E-03
0.31E+07  0.23E-03    0.34E+07  0.15E-03    0.37E+07  0.90E-04    0.41E+07  0.53E-04
0.44E+07  0.32E-04    0.48E+07  0.20E-04"""
        fragment_column_density = """
0.970E+03   4.77E+11    1.088E+03   4.69E+11    1.221E+03   4.61E+11    1.370E+03   4.53E+11
1.537E+03   4.44E+11    1.724E+03   4.34E+11    1.935E+03   4.22E+11    2.171E+03   4.15E+11
2.436E+03   4.00E+11    2.733E+03   3.85E+11    3.066E+03   3.74E+11    3.441E+03   3.67E+11
3.860E+03   3.57E+11    4.332E+03   3.47E+11    4.860E+03   3.36E+11    5.453E+03   3.27E+11
6.118E+03   3.18E+11    6.865E+03   3.08E+11    7.703E+03   2.99E+11    8.643E+03   2.89E+11
9.697E+03   2.79E+11    1.088E+04   2.70E+11    1.221E+04   2.61E+11    1.370E+04   2.51E+11
1.537E+04   2.42E+11    1.724E+04   2.32E+11    1.935E+04   2.23E+11    2.171E+04   2.13E+11
2.436E+04   2.04E+11    2.733E+04   1.95E+11    3.066E+04   1.86E+11    3.441E+04   1.76E+11
3.860E+04   1.67E+11    4.332E+04   1.58E+11    4.860E+04   1.49E+11    5.453E+04   1.41E+11
6.118E+04   1.32E+11    6.865E+04   1.23E+11    7.703E+04   1.15E+11    8.643E+04   1.07E+11
9.697E+04   9.87E+10    1.088E+05   9.09E+10    1.221E+05   8.33E+10    1.370E+05   7.60E+10
1.537E+05   6.89E+10    1.724E+05   6.21E+10    1.935E+05   5.56E+10    2.171E+05   4.95E+10
2.436E+05   4.37E+10    2.733E+05   3.83E+10    3.066E+05   3.32E+10    3.441E+05   2.86E+10
3.860E+05   2.43E+10    4.332E+05   2.05E+10    4.860E+05   1.71E+10    5.453E+05   1.40E+10
6.118E+05   1.14E+10    6.865E+05   9.11E+09    7.703E+05   7.18E+09    8.643E+05   5.57E+09
9.697E+05   4.25E+09    1.088E+06   3.19E+09    1.221E+06   2.34E+09    1.370E+06   1.68E+09
1.537E+06   1.19E+09    1.724E+06   8.29E+08    1.935E+06   5.67E+08    2.171E+06   3.78E+08
2.436E+06   2.45E+08    2.733E+06   1.54E+08    3.066E+06   9.34E+07    3.441E+06   5.39E+07
"""

        # convert strings to arrays
        x = np.fromstring(fragment_volume_density, sep=" ")
        n0_rho = x[::2] * u.km
        n0 = x[1::2] / u.cm**3
        # absolute tolerance: 2 significant figures
        n0_atol = 1.1 * 10 ** (np.floor(np.log10(n0.value)) - 1) * n0.unit
        n0_atol_revised = n0_atol * 7

        x = np.fromstring(fragment_column_density, sep=" ")
        sigma0_rho = x[::2] * u.km
        sigma0 = x[1::2] / u.cm**2
        # absolute tolerance: 3 significant figures
        sigma0_atol = 1.1 * \
            10 ** (np.floor(np.log10(sigma0.value)) - 2) * sigma0.unit
        sigma0_atol_revised = sigma0_atol * 40

        # evaluate the model
        Q0 = 1e27 / u.s
        coma = VectorialModel(base_q=Q0, parent=parent, fragment=fragment)
        n = [coma.volume_density(r) for r in n0_rho]
        sigma = [coma.column_density(r) for r in sigma0_rho]

        # compare results
        assert u.isclose(
            coma.vmr.collision_sphere_radius,
            collision_sphere_radius,
            atol=collision_sphere_radius_atol,
        )

        assert np.isclose(
            coma.vmr.num_fragments_theory,
            N_fragments_theory,
            rtol=N_fragments_theory_rtol,
        )

        assert np.isclose(
            coma.vmr.num_fragments_grid, N_fragments, rtol=N_fragments_rtol
        )

        assert u.allclose(n, n0, atol=n0_atol_revised)

        assert u.allclose(sigma, sigma0, atol=sigma0_atol_revised)
