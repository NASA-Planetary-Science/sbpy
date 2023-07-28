# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import astropy.units as u
from .. import Haser, photo_timescale
from ....data import Ephem, Phys
from .. import (LTE, einstein_coeff, intensity_conversion, beta_factor,
                total_number, from_Haser)

'''
import os


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


class MockResponseSpec(object):

    def __init__(self, filename):
        self.filename = data_path(filename)

    @property
    def text(self):
        with open(self.filename) as f:
            return f.read()
'''


class TestProductionRate:

    def test_simple_prodrate(self):

        temp_estimate = 47. * u.K
        vgas = 0.8 * u.km / u.s
        aper = 30 * u.m
        b = 1.13
        mol_tag = 27001
        transition_freq = (265.886434 * u.GHz).to('MHz')
        mol_tag = 'co'
        t_freq = (157.178987 * u.GHz).to('MHz')
        lgint300 = 7.4336139*10**(-5) * u.MHz * u.nm * u.nm
        part300 = 9473.271845042498
        partition = 308.5479509730844
        gu = 11.0
        energy_J = 6.61822698*10**(-22) * u.J
        elo_J = 5.57674801*10**(-22) * u.J
        df = 3
        integrated_flux = 0.672 * u.K * u.km / u.s

        quantities = [t_freq, temp_estimate, lgint300, part300, partition, gu, energy_J,
                      elo_J, df, mol_tag]

        names = ['t_freq', 'temp', 'lgint300', 'partfn300', 'partfn',
                 'dgup', 'eup_J', 'elo_J', 'degfreedom', 'mol_tag']

        mol_data = Phys.from_dict(dict(zip(names, quantities)))

        intl = intensity_conversion(mol_data)

        mol_data.apply([intl.value] * intl.unit, name='intl')

        au = einstein_coeff(mol_data)

        mol_data.apply([au.value] * au.unit, name='eincoeff')

        r = 1.06077567 * u.AU

        delta = 0.14633757 * u.AU

        quanteph = [r, delta]

        nameseph = ['r', 'delta']

        ephemobj = Ephem.from_dict(dict(zip(nameseph, quanteph)))

        lte = LTE()

        q = lte.from_Drahus(integrated_flux, mol_data, ephemobj, vgas, aper, b=b)

        q = np.log10(q.value)

        assert np.isclose(q, 27.3167, rtol=0.01)

    def test_Haser_prodrate(self):
        """Test a set of dummy values."""

        pytest.importorskip("scipy")

        temp_estimate = 47. * u.K
        vgas = 0.8 * u.km / u.s
        aper = 30 * u.m
        b = 1.13
        mol_tag = 'CO'
        t_freq = (157.178987 * u.GHz).to('MHz')
        lgint300 = 7.4336139*10**(-5) * u.MHz * u.nm * u.nm
        part300 = 9473.271845042498
        partition = 308.5479509730844
        gu = 11.0
        energy_J = 6.61822698*10**(-22) * u.J
        elo_J = 5.57674801*10**(-22) * u.J
        df = 3
        integrated_flux = 0.672 * u.K * u.km / u.s

        quantities = [t_freq, temp_estimate, lgint300, part300, partition, gu, energy_J,
                      elo_J, df, mol_tag]

        names = ['t_freq', 'temp', 'lgint300', 'partfn300', 'partfn',
                 'dgup', 'eup_J', 'elo_J', 'degfreedom', 'mol_tag']

        mol_data = Phys.from_dict(dict(zip(names, quantities)))

        intl = intensity_conversion(mol_data)

        mol_data.apply([intl.value] * intl.unit, name='intl')

        au = einstein_coeff(mol_data)

        mol_data.apply([au.value] * au.unit, name='eincoeff')

        mol_data.apply([1.] * u.AU * u.AU * u.s, name='beta')
        mol_data.apply([1.] * u.km / (u.m * u.m * u.m), name='cdensity')
        mol_data.apply([1.], name='total_number')

        lte = LTE()

        Q_estimate = 3.594*10**(28) / u.s

        parent = photo_timescale('CO') * vgas

        coma = Haser(Q_estimate, vgas, parent)

        r = 1.06077567 * u.AU

        delta = 0.14633757 * u.AU

        quanteph = [r, delta]

        nameseph = ['r', 'delta']

        ephemobj = Ephem.from_dict(dict(zip(nameseph, quanteph)))

        beta = beta_factor(mol_data, ephemobj)

        mol_data['beta'] = beta

        cdensity = lte.cdensity_Bockelee(integrated_flux, mol_data)

        mol_data['cdensity'] = cdensity

        tnum = total_number(mol_data, aper, b)

        mol_data['total_number'] = tnum

        Q = from_Haser(coma, mol_data, aper=aper)

        qH = np.log10(Q.value)

        assert np.isclose(qH[0], 25.2862, rtol=0.01)
