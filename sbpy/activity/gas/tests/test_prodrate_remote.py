# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table

import pytest

from .. import (Haser, photo_timescale, LTE, NonLTE, einstein_coeff,
                intensity_conversion, beta_factor, total_number, from_Haser)
from ....data import Ephem, Phys


class MockPyradex:
    """
    Class to be the mock return value of NonLTE.from_pyradex
    """

    def __init__(self):
        """
        Define a testing dictionary
        """
        self.value = {0: 1.134e14}

# monkeypatched NonLTE.from_pyradex


@pytest.fixture
def mock_nonlte(monkeypatch):
    """
    from_pyradex.value mocked to return dictionary.
    """

    def mock_cdensity(*args, **kwargs):
        """
        Define a testing Quantity
        """
        integrated_flux = args[1]

        if integrated_flux == 0.26 * u.K * u.km / u.s:
            return u.Quantity([2.17469686e+14], 1/u.cm**2)
        elif integrated_flux == 0.27 * u.K * u.km / u.s:
            return u.Quantity([2.25833905e+14], 1/u.cm**2)
        elif integrated_flux == 0.28 * u.K * u.km / u.s:
            return u.Quantity([2.34198124e+14], 1/u.cm**2)
        elif integrated_flux == 1.234 * u.K * u.km / u.s:
            return u.Quantity([1.134e14], 1/u.cm**2)

    monkeypatch.setattr(NonLTE, "from_pyradex", mock_cdensity)


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


@pytest.mark.remote_data
def test_remote_prodrate_simple_hcn():
    pytest.importorskip("astroquery", minversion="0.4.7")

    hcn = Table.read(data_path('HCN.csv'), format="ascii.csv")

    temp_estimate = 47. * u.K
    target = '103P'
    vgas = 0.8 * u.km / u.s
    aper = 30 * u.m  # The aperture for telescope used (Drahus et al. 2012)
    b = 1.13  # Value taken from (Drahus et al. 2012)
    mol_tag = 27001
    transition_freq = (265.886434 * u.GHz).to('MHz')
    q_found = []
    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    intl = intensity_conversion(mol_data)
    mol_data.apply([intl.value] * intl.unit, name='intl')
    au = einstein_coeff(mol_data)
    mol_data.apply([au.value] * au.unit, name='eincoeff')

    for i in range(0, 28):

        time = Time(hcn['Time'][i], format='iso')
        integrated_flux = hcn['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(
            target, epochs=time, id_type='designation',
            closest_apparition=True)

        lte = LTE()

        q = lte.from_Drahus(integrated_flux, mol_data,
                            ephemobj, vgas, aper, b=b)

        q = np.log10(q.value)

        q_found.append(q)

    q_pred = list(hcn['log(Q)'])

    np.testing.assert_almost_equal(q_pred, q_found, decimal=1.3)

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 0.2345)


@pytest.mark.remote_data
def test_remote_prodrate_simple_ch3oh():
    pytest.importorskip("astroquery", minversion="0.4.7")

    ch3oh = Table.read(data_path('CH3OH.csv'), format="ascii.csv")
    temp_estimate = 47. * u.K
    target = '103P'
    vgas = 0.8 * u.km / u.s
    aper = 30 * u.m  # The aperture for telescope used (Drahus et al. 2012)
    b = 1.13  # Value taken from (Drahus et al. 2012)
    mol_tag = 32003
    transition_freq = (157.178987 * u.GHz).to('MHz')
    q_found = []
    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    intl = intensity_conversion(mol_data)
    mol_data.apply([intl.value] * intl.unit, name='intl')
    au = einstein_coeff(mol_data)
    mol_data.apply([au.value] * au.unit, name='eincoeff')

    for i in range(0, 20):

        time = Time(ch3oh['Time'][i], format='iso')
        integrated_flux = ch3oh['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(target, epochs=time,
                                       id_type='designation',
                                       closest_apparition=True)

        lte = LTE()

        q = lte.from_Drahus(integrated_flux, mol_data,
                            ephemobj, vgas, aper, b=b)

        q = np.log10(q.value)

        q_found.append(q)

    q_pred = list(ch3oh['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 0.35)


# Issue #296
# MSK: disabling test as CO and HCN are no longer present in LAMDA database(?)
@pytest.mark.skip
@pytest.mark.remote_data
def test_einstein():
    pytest.importorskip("astroquery", minversion="0.4.7")

    transition_freq_list = [(1611.7935180 * u.GHz).to('MHz'),
                            (177.26111120 * u.GHz).to('MHz')]
    mol_tag_list = [28001, 27001]
    temp_estimate = 300. * u.K

    result = []
    catalog_result = []

    cat = JPLSpec.get_species_table()

    for i in range(2):

        transition_freq = transition_freq_list[i]
        mol_tag = mol_tag_list[i]

        mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
        intl = intensity_conversion(mol_data)
        mol_data.apply([intl.value] * intl.unit, name='intl')

        au = einstein_coeff(mol_data)

        result.append(au.value)

        mol = cat[cat['TAG'] == mol_tag]

        mol_name = mol['NAME'].data[0]

        lam_search = Lamda.query(mol=mol_name.lower())

        lam_result = lam_search[1]

        tran = transition_freq.to('GHz').value

        lam_found = lam_result[lam_result['Frequency'] == tran]

        au_cat = lam_found['EinsteinA']

        au_cat = au_cat.data[0]

        catalog_result.append(au_cat)

    err = (abs((np.array(catalog_result) - np.array(result)) /
               np.array(catalog_result) * 100))

    assert np.all(err < 23.5)


@pytest.mark.remote_data
def test_Haser_prodrate():
    pytest.importorskip("astroquery", minversion="0.4.7")

    co = Table.read(data_path('CO.csv'), format="ascii.csv")

    lte = LTE()
    Q_estimate = 2.8*10**(28) / u.s
    transition_freq = (230.53799 * u.GHz).to('MHz')
    aper = 10 * u.m
    mol_tag = 28001
    temp_estimate = 25. * u.K
    vgas = 0.5 * u.km / u.s
    target = 'C/2016 R2'
    b = 0.74
    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    intl = intensity_conversion(mol_data)
    mol_data.apply([intl.value] * intl.unit, name='intl')
    au = einstein_coeff(mol_data)
    mol_data.apply([au.value] * au.unit, name='eincoeff')
    mol_data.apply([1.] * u.AU * u.AU * u.s, name='beta')
    mol_data.apply([1.] / (u.m * u.m), name='cdensity')
    mol_data.apply([1.], name='total_number')

    q_found = []

    parent = photo_timescale('CO') * vgas
    coma = Haser(Q_estimate, vgas, parent)

    for i in range(0, 5):

        time = Time(co['Time'][i], format='iso')
        integrated_flux = co['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(target, epochs=time)
        beta = beta_factor(mol_data, ephemobj)
        mol_data['beta'] = beta
        cdensity = lte.cdensity_Bockelee(integrated_flux, mol_data)
        mol_data['cdensity'] = cdensity
        tnum = total_number(mol_data, aper, b)
        mol_data['total_number'] = tnum

        Q = from_Haser(coma, mol_data, aper=aper)

        q_found.append(np.log10(Q.value)[0])

    q_pred = list(co['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 2.5)


'''
Last test run: 08/01/2019 11:15:00 , sbpy version: v0.2dev259, python 3.6.8
Author: Giannina Guzman
Tester: Giannina Guzman
Tested: locally, needs pyradex to be installed
Status: Passed
See https://github.com/keflavich/pyradex for installment
'''


@pytest.mark.remote_data
def test_Haser_pyradex():
    pytest.importorskip("pyradex")
    pytest.importorskip("astroquery", minversion="0.4.7")

    co = Table.read(data_path('CO.csv'), format="ascii.csv")

    nonlte = NonLTE()
    lte = LTE()
    Q_estimate = 2.8*10**(28) / u.s
    transition_freq = (230.53799 * u.GHz).to('MHz')
    aper = 10 * u.m
    mol_tag = 28001
    temp_estimate = 25. * u.K
    vgas = 0.5 * u.km / u.s
    target = 'C/2016 R2'
    b = 0.74
    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    intl = intensity_conversion(mol_data)
    mol_data.apply([intl.value] * intl.unit, name='intl')
    au = einstein_coeff(mol_data)
    mol_data.apply([au.value] * au.unit, name='eincoeff')
    mol_data.apply([1.] * u.AU**2 * u.s, name='beta')
    mol_data.apply([1.] / u.m**2, name='cdensity')
    mol_data.apply([1.], name='total_number')

    q_found = []

    parent = photo_timescale('CO') * vgas
    coma = Haser(Q_estimate, vgas, parent)

    for i in range(0, 5):

        time = Time(co['Time'][i], format='iso')
        integrated_flux = co['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(target, epochs=time)
        beta = beta_factor(mol_data, ephemobj)
        mol_data['beta'] = beta
        cdensity_bockelee = lte.cdensity_Bockelee(integrated_flux, mol_data)
        mol_data['cdensity'] = cdensity_bockelee
        cdensity = nonlte.from_pyradex(integrated_flux, mol_data)
        mol_data['cdensity'] = cdensity
        tnum = total_number(mol_data, aper, b)
        mol_data['total_number'] = tnum

        Q = from_Haser(coma, mol_data, aper=aper)

        q_found.append(np.log10(Q.value)[0])

    q_pred = list(co['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 0.35)


@pytest.mark.remote_data
def test_intensity_conversion():
    # test untested case for intensity conversion function
    pytest.importorskip("astroquery", minversion="0.4.7")

    temp_estimate = 47. * u.K
    # vgas = 0.8 * u.km / u.s
    mol_tag = 27001
    transition_freq = (265.886434 * u.GHz).to('MHz')
    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    mol_data['eup_J'] = 3.52359898e-20 * u.J
    mol_data['elo_J'] = 1.76181853e-20 * u.J
    intl = intensity_conversion(mol_data)

    assert np.isclose(intl.value, 6.186509000388917e-11)


@pytest.mark.remote_data
def test_einsteincoeff_case():
    # test untested case for einstein coefficient
    pytest.importorskip("astroquery", minversion="0.4.7")

    temp_estimate = 47. * u.K
    # vgas = 0.8 * u.km / u.s
    mol_tag = 27001
    transition_freq = (265.886434 * u.GHz).to('MHz')
    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    mol_data['t_freq'] = 2658864.34 * u.MHz
    intl = intensity_conversion(mol_data)
    mol_data.apply([intl.value] * intl.unit, name='intl')
    au = einstein_coeff(mol_data)

    assert np.isclose(round(au.value, 4), 0.0086)


@pytest.mark.remote_data
def test_betafactor_case():
    # test untested case for beta beta_factor

    r = 1.06077567 * u.AU
    delta = 0.14633757 * u.AU
    quanteph = [r, delta]
    nameseph = ['r', 'delta']
    ephemobj = Ephem.from_dict(dict(zip(nameseph, quanteph)))
    mol_name = 'CN'
    namephys = ['mol_tag']
    quantphys = [mol_name]
    mol_data = Phys.from_dict(dict(zip(namephys, quantphys)))
    beta = beta_factor(mol_data, ephemobj)

    assert np.isclose(beta.value[0], 354452.18195014383)


'''
Last test run: 08/02/2019 09:32:00 , sbpy version: v0.2dev259, python 3.6.8
Author: Giannina Guzman
Tester: Giannina Guzman
Tested: locally, needs pyradex to be installed
Status: Passed
See https://github.com/keflavich/pyradex for installment
'''


@pytest.mark.remote_data
def test_pyradex_case():
    pytest.importorskip("pyradex")
    pytest.importorskip("astroquery", minversion="0.4.7")

    transition_freq = (177.196 * u.GHz).to(u.MHz)
    mol_tag = 29002
    cdensity_guess = (1.89*10.**(14) / (u.cm * u.cm))
    temp_estimate = 20. * u.K
    temp_back = 2.8 * u.K

    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    mol_data.apply([cdensity_guess.value] *
                   cdensity_guess.unit, name='cdensity')
    mol_data.apply([temp_back.value] * temp_back.unit, name='temp_back')
    mol_data.apply(['HCO+@xpol'], name='lamda_name')
    nonLTE = NonLTE()
    cdensity = nonLTE.from_pyradex(1.234 * u.K * u.km / u.s, mol_data,
                                   iter=100, collider_density={'H2': 900})

    assert np.isclose(cdensity.value[0], 1.134e14)


@pytest.mark.remote_data
def test_Haser_prodrate_pyradex(mock_nonlte):
    pytest.importorskip("pyradex")
    pytest.importorskip("astroquery", minversion="0.4.7")

    co = Table.read(data_path('CO.csv'), format="ascii.csv")

    nonlte = NonLTE()
    lte = LTE()
    Q_estimate = 2.8*10**(28) / u.s
    transition_freq = (230.53799 * u.GHz).to('MHz')
    aper = 10 * u.m
    mol_tag = 28001
    temp_estimate = 25. * u.K
    vgas = 0.5 * u.km / u.s
    target = 'C/2016 R2'
    b = 0.74
    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    intl = intensity_conversion(mol_data)
    mol_data.apply([intl.value] * intl.unit, name='intl')
    au = einstein_coeff(mol_data)
    mol_data.apply([au.value] * au.unit, name='eincoeff')
    mol_data.apply([1.] * u.AU * u.AU * u.s, name='beta')
    mol_data.apply([1.] / (u.m * u.m), name='cdensity')
    mol_data.apply([1.], name='total_number')

    q_found = []

    parent = photo_timescale('CO') * vgas
    coma = Haser(Q_estimate, vgas, parent)

    for i in range(0, 5):

        time = Time(co['Time'][i], format='iso')
        integrated_flux = co['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(target, epochs=time)
        beta = beta_factor(mol_data, ephemobj)
        mol_data['beta'] = beta
        cdensity_bockelee = lte.cdensity_Bockelee(integrated_flux, mol_data)
        mol_data['cdensity'] = cdensity_bockelee
        cdensity = nonlte.from_pyradex(integrated_flux, mol_data)
        mol_data['cdensity'] = cdensity
        tnum = total_number(mol_data, aper, b)
        mol_data['total_number'] = tnum

        Q = from_Haser(coma, mol_data, aper=aper)

        q_found.append(np.log10(Q.value)[0])

    q_pred = list(co['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 0.35)


@pytest.mark.remote_data
def test_pyradex_cdensity(mock_nonlte):
    pytest.importorskip("astroquery", minversion="0.4.7")

    transition_freq = (177.196 * u.GHz).to(u.MHz)
    mol_tag = 29002
    cdensity_guess = (1.89*10.**(14) / (u.cm * u.cm))
    temp_estimate = 20. * u.K
    temp_back = 2.8 * u.K

    mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    mol_data.apply([cdensity_guess.value] *
                   cdensity_guess.unit, name='cdensity')
    mol_data.apply([temp_back.value] * temp_back.unit, name='temp_back')
    mol_data.apply(['HCO+@xpol'], name='lamda_name')
    nonLTE = NonLTE()
    cdensity = nonLTE.from_pyradex(1.234 * u.K * u.km / u.s, mol_data,
                                   iter=100, collider_density={'H2': 900})

    assert np.isclose(cdensity.value[0], 1.134e14)
    assert cdensity.unit == u.cm**-2
