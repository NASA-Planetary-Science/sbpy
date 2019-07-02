import os

import copy
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
from astropy.table import Table
from astroquery.lamda import Lamda
from astroquery.jplspec import JPLSpec
from .. import Haser, photo_timescale
from ....data import Ephem, Phys
from .. import LTE, einstein_coeff, intensity_conversion, beta_factor, total_number_nocd


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


@remote_data
def test_remote_prodrate_simple_hcn():

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
    mol_data.add_column([intl.value] * intl.unit,
                        name='Integrated line intensity at desired temp')
    au = einstein_coeff(mol_data)
    mol_data.add_column([au.value] * au.unit, name='eincoeff')

    for i in range(0, 28):

        time = Time(hcn['Time'][i], format='iso')
        integrated_flux = hcn['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(target, epochs=time.jd, id_type='id')

        lte = LTE()

        q = lte.from_Drahus(integrated_flux, mol_data, ephemobj, vgas, aper, b=b)

        q = np.log10(q.value)

        q_found.append(q)

    q_pred = list(hcn['log(Q)'])

    np.testing.assert_almost_equal(q_pred, q_found, decimal=1.3)

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 0.2345)


@remote_data
def test_remote_prodrate_simple_ch3oh():

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
    mol_data.add_column([intl.value] * intl.unit,
                        name='Integrated line intensity at desired temp')
    au = einstein_coeff(mol_data)
    mol_data.add_column([au.value] * au.unit, name='eincoeff')

    for i in range(0, 20):

        time = Time(ch3oh['Time'][i], format='iso')
        integrated_flux = ch3oh['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(target, epochs=time.jd, id_type='id')

        lte = LTE()

        q = lte.from_Drahus(integrated_flux, mol_data, ephemobj, vgas, aper, b=b)

        q = np.log10(q.value)

        q_found.append(q)

    q_pred = list(ch3oh['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 0.35)


@remote_data
def test_einstein():

    transition_freq_list = [(1611.7935180 * u.GHz).to('MHz'),
                            (177.26111120 * u.GHz).to('MHz')]
    mol_tag_list = [28001, 27001]
    temp_estimate = 300. * u.K

    result = []
    catalog_result = []

    cat = JPLSpec.get_species_table()

    for i in range(0, 2):

        transition_freq = transition_freq_list[i]
        mol_tag = mol_tag_list[i]

        mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
        intl = intensity_conversion(mol_data)
        mol_data.add_column([intl.value] * intl.unit,
                            name='Integrated line intensity at desired temp')

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


@remote_data
def test_Haser_prodrate():

    co = Table.read(data_path('CO.csv'), format="ascii.csv")

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
    mol_data.add_column([intl.value] * intl.unit,
                        name='Integrated line intensity at desired temp')
    au = einstein_coeff(mol_data)
    mol_data.add_column([au.value] * au.unit, name='eincoeff')
    mol_data.add_column([1.] * u.AU * u.AU * u.s, name='beta')
    mol_data.add_column([1.], name='total_number_nocd')

    q_found = []

    parent = photo_timescale('CO') * vgas
    coma = Haser(Q_estimate, vgas, parent)

    for i in range(0, 5):

        time = Time(co['Time'][i], format='iso')
        integrated_flux = co['T_B'][i] * u.K * u.km / u.s
        ephemobj = Ephem.from_horizons(target, epochs=time.jd)
        beta = beta_factor(mol_data, ephemobj)
        mol_data['beta'] = beta
        tnum = total_number_nocd(integrated_flux, mol_data, aper, b)

        mol_data['total_number_nocd'] = tnum

        lte = LTE()

        Q = lte.from_Haser(coma, mol_data, aper=aper)

        q_found.append(np.log10(Q.value)[0])

    q_pred = list(co['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 2.5)
