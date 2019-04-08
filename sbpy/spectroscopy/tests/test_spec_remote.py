import os

import numpy as np
import astropy.units as u
from astropy.tests.helper import remote_data
from astropy.table import Table
from astroquery.lamda import Lamda
from astroquery.jplspec import JPLSpec
from ...activity.gas import Haser, photo_timescale
from .. import Spectrum, einstein_coeff


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
    dispersionaxis = 1
    unit = u.Hz

    for i in range(0, 28):

        time = hcn['Time'][i]
        spectra = hcn['T_B'][i] * u.K * u.km / u.s

        s = Spectrum(spectra, dispersionaxis, unit)

        q = s.prodrate_np(spectra, temp_estimate, transition_freq,
                          mol_tag, time, target, vgas, aper, b=b,
                          id_type='id')

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
    transition_freq = [(157.178987 * u.GHz).to('MHz'),
                       (157.246062 * u.GHz).to('MHz'),
                       (157.270832 * u.GHz).to('MHz'),
                       (157.272338 * u.GHz).to('MHz'),
                       (157.276019 * u.GHz).to('MHz')]
    q_found = []
    dispersionaxis = 1
    unit = u.Hz

    for i in range(0, 20):

        time = ch3oh['Time'][i]
        spectra = ch3oh['T_B'][i] * u.K * u.km / u.s

        s = Spectrum(spectra, dispersionaxis, unit)

        q = s.prodrate_np(spectra, temp_estimate, transition_freq,
                          mol_tag, time, target, vgas, aper, b=b,
                          id_type='id')

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

        au = einstein_coeff(temp_estimate, transition_freq, mol_tag)

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

    spec = Spectrum(0.26 * u.K * u.km / u.s, 1, u.Hz)

    Q_estimate = 2.8*10**(28) / u.s
    transition_freq = (230.53799 * u.GHz).to('MHz')
    aper = 10 * u.m
    mol_tag = 28001
    temp_estimate = 25. * u.K
    vgas = 0.5 * u.km / u.s
    target = 'C/2016 R2'
    b = 0.74

    q_found = []

    for i in range(0, 5):

        time = co['Time'][i]
        spectra = co['T_B'][i] * u.K * u.km / u.s

        parent = photo_timescale('CO') * vgas

        coma = Haser(Q_estimate, vgas, parent)

        Q = spec.production_rate(coma, spectra, temp_estimate,
                                 transition_freq, mol_tag, time, target,
                                 aper=aper, b=b)

        q_found.append(np.log10(Q.value)[0])

    q_pred = list(co['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    assert np.all(err < 2.5)
