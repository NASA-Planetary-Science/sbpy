import astropy.units as u
from .. import Spectrum
from astropy.tests.helper import remote_data
from astropy.table import Table
import numpy as np

import os


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


@remote_data
def test_remote_prodrate_simple_hcn():

    hcn = Table.read(data_path('HCN.csv'), format="ascii.csv")

    temp_estimate = 33. * u.K

    target = '900918'

    vgas = 0.8 * u.km / u.s

    diameter = 30 * u.m  # The diameter for telescope used (Drahus et al. 2012)

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
                          mol_tag, time, target, vgas, diameter, b=b,
                          id_type='id')

        q = np.log10(q.value)

        q_found.append(q)

    q_pred = list(hcn['log(Q)'])

    np.testing.assert_almost_equal(q_pred, q_found, decimal=1.3)

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    for i in range(0, len(err)):

        assert err[i] < 0.135


@remote_data
def test_remote_prodrate_simple_ch3oh():

    ch3oh = Table.read(data_path('CH3OH.csv'), format="ascii.csv")

    temp_estimate = 33. * u.K

    target = '900918'

    vgas = 0.8 * u.km / u.s

    diameter = 30 * u.m  # The diameter for telescope used (Drahus et al. 2012)

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

    for i in range(0, 21):

        time = ch3oh['Time'][i]

        spectra = ch3oh['T_B'][i] * u.K * u.km / u.s

        s = Spectrum(spectra, dispersionaxis, unit)

        q = s.prodrate_np(spectra, temp_estimate, transition_freq,
                          mol_tag, time, target, vgas, diameter, b=b,
                          id_type='id')

        q = np.log10(q.value)

        q_found.append(q)

    q_pred = list(ch3oh['log(Q)'])

    err = abs((np.array(q_pred) - np.array(q_found)) / np.array(q_pred) * 100)

    for i in range(0, len(err)):

        assert err[i] < 0.39
