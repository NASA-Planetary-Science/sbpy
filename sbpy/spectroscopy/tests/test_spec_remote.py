import astropy.units as u
from ..spectroscopy import prodrate_np
from astropy.tests.helper import remote_data
from astropy.table import Table
import numpy as np


@remote_data
def test_remote_prodrate_simple():

    hcn = Table.read('data/HCN.csv', format="ascii.csv")

    temp_estimate = 33. * u.K

    target = '900918'

    vgas = 0.8 * u.km / u.s

    diameter = 30 * u.m  # The diameter for telescope used (Drahus et al. 2012)

    b = 1.13  # Value taken from (Drahus et al. 2012)

    mol_tag = 27001

    transition_freq = 265.886434 * u.GHz

    q_found = []

    for i in range(0, 28):

        time = hcn['Time'][i]

        spectra = hcn['T_B'][i] * u.K * u.km / u.s

        q = prodrate_np(spectra, temp_estimate, transition_freq, mol_tag,
                        time, target, vgas, diameter, b=b, id_type='id')

        q = np.log10(q.value)

        q_found.append(q)

    q_pred = list(hcn['log(Q)'])

    q_pred_upperbound = list(hcn['log(Q)']+hcn['error_top'])

    q_pred_lowerbound = list(hcn['log(Q)']-hcn['error_bottom'])

    np.testing.assert_almost_equal(q_pred, q_found, decimal=1.3)

    assert q_pred_lowerbound < q_found < q_pred_upperbound
