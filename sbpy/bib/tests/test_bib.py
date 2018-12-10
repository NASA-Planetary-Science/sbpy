# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import pytest
import time
from ...thermal import NEATM
from .. import (track, stop, register, reset, to_text, to_bibtex,
                to_aastex, to_icarus, to_mnras, Tracking)


# get file path of a static data file for testing
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)


# skip function tests utilizing ads.ExportQuery for now as it is unstable

# @pytest.mark.remote_data
# def test_text():
#     reset()
#     track()
#     neatm = NEATM()
#     assert (['sbpy.thermal.NEATM:', 'method:', 'Harris', '1998',
#              'Icarus', 'Vol',  '131', '2', '291'] ==
#             to_text().replace(',', '').split())
#     reset()
#     stop()
#     time.sleep(1)


# @pytest.mark.remote_data
# def test_bibtex():
#     reset()
#     track()
#     neatm = NEATM()
#     with open(data_path('neatm.bib')) as bib_file:
#         assert to_bibtex().strip() == bib_file.read().strip()
#     reset()
#     stop()
#     time.sleep(1)


# @pytest.mark.remote_data
# def test_aastex():
#     reset()
#     track()
#     register('faketask', {'fakesubtask': '2018ApJS..238...22H'})
#     with open(data_path('hora.aas')) as aas_file:
#         assert to_aastex().strip() == aas_file.read().strip()
#     reset()
#     stop()
#     time.sleep(1)


# @pytest.mark.remote_data
# def test_icarus():
#     reset()
#     track()
#     register('faketask', {'fakesubtask': '1996DPS....28.2504G'})
#     with open(data_path('giorgini.icar')) as icar_file:
#         assert to_icarus().strip() == icar_file.read().strip()
#     reset()
#     stop()
#     print(to_text().split())


def test_Tracking():
    reset()

    with Tracking():
        register('test1', {'track_this': 'bibcode1'})
        register('test1', {'track_this': 'bibcode2'})
        register('test1', {'track_this_too': 'bibcode'})
        register('test2', {'track_this': 'bibcode'})
        register('test3', {'track_this': 'bibcode',
                           'and_track_that': 'bibcode'})

    register('test', {'do not track this': 'bibcode'})
    assert set(['test1:', 'track_this:', 'bibcode1', 'bibcode2',
                'track_this_too:', 'bibcode', 'test2:', 'track_this:',
                'bibcode', 'test3:', 'track_this:', 'bibcode',
                'and_track_that:', 'bibcode']) == set(to_text().split())
    # different builds will have different orders for bibcode 1 and 2, to
    # avoid the build failing because of this we use sets


def test_to_text_filter():

    with Tracking():
        register('test1', {'track_this': 'bibcode1'})
        register('test1', {'software': 'bibcode2'})
        register('test1', {'track_this_too': 'bibcode'})
        register('test2', {'software': 'bibcode'})
        register('test3', {'track_this': 'bibcode',
                           'software': 'bibcode'})

    assert set(['test1:', 'software:', 'bibcode2',
                'test2:', 'software:', 'bibcode',
                'test3:', 'software:',
                'bibcode']) == set(to_text(filter='software').split())
    # different builds will have different orders for bibcode 1 and 2, to
    # avoid the build failing because of this we use sets


def test_Tracking_issue_64():
    from sbpy.activity import photo_lengthscale
    reset()
    with Tracking():
        gamma_H2O = photo_lengthscale('H2O')
        gamma_OH = photo_lengthscale('OH')
    words = to_text().split()
    assert 'OH' in words
    assert 'H2O' in words
