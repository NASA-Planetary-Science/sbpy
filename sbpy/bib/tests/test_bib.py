# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import pytest
from ...thermal import NEATM
from .. import track, stop, register, reset, to_text, to_bibtex, Tracking


# get file path of a static data file for testing
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

# deactivate remote tests for now: not sure about how to handle token

# @pytest.mark.remote_data
# def test_text():
#     reset()
#     track()
#     neatm = NEATM()
#     print(to_text())
#     assert ['sbpy.thermal.NEATM:', 'method:', 'Harris', '1998,',
#             '1998Icar..131..291H'] == to_text().split()
#     reset()
#     stop()


# @pytest.mark.remote_data
# def test_bibtex():
#     reset()
#     track()
#     neatm = NEATM()
#     with open(data_path('neatm.bib')) as bib_file:
#         assert to_bibtex() == bib_file.read()
#     reset()
#     stop()


def test_Tracking():
    reset()

    with Tracking():
        register('test1', {'track_this': 'bibcode1'})
        register('test1', {'track_this': 'bibcode2'})
        register('test1', {'track_this_too': 'bibcode'})
        register('test2', {'track_this': 'bibcode'})

    register('test', {'do not track this': 'bibcode'})
    assert set(['test1:', 'track_this:', 'bibcode1', 'bibcode2', 'track_this_too:',
                'bibcode', 'test2:', 'track_this:',
                'bibcode']) == set(to_text().split())  #different builds will have different orders for bibcode 1 and 2, to avoid the build failing because of this we use sets
    print(to_text().split())


def test_Tracking_issue_64():
    from sbpy.activity import photo_lengthscale
    reset()
    with Tracking():
        gamma_H2O = photo_lengthscale('H2O')
        gamma_OH = photo_lengthscale('OH')
    words = to_text().split()
    assert 'OH' in words
    assert 'H2O' in words
