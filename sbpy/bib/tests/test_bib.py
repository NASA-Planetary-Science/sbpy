# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import pytest
from ..core import *


# skip function tests utilizing ads.ExportQuery for now as it is unstable

# get file path of a static data file for testing
# def data_path(filename):
#     data_dir = os.path.join(os.path.dirname(__file__), 'data')
#     return os.path.join(data_dir, filename)

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


def test_register_single():
    reset()
    with Tracking():
        register('test1', {'track_this': 'bibcode1'})

    assert (set(['sbpy:', 'software:', '2019JOSS....4.1426M',
                 'test1:', 'track_this:', 'bibcode1'])
            == set(show().split()))
    stop()
    reset()


def test_register_list():
    reset()
    with Tracking():
        register('test1', {'track_this': ['bibcode1', 'bibcode2']})

    assert (set(['sbpy:', 'software:', '2019JOSS....4.1426M',
                 'test1:', 'track_this:', 'bibcode1', 'bibcode2'])
            == set(show().split()))
    stop()
    reset()


def test_register_double():
    reset()
    with Tracking():
        register('test1', {'track_this': ['bibcode1', 'bibcode2']})
        register('test1', {'track_this': ['bibcode2']})
        register('test1', {'track_this': ['bibcode3']})

    assert show().count('bibcode2') == 1
    stop()
    reset()


def test_Tracking():
    reset()
    with Tracking():
        assert status()
        register('test1', {'track_this': 'bibcode1'})
        register('test1', {'track_this': 'bibcode2'})
        register('test1', {'track_this_too': 'bibcode'})
        register('test2', {'track_this': 'bibcode'})
        register('test3', {'track_this': 'bibcode',
                           'and_track_that': 'bibcode'})
    assert not status()

    register('test', {'do not track this': 'bibcode'})
    assert set(['sbpy:', 'software:', '2019JOSS....4.1426M',
                'test1:', 'track_this:', 'bibcode1', 'bibcode2',
                'track_this_too:', 'bibcode', 'test2:', 'track_this:',
                'bibcode', 'test3:', 'track_this:', 'bibcode',
                'and_track_that:', 'bibcode']) == set(show().split())
    # different builds will have different orders for bibcode 1 and 2, to
    # avoid the build failing because of this we use sets
    stop()
    reset()


def test_Tracking_issue_64():
    from sbpy.activity import photo_lengthscale
    reset()
    with Tracking():
        gamma_H2O = photo_lengthscale('H2O')
        gamma_OH = photo_lengthscale('OH')
    words = show().split()
    assert 'OH' in words
    assert 'H2O' in words
    stop()
    reset()


def test_Tracking_reporter(capsys):
    reset()
    with Tracking(reporter=show):
        register('test1', {'track_this': 'bibcode1'})
    captured = capsys.readouterr()
    assert (set(['sbpy:', 'software:', '2019JOSS....4.1426M',
                 'test1:', 'track_this:', 'bibcode1'])
            == set(captured.out.split()))
    stop()
    reset()


def test_cite_function():

    @cite({'method': '1687pnpm.book.....N'})
    def force(mass, acceleration):
        return mass * acceleration

    reset()
    track()
    force(1, 2)
    assert (set(
        ['sbpy:', 'software:', '2019JOSS....4.1426M',
         'sbpy.bib.tests.test_bib.test_cite_function.<locals>.force:',
         'method:', '1687pnpm.book.....N'])
        == set(show().split()))
    stop()
    reset()


def test_cite_function_twice():

    @cite({'method': '1687pnpm.book.....N'})
    @cite({'interpretation': 'philosophical reference'})
    def force(mass, acceleration):
        return mass * acceleration

    reset()
    track()
    force(1, 2)
    assert (set(
        ['sbpy:', 'software:', '2019JOSS....4.1426M',
         'sbpy.bib.tests.test_bib.test_cite_function_twice.<locals>.force:',
         'method:', '1687pnpm.book.....N', 'interpretation:',
         'philosophical', 'reference'])
        == set(show().split()))
    stop()
    reset()


def test_cite_class_method():
    reset()

    class Physics:
        @staticmethod
        @cite({'method': '1687pnpm.book.....N'})
        def force(mass, acceleration):
            return mass * acceleration

    with Tracking():
        p = Physics()
        p.force(1, 2)

    assert (set([
        'sbpy:', 'software:', '2019JOSS....4.1426M',
        'sbpy.bib.tests.test_bib.test_cite_class_method'
        '.<locals>.Physics.force:',
        'method:', '1687pnpm.book.....N'])
        == set(show().split()))
    reset()


def test_cite_class():
    reset()

    @cite({'method': '1687pnpm.book.....N'})
    class Force:
        def __call__(self, mass, acceleration):
            return mass * acceleration

    with Tracking():
        f = Force()
        f(1, 2)

    assert (set([
        'sbpy:', 'software:', '2019JOSS....4.1426M',
        'sbpy.bib.tests.test_bib.test_cite_class'
        '.<locals>.Force:', 'method:', '1687pnpm.book.....N'])
        == set(show().split()))
    reset()


def test_filter():
    reset()
    with Tracking():
        register('test1', {'track_this': 'bibcode1'})
        register('test1', {'software': 'bibcode2'})
        register('test1', {'track_this_too': 'bibcode'})
        register('test2', {'software': 'bibcode'})
        register('test3', {'track_this': 'bibcode',
                           'software': 'bibcode'})

    assert set(['sbpy:', 'software:', '2019JOSS....4.1426M',
                'test1:', 'software:', 'bibcode2',
                'test2:', 'software:', 'bibcode',
                'test3:', 'software:',
                'bibcode']) == set(show(filter='software').split())
    # different builds will have different orders for bibcode 1 and 2, to
    # avoid the build failing because of this we use sets
    stop()
    reset()
