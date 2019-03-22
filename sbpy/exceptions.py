# Licensed under a 3-clause BSD style license - see LICENSE.rst


class SbpyException(Exception):
    pass


class SinglePointSpectrumError(SbpyException):
    '''Single point provided, but multiple values expected.'''
    pass
