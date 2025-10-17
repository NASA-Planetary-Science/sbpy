# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
sbpy: The Small Bodies Python package
"""

from importlib.metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass
