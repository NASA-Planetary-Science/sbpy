# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy exceptions and warnings

General exceptions/warnings for all of sbpy are specified below.
Exceptions and warnings that are sub-module specific should be in the
respective sub-module.

"""


class SbpyException(Exception):
    "Exception base class for all sbpy exceptions."


class RequiredPackageUnavailable(SbpyException):
    "Required package is not available."


class SbpyWarning(Warning):
    "Warning base class for all sbpy warnings."


class OptionalPackageUnavailable(SbpyWarning):
    "Optional package is not available."


class TestingNeeded(SbpyWarning):
    "More testing is needed to understand the issue."
