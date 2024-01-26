# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy utils module

created on March 12, 2019

"""

__all__ = ["required_packages", "optional_packages"]

from importlib import import_module
from warnings import warn
from ..exceptions import RequiredPackageUnavailable, OptionalPackageUnavailable


def required_packages(*packages, message=None):
    """Verifies the arguments are valid packages.


    Parameters
    ----------
    *modules : str
        The names of packages to test.

    message : str
        An additional message to show the user, e.g., the reason for the
        requirement.

    Raises
    ------

    `~sbpy.exceptions.RequiredPackageUnavailable`


    Examples
    --------

    >>> from sbpy.utils import required_packages
    >>> required_packages("unavailable_package")
    Traceback (most recent call last):
    ...
    sbpy.exceptions.RequiredPackageUnavailable: `unavailable_package` is required.

    """

    for package in packages:
        try:
            import_module(package)
        except ModuleNotFoundError as exc:
            _message = "" if message is None else "  " + message
            raise RequiredPackageUnavailable(
                f"`{package}` is required.{_message}"
            ) from None


def optional_packages(*packages, message=None):
    """Decorator that warns if the arguments are not valid packages.


    Parameters
    ----------
    *modules : str
        The names of packages to test.

    message : str
        An additional note to show the user, e.g., a description of the fallback
        behavior.


    Returns
    -------
    success : bool
        ``True`` if all optional packages are available.


    Warnings
    --------

    `~sbpy.exceptions.OptionalPackageUnavailable`


    Examples
    --------

    >>> from sbpy.utils import optional_packages
    >>> optional_packages("unavailable_package")  # doctest: +SKIP
    OptionalPackageUnavailable: Optional package `unavailable_package` is unavailable.

    """
    # the doctest line is skipped to avoid polluting the testing suite with a warning

    for package in packages:
        try:
            import_module(package)
        except ModuleNotFoundError:
            _message = "" if message is None else "  " + message
            warn(
                f"`{package}`.{_message}",
                OptionalPackageUnavailable,
            )
            return False
    return True
