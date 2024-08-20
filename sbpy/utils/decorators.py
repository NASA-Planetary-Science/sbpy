# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Common sbpy method/function decorators."""

__all__ = ["requires", "optionally_uses"]

from typing import Callable
from functools import wraps
from astropy.time import Time
from . import core
from ..exceptions import RequiredPackageUnavailable


def time_input(func: Callable, **kwargs):
    """Decorator that validates time inputs.


    Examples
    --------

    from sbpy.utils import time_input
    @time_input(t=Time)
    def myfunction(t):
        return t.mjd

    from sbpy.utils import time_input
    @time_input(t=Time)
    def myfunction(t):
        return t.mjd


    """


def requires(*packages, message=None):
    """Decorator that verifies the arguments are valid packages.


    Parameters
    ----------
    *packages : str
        The names of packages to test.


    Raises
    ------

    `~sbpy.exceptions.RequiredPackageUnavailable`


    Examples
    --------

    >>> from sbpy.utils.decorators import requires
    >>> @requires("unavailable_package")
    ... def f():
    ...     pass
    >>> f()
    Traceback (most recent call last):
    ...
    sbpy.exceptions.RequiredPackageUnavailable: `unavailable_package` is required. (sbpy.utils.decorators.f)

    """

    def decorator(wrapped_function):
        function_name = ".".join(
            (wrapped_function.__module__, wrapped_function.__qualname__)
        )
        _message = ("" if message is None else f"{message} ") + f"({function_name})"

        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            try:
                core.required_packages(*packages, message=_message)
            except RequiredPackageUnavailable as exc:
                # trim a couple levels of the traceback to clean up error messages
                raise exc.with_traceback(exc.__traceback__.tb_next.tb_next)
            return wrapped_function(*func_args, **func_kwargs)

        return wrapper

    return decorator


def optionally_uses(*packages, message=None):
    """Decorator that warns if the arguments are not valid modules.


    Parameters
    ----------
    *packages : str
        The names of packages to test.

    message : str
        An additional message to show the user, e.g., a description of the
        fallback behavior.


    Warnings
    --------

    `~sbpy.exceptions.OptionalPackageUnavailable`


    Examples
    --------

    >>> from sbpy.utils.decorators import optionally_uses
    >>> @optionally_uses("unavailable_package")
    ... def f():
    ...     pass
    >>> f()  # doctest: +SKIP
    sbpy/utils/decorators.py::sbpy.utils.decorators.optional
    ...
    OptionalPackageUnavailable: Optional package `unavailable_package` is unavailable. (sbpy.utils.decorators.f)

    """
    # the doctest line is skipped to avoid polluting the testing suite with a warning

    def decorator(wrapped_function):
        function_name = ".".join(
            (wrapped_function.__module__, wrapped_function.__qualname__)
        )
        _message = ("" if message is None else f"{message} ") + f"({function_name})"

        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            core.optional_packages(*packages, message=_message)
            return wrapped_function(*func_args, **func_kwargs)

        return wrapper

    return decorator
