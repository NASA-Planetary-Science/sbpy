# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Common sbpy method/function decorators."""

__all__ = ["requires", "optional"]

from warnings import warn
from functools import wraps
from ..exceptions import RequiredPackageUnavailable, OptionalPackageUnavailable

def requires(*packages):
    """Decorator that verifies the arguments are valid packages.

    Parameters
    ----------
    *modules : str
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
    sbpy.exceptions.RequiredPackageUnavailable: unavailable_package is required for sbpy.utils.decorators.f.

    """

    def decorator(wrapped_function):
        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            for package in packages:
                try:
                    __import__(package)
                except ImportError:
                    function_name = '.'.join((wrapped_function.__module__,
                                              wrapped_function.__qualname__))
                    raise RequiredPackageUnavailable(
                        f"{package} is required for {function_name}."
                    )
            return wrapped_function(*func_args, **func_kwargs)
        
        return wrapper

    return decorator


def optional(*packages, note=None):
    """Decorator that warns if the arguments are not valid packages.

    Parameters
    ----------
    *modules : str
        The names of packages to test.

    note : str
        An additional note to show the user, e.g., a description of the fallback
        behavior.


    Warnings
    --------

    `~sbpy.exceptions.OptionalPackageUnavailable`


    Examples
    --------

    >>> from sbpy.utils.decorators import requires
    >>> @optional("unavailable_package")
    ... def f():
    ...     pass
    >>> f()  # doctest: +IGNORE_OUTPUT
    sbpy/utils/decorators.py::sbpy.utils.decorators.optional
    ...
    OptionalPackageUnavailable: Optional package unavailable_package is not available for sbpy.utils.decorators.f.

    """

    def decorator(wrapped_function):
        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            for package in packages:
                try:
                    __import__(package)
                except ImportError:
                    function_name = '.'.join((wrapped_function.__module__,
                                              wrapped_function.__qualname__))
                    warn(OptionalPackageUnavailable(
                        f"Optional package {package} is not available "
                        f"for {function_name}."
                        + ("" if note is None else note)
                    ))
            return wrapped_function(*func_args, **func_kwargs)
        
        return wrapper

    return decorator
