# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Common sbpy method/function decorators."""

from functools import wraps
from ..exceptions import RequiredPackageUnavailable

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

    >>> try:
    ...     import unavailable_package
    ... except ImportError:
    ...     unavailable_package = None
    >>>
    >>> from sbpy.utils.decorators import requires
    >>>
    >>> @requires(unavailable_package=unavailable_package)
    ... def f():
    ...     pass
    >>>
    >>> f()
    Traceback (most recent call last):
    ...
    sbpy.exceptions.RequiredPackageUnavailable: unavailable_package is required for f.

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
