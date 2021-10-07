# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""sbpy Data Decorators

Based on astropy's units decorator: `~astropy.units.quantity_input`.

"""

__all__ = [
    'quantity_to_dataclass',
    'dataclass_input'
]

import inspect
from functools import wraps
from astropy.table import Table, QTable
from astropy.time import Time
import astropy.units as u
from astropy.units.decorators import _get_allowed_units, _validate_arg_value
from .core import DataClass, DataClassError
from . import Conf


def quantity_to_dataclass(**kwargs):
    """Decorator that converts astropy quantities to sbpy data classes.

    Use this decorator when your function is based on a single field
    in an sbpy `~sbpy.data.DataClass`.


    Examples
    --------

    This function accepts `~sbpy.data.Ephem` objects, but only uses
    heliocentric distance:

    >>> import astropy.units as u
    >>> import sbpy.data as sbd
    >>>
    >>> @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'rh'))
    ... def temperature(eph):
    ...     return 278 * u.K / (eph['rh'] / u.au)**0.5
    >>>
    >>> print(temperature(1 * u.au))    # doctest: +FLOAT_CMP
    [278.] K
    >>> eph = sbd.Ephem.from_dict({'rh': 1 * u.au})
    >>> print(temperature(eph))         # doctest: +FLOAT_CMP
    [278.] K

    This decorator also validates the dimensions of function parameters
    against the default dimensions as listed in the Field Name List
    (https://sbpy.readthedocs.io/en/latest/sbpy/data/fieldnames.html#id1).
    Users can provide equivalencies through an optional parameter
    `equivalencies=` to be used in unit checking.  Equivalencies for
    dimensionless angle and temperature are automatically enabled.

    A `~astropy.units.UnitsError` will be raised if the unit attribute of the
    argument is not equivalent to default unit.    If the default dimension is
    not `None`, and the argument has no unit attribute, and  i.e. it is not a
    Quantity object, a `ValueError` will be raised.
    """
    equivalencies = kwargs.pop('equivalencies', [])

    def decorator(wrapped_function):
        decorator_kwargs = kwargs  # for clarity

        # Extract the function signature for the function we are wrapping.
        wrapped_signature = inspect.signature(wrapped_function)

        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            # Bind the arguments of our new function to the signature
            # of the original.
            bound_args = wrapped_signature.bind(*func_args, **func_kwargs)

            for param in wrapped_signature.parameters.values():
                # is this a parameter that we might want to replace?
                if param.name not in decorator_kwargs:
                    # no
                    continue

                # bind relied on a default value
                if (param.name not in bound_args.arguments and
                        param.default is not param.empty):
                    bound_args.arguments[param.name] = param.default

                # check passed argument, update as needed
                arg = bound_args.arguments[param.name]

                # argument value is None, and the default value is
                # None, pass through the None
                if arg is None and param.default is None:
                    continue

                # get requested DataClass and field name
                dataclass, field = None, None
                for v in decorator_kwargs[param.name][:2]:
                    if isinstance(v, str):
                        field = v
                    elif issubclass(v, DataClass):
                        dataclass = v

                if any((dataclass is None, field is None)):
                    raise ValueError(
                        'quantity_to_dataclass decorator requires a '
                        'DataClass object and a field name as a string.')

                # check field name
                if field not in Conf.fieldname_idx:
                    raise DataClassError("argument '{}' to function '{}' has"
                        " an invalid field name '{}' for {}"
                        " object".format(param.name, wrapped_function.__name__,
                        field, dataclass))

                if not isinstance(arg, dataclass):
                    # Argument is not a DataClass.  Make it so.
                    new_arg = dataclass.from_dict({field: arg})
                    bound_args.arguments[param.name] = new_arg

            return wrapped_function(*bound_args.args, **bound_args.kwargs)
        return wrapper
    return decorator


class DataClassInput:
    @classmethod
    def as_decorator(cls, func=None, **kwargs):
        """Decorator that converts parameters to `DataClass`.

        sbpy methods use ``DataClass`` objects whenever possible.  But for
        convenience, we may let users pass other objects that are
        internally converted:

            * dictionary,
            * file name,
            * `astropy.table.Table` or `~astropy.table.QTable`.


        Examples
        --------

        >>> import astropy.units as u
        >>> import sbpy.data as sbd
        >>>
        >>> @sbd.dataclass_input(eph=sbd.Ephem)
        ... def myfunction(eph):
        ...     return eph['rh']**2 * eph['delta']**2
        >>>
        >>> dictionary = {'rh': 2 * u.au, 'delta': 1 * u.au}
        >>> print(myfunction(dictionary))    # doctest: +FLOAT_CMP
        [4.0] AU4
        >>>
        >>> from astropy.table import QTable
        >>> qtable = QTable([[2] * u.au, [1] * u.au], names=('rh', 'delta'))
        >>> print(myfunction(qtable))        # doctest: +FLOAT_CMP
        [4.0] AU4

        Data classes may also be specified with function annotations:
        >>> import sbpy.data as sbd
        >>>
        >>> @sbd.dataclass_input
        ... def myfunction(eph: sbd.Ephem):
        ...     return eph['rh']**2 * eph['delta']**2

        """
        self = cls(**kwargs)
        if func is not None and not kwargs:
            return self(func)
        else:
            return self

    def __init__(self, func=None, **kwargs):
        self.decorator_kwargs = kwargs

    def __call__(self, wrapped_function):
        # Extract the function signature for the function we are
        # wrapping.
        wrapped_signature = inspect.signature(wrapped_function)

        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            # Bind the arguments of our new function to the signature
            # of the original.
            bound_args = wrapped_signature.bind(*func_args,
                                                **func_kwargs)

            for param in wrapped_signature.parameters.values():
                # is this a parameter that we might want to replace?
                if param.name in self.decorator_kwargs:
                    target = self.decorator_kwargs[param.name]
                else:
                    target = param.annotation

                # not in decorator_kwargs and not annotated
                if target is inspect.Parameter.empty:
                    continue

                # bind relied on a default value
                if (param.name not in bound_args.arguments and
                        param.default is not param.empty):
                    bound_args.arguments[param.name] = param.default

                arg = bound_args.arguments[param.name]

                # argument value is None, and the default value is
                # None, pass through the None
                if arg is None and param.default is None:
                    continue

                # not a DataClass?  carry on.
                try:
                    if issubclass(target, DataClass):
                        dataclass = target
                    else:
                        continue
                except TypeError:
                    continue

                # check passed argument, update as needed
                if isinstance(arg, dict):
                    new_arg = dataclass.from_dict(arg)
                elif isinstance(arg, (Table, QTable)):
                    new_arg = dataclass.from_table(arg)
                elif isinstance(arg, str):
                    new_arg = dataclass.from_file(arg)
                else:
                    continue

                bound_args.arguments[param.name] = new_arg

            return wrapped_function(*bound_args.args,
                                    **bound_args.kwargs)
        return wrapper


dataclass_input = DataClassInput.as_decorator
