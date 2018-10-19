
List of Alternative Field Names
===============================

The following table lists alternative field names accepted by `sbpy`
when accessing `~sbpy.data.DataClass` objects, i.e.,
`~sbpy.data.Ephem`, `~sbpy.data.Orbit`, or `~sbpy.data.Phys` objects.

As an example, heliocentric distance can be addressed a ``'r'`` or
``'heldist'``:

    >>> from sbpy.data import Ephem
    >>> ceres = Ephem.from_horizons('Ceres')
    >>> print(ceres['r']) # doctest: +IGNORE_OUTPUT
    [2.69866993] AU
    >>> print(ceres['heldist'])
    [2.69866993] AU

The list of alternative field names is always up to date, but not
complete. The source list is located as
``sbpy.data.conf.fieldnames``. If you think an important alternative
is missing, please suggest by opening an issue. However, keep in mind
that each alternative field name has to be *unique* and *unambiguous*.

