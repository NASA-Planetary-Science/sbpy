=============
 Using Names
=============

`~sbpy.data.Names` is different from the other classes in `~sbpy.data`
in that it does not use `~sbpy.data.DataClass` as a base class. Instead,
`~sbpy.data.Names` does not contain any data, it merely serves as an
umbrella for functions to identify asteroid and comet names, numbers,
and designations.

In order to distinguish if a string designates a comet or an asteroid,
you can use the following code:

    >>> from sbpy.data import Names
    >>> Names.asteroid_or_comet('(1) Ceres')
    'asteroid'
    >>> Names.asteroid_or_comet('2P/Encke')
    'comet'

The module basically uses regular expressions to match the input
strings and find patterns that agree with asteroid and comet names,
numbers, and designations. There are separate tasks to identify
asteroid and comet identifiers:

    >>> Names.parse_asteroid('(228195) 6675 P-L') # doctest: +SKIP
    {'number': 228195, 'desig': '6675 P-L'}
    >>> Names.parse_asteroid('C/2001 A2-A (LINEAR)') # doctest: +SKIP
    ... sbpy.data.names.TargetNameParseError: C/2001 A2-A (LINEAR) does not appear to be an asteroid identifier
    >>> Names.parse_comet('12893') # doctest: +SKIP
    ... sbpy.data.names.TargetNameParseError: 12893 does not appear to be a comet name
    >>> Names.parse_comet('73P-C/Schwassmann Wachmann 3 C') # doctest: +SKIP
    {'type': 'P', 'number': 73, 'fragment': 'C', 'name': 'Schwassmann Wachmann 3 C'}

In order to be able to distinguish between asteroid and comet
identifiers, `sbpy` follows the MPC guideline in that it requires
comet identifiers to include the comet type in either in combination
with a number (e.g., ``'259P'``), a name (e.g., ``'P/Halley'``), or
both (e.g., ``'2P/Encke'``). For instance, the identifier ``'Halley'``
would be identified as an asteroid, as it lacks a comet type
identifier. Hence, some caution is advised when using these routines -
identification might not be unambiguous.

Sorting names with a natural sort order
---------------------------------------

Sorting with Python's built-in functions might not return the desired
order:

    >>> comets = ['9P/Tempel 1',
    ...           '101P/Chernykh',
    ...           '10P/Tempel 2',
    ...           '2P/Encke']
    >>> sorted(comets)
    ['101P/Chernykh', '10P/Tempel 2', '2P/Encke', '9P/Tempel 1']

101P and 10P are placed at the start of the list because Python is
performing a string comparison, which is character-by-character, and
``'1' < '2'``.  With `sbpy`'s ``natural_sort_key``, numerical
comparisons are made whenever possible:

    >>> from sbpy.data import natural_sort_key
    >>> sorted(comets, key=natural_sort_key)
    ['2P/Encke', '9P/Tempel 1', '10P/Tempel 2', '101P/Chernykh']

Packed Numbers and Designations
-------------------------------

`~sbpy.data.Names.from_packed` and `~sbpy.data.Names.to_packed`
provide functionality to convert between packed designations and
numbers and unpacked ones:

    >>> Names.from_packed('J95A01A')
    '1995 AA1'
    >>> Names.from_packed('G3693')
    163693
    >>> Names.to_packed('1995 AA1')
    'J95A01A'
    >>> Names.to_packed('163693')
    'G3693'
