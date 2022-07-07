===========================================
Small-Body Name Parsing (`sbpy.data.Names`)
===========================================

`~sbpy.data.Names` is different from the other classes in `~sbpy.data`
in that it does not use `~sbpy.data.DataClass` as a base class. Instead,
`~sbpy.data.Names` does not contain any data, it merely serves as an
umbrella for functions to identify asteroid and comet names, numbers,
and designations.


Cometary and Asteroidal Name Parsing
------------------------------------

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

    >>> Names.parse_asteroid('(228195) 6675 P-L')
    {'number': 228195, 'desig': '6675 P-L'}
    >>> Names.parse_comet('73P-C/Schwassmann Wachmann 3 C')
    {'type': 'P', 'number': 73, 'fragment': 'C', 'name': 'Schwassmann Wachmann 3 C'}

These methods will raise exceptions when the name cannot be parsed as expected:

    >>> Names.parse_asteroid('C/2001 A2-A (LINEAR)')
    Traceback (most recent call last):
    ...
    sbpy.data.names.TargetNameParseError: C/2001 A2-A (LINEAR) does not appear to be an asteroid identifier
    >>> Names.parse_comet('12893')
    Traceback (most recent call last):
    ...
    sbpy.data.names.TargetNameParseError: 12893 does not appear to be a comet name

In order to be able to distinguish between asteroid and comet identifiers,
`sbpy` follows the MPC guideline in that it requires comet identifiers to
include the comet type in either in combination with a number (e.g.,
``'259P'``), a name (e.g., ``'P/Halley'``), or both (e.g., ``'2P/Encke'``). For
instance, the identifier ``'Halley'`` would be identified as an asteroid, as it
lacks a comet type identifier. Hence, some caution is advised when using these
routines - identification might not be unambiguous.


A/ objects: asteroids in cometary orbits
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Small bodies designated with an A/ prefix have cometary orbits, but appear
asteroidal [MPEC2018H54]_.  ``sbpy`` considers them to be asteroids:

    >>> Names.asteroid_or_comet('A/2018 V3')
    'asteroid'


Interstellar objects
^^^^^^^^^^^^^^^^^^^^

Interstellar object designations, which start with an I/, do not
give any insight into the nature of the object.  For example, 1I/ʻOumuamua was
asteroidal in appearance but 2I/Borisov was cometary.  ``sbpy`` raises an
exception for I/ objects:

    >>> Names.asteroid_or_comet('1I/ʻOumuamua')
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/disks/data0/astro/Projects/sbpy/sbpy/data/names.py", line 597, in asteroid_or_comet
        raise TargetNameParseError('Target nature unclear.')
    sbpy.data.names.TargetNameParseError: Target nature unclear.


.. [MPEC2018H54] Williams, G. V. 2018.  A/ Objects.  MPEC 2018-H54.  https://minorplanetcenter.net/mpec/K18/K18H54.html


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
