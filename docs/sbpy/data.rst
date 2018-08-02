Data Module (`sbpy.data`)
=========================

Introduction
------------

`sbpy.data` provides classes for dealing with orbital elements
(`~sbpy.data.Orbit`), ephemerides (`~sbpy.data.Ephem`), and physical
properties (`~sbpy.data.Phys`). `~sbpy.data.Ephem`,
`~sbpy.data.Orbit`, and `~sbpy.data.Phys` objects act as containers
for such parameters and can (and should) be used to provide these to
functions in `sbpy`. Each of these classes is based on the
`~sbpy.data.DataClass` base class, which internally uses an
`~astropy.table.QTable` object and provides the same functionality and
features as the latter.

Furthermore, `~sbpy.data` also provides additional interfaces to a number of
different services and `~sbpy.data.Names` provides functions
related to naming conventions for asteroids and comets.


How to use Ephem, Orbit, and Phys objects
-----------------------------------------

All of the data objects dealt with in `sbpy.data` share the same
common base class: `sbpy.data.DataClass`. `~sbpy.data.DataClass`
defines the basic functionality and makes sure that all `sbpy.data`
objects can used in the exact same way.

In plain words, this means that in the following examples you can
replace `~sbpy.data.DataClass`, `~sbpy.data.Ephem`,
`~sbpy.data.Orbit`, and `~sbpy.data.Phys` object with each other. In
order to show some useful use cases, we will iterate between these
types, but keep in mind: they all work the exact same way.

`~sbpy.data.DataClass` uses `~astropy.table.QTable` objects under the
hood. You can think of those as tables - consisting of columns and
rows - that have `~astropy.units` attached to them, allowing you to
propagate these units through your code. Each `~sbpy.data.DataClass`
object can hold as many data as you want, where each datum can be a
different object or the same object at a different epoch.


Building an object
^^^^^^^^^^^^^^^^^^

While `~sbpy.data.Ephem`, `~sbpy.data.Orbit`, and `~sbpy.data.Phys`
provide a range of convience functions to build objects containing
data, for instance from online data archives, it is easily possible to
build these objects from scratch. This can be done for input data
stored in dictionaries (`~sbpy.data.DataClass.from_dict`), lists or
arrays (`~sbpy.data.DataClass.from_array`), `~astropy.table.Table`
objects (`~sbpy.data.DataClass.from_table`), or from data files
(`~sbpy.data.DataClass.from_file`).

Depending on how your input data are organized, you can use different
options in different cases:

1. Assume that you want to build a `~sbpy.data.Orbit` object to
   propagate this orbit and obtain ephemerides. Since you are dealing
   with a single orbit, the most convenient solution might be to use a
   dictionary to build your object:

    >>> from sbpy.data import Orbit
    >>> import astropy.units as u
    >>> elements = {'a':1.234*u.au, 'e':0.1234, 'i':12.34*u.deg,
    ...             'argper': 123.4*u.deg, 'node': 45.2*u.deg,
    ...             'epoch': 2451200.5*u.d, 'true_anom':23.1*u.deg}
    >>> orb = Orbit.from_dict(elements)
    >>> print(orb)  # doctest:+ELLIPSIS
    <sbpy.data.orbit.Orbit object at ...>

2. Now assume that you want to build an `~sbpy.data.Ephem` object
   holding RA, Dec, and observation midtime for some target that you
   observed. In this case, you could provide a list of three
   dictionaries to `~sbpy.data.DataClass.from_dict`, which means a lot
   of typing. Instead, you can use `~sbpy.data.DataClass.from_array`,
   which allows to provide your input data in the form of a list,
   tuple, or `~numpy.ndarray`:

    >>> from sbpy.data import Ephem
    >>> import astropy.units as u
    >>> from numpy import array
    >>> ra = [10.223423, 10.233453, 10.243452]*u.deg
    >>> dec = [-12.42123, -12.41562, -12.40435]*u.deg
    >>> epoch = (2451523.5 + array([0.1234, 0.2345, 0.3525]))*u.d
    >>> obs = Ephem.from_array([ra, dec, epoch], names=['ra', 'dec', 't'])
    >>> print(obs)  # doctest:+ELLIPSIS   
    <sbpy.data.ephem.Ephem object at ...>

3. If your data are already available as a `~astropy.table.Table` or
   `~astropy.table.QTable`, you can simply convert it into a
   `~sbpy.data.DataClass` object using
   `~sbpy.data.DataClass.from_table`.

4. You can also read in the data from a file that should be properly
   formatted (e.g., it should have a headline with the same number of
   elements as there are columns) using
   `~sbpy.data.DataClass.from_file`. This function merely serves as a
   wrapper for `~astropy.table.Table.read` and uses the same
   parameters as the latter function. You can read in an ASCII file
   using the following lines:

   >>> from sbpy.data import Ephem
   >>> data = Ephem.from_file('data.txt', format='ascii') # doctest: +SKIP

   Please not that `~sbpy.data.DataClass.from_file` is not able to
   identify units automatically. If you want to take advantage for
   `~astropy.units` you will have to assign these units manually later
   on.


Accessing an object
^^^^^^^^^^^^^^^^^^^

In order to obtain a list of column names in a `~sbpy.data.DataClass` object, you can use `~sbpy.data.DataClass.column_names`:

    >>> obs.column_names
    <TableColumns names=('ra','dec','t')>

Each of these columns can be accessed easily, for instance:

    >>> obs['ra']
    [10.223423 10.233453 10.243452] deg

Similarly, if you are interested in the first set of observations in
``obs``, you can use:

    >>> obs[0]
        ra       dec         t      
       deg       deg         d      
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234

which returns you a table with only the requested subset of the
data. In order to retrieve RA from the second observation, you can
combine both examples and do:

    >>> obs[1]['ra']
    10.233453 deg

Just like in any `~astropy.table.Table` or `~astropy.table.QTable` object, you can use slicing to obtain subset tables from your data, for instance:

    >>> obs['ra', 'deg']
        ra       dec   
       deg       deg   
    --------- ---------
    10.223423 -12.42123
    10.233453 -12.41562
    10.243452 -12.40435

    >>> obs[obs['ra'] <= 10.233453*u.deg]
        ra       dec         t      
       deg       deg         d      
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234
    10.233453 -12.41562 2451523.7345

The latter uses a condition to filter data (only those observations
with RA less than or equal to 10.233453 degrees; note that it is
necessary here to apply ``u.deg`` to the value that all the RAs are
compared against) but selects all the columns in the original table.

If you ever need to access the actual `~astropy.table.QTable` object
that is inside each `~sbpy.data.DataClass` object, you can access it
as ``obs.table``.
    
Modifying an object
^^^^^^^^^^^^^^^^^^^

`~sbpy.data.DataClass` offers some convenience functions for object
modifications. It is trivial to add additional rows and columns to
these objects in the form of lists, arrays, or dictionaries. 

Let's assume you want to add some more observations to your ``obs``
object:

    >>> obs.add_rows([[10.255460*u.deg, -12.39460*u.deg, 2451523.94653*u.d],
    ...               [10.265425*u.deg, -12.38246*u.deg, 2451524.0673*u.d]])
    5
    >>> obs.table
        ra       dec          t      
       deg       deg          d      
    --------- --------- -------------
    10.223423 -12.42123  2451523.6234
    10.233453 -12.41562  2451523.7345
    10.243452 -12.40435  2451523.8525
     10.25546  -12.3946 2451523.94653
    10.265425 -12.38246  2451524.0673

or if you want to add a column to your object:

    >>> obs.add_column(['V', 'V', 'R', 'i', 'g'], name='filter')
    4
    >>> obs.table
        ra       dec          t       filter
       deg       deg          d             
    --------- --------- ------------- ------
    10.223423 -12.42123  2451523.6234      V
    10.233453 -12.41562  2451523.7345      V
    10.243452 -12.40435  2451523.8525      R
     10.25546  -12.3946 2451523.94653      i
    10.265425 -12.38246  2451524.0673      g
    
A few things to be mentioned here:

* Note how both functions return the number of rows or columns in the
  updated object.
* If you are adding rows, the elements in the rows will be assigned to
  the column in the corresponding order of the table columns. The
  `~astropy.units` of the row elements have to be of the same
  dimension as the table columns (e.g., one of the table column units
  is degrees, then the corresponding row element has to define an
  angular distance: ``u.deg`` or ``u.rad``).
* Naturally, the number of columns and rows of the rows and columns
  to be added has to be identical to the numbers in the data table.

If you are trying to add a single row to your object data table, using a dictionary might be the most convenient solution:

    >>> obs.add_rows({'ra':10.255460*u.deg, 'dec': -12.39460*u.deg,
    ...               't': 2451524.14653*u.d, 'filter': 'z'})
    6

    
When adding a large number of rows to your object, it might be most
convenient to first convert all the new rows into new
`~sbpy.data.DataClass` object and then append that using
`~sbpy.data.DataClass.add_rows`:

    >>> obs2 = Ephem.from_array([[10.4545, 10.5656]*u.deg,
    ...                          [-12.1212, -12.0434]*u.deg,
    ...                          [2451524.14653, 2451524.23541]*u.d,
    ...                          ['r', 'z']],
    ...                         names=['ra', 'dec', 't', 'filter'])
    >>> obs.add_rows(obs2)
    7

Individual elements, entire rows, and columns can be modified by
directly addressing them:

    >>> print(obs['ra'])
    [10.223423 10.233453 10.243452 10.25546  10.265425 10.4545   10.5656  ] deg
    >>> obs['ra'][:] = obs['ra'] + 0.1*u.deg
    >>> print(obs['ra'])
    [10.323423 10.333453 10.343452 10.35546  10.365425 10.5545   10.6656  ] deg

Note the specific syntax in this case (``obs['ra'][:] = ...``) that
is required by `~astropy.table.Table` if you want to replace
an entire column.
    
More complex data table modifications are possible by directly
accessing the underlying `~astropy.table.QTable` object.

Writing object data to a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`~sbpy.data.DataClass` objects can be written to files using
`~sbpy.data.DataClass.to_file`:

    >>> obs.to_file('observations.dat')

By default, the data are written in ASCII format, but other formats
are available, too (cf. `~astropy.table.Table.write`).

  
How to use Orbit
----------------
tbd


How to use Ephem
----------------
tbd


How to use Phys
---------------
tbd


How to use Names
----------------

`~sbpy.data.Names` is different from the other classes in `~sbpy.data`
in that it does not use `~sbpy.data.DataClass` as a base class. Instead,
`~sbpy.data.Names` does not contain any data, it merely serves as an
umbrella for functions to identify asteroid and comet names, numbers,
and designations.

In order to distinguish if a string designates a comet or an asteroid,
you can use the following code:

    >>> from sbpy.data import Names
    >>> print(Names.asteroid_or_comet('(1) Ceres'))
    asteroid
    >>> print(Names.asteroid_or_comet('2P/Encke'))
    comet

The module basically uses regular expressions to match the input
strings and find patterns that agree with asteroid and comet names,
numbers, and designations. There are separate tasks to identify
asteroid and comet identifiers:

    >>> print(Names.parse_asteroid('(228195) 6675 P-L'))
    {'number': 228195, 'desig': '6675 P-L'}
    >>> print(Names.parse_asteroid('C/2001 A2-A (LINEAR)')) # doctest: _ELLIPSIS
    ...sbpy.data.names.TargetNameParseError: C/2001 A2-A (LINEAR) does not appear to be an asteroid identifier
    >>> print(Names.parse_comet('12893')) # doctest: +ELLIPSIS
    ...sbpy.data.names.TargetNameParseError: 12893 does not appear to be a comet name
    >>> print(Names.parse_comet('73P-C/Schwassmann Wachmann 3 C	'))
    {'type': 'P', 'number': 73, 'fragment': 'C', 'name': 'Schwassmann Wachmann 3 C'}
    
Note that these examples are somewhat idealized. Consider the
following query:

    >>> print(Names.parse_comet('12893 Mommert (1998 QS55)'))
    {'name': 'Mommert ', 'desig': '1998 QS55'}

Although this target identifier clearly denotes an asteroid, the
routine finds a comet name and a comet designation. The reason for
this is that some comets are discovered as asteroids and hence obtain
asteroid-like designations that stick to them; similarly, comet names
cannot be easily distinguished from asteroids names, unless one knows
all comet and asteroid names. Hence, some caution is advised when
using these routines - identification might not be unambiguous.

    


Reference/API
-------------
.. automodapi:: sbpy.data
    :no-heading:
		


 
