.. _data containers:

=======================================
Data Containers (`sbpy.data.DataClass`)
=======================================

`sbpy` relies heavily on the use of `~sbpy.data.DataClass` data
containers that are used to encapsulate data and to propagate them
through your workflow.

`sbpy` provides data containers for orbital elements
(`~sbpy.data.Orbit`), ephemerides (`~sbpy.data.Ephem`), observational
data (`~sbpy.data.Obs`), and physical properties
(`~sbpy.data.Phys`) data. 



What are Data Containers?
=========================

`~sbpy.data.DataClass` - and hence all the data containers presented
here - uses a `~astropy.table.QTable` object under the hood. You can
think of those as **tables** - consisting of **fields** (or columns)
and **rows** - that have `~astropy.units` attached to them, allowing
you to propagate these units through your programs. **We strongly urge
the user to make use of** `~astropy.units` in the definition of data
containers in order to minimize confusion and tap the full potential
of `sbpy`. Finally, `~sbpy.data.DataClass` objects have a ``meta``
attribute that enables the user to label these objects with
unstructured meta data.

The user is free to add any fields they want to a
`~sbpy.data.DataClass` object. However, in order to enable the
seamless use of `sbpy` functions and methods, we require the user to
pick among a few common field names for different properties as listed
:ref:`here <field name list>`. `~sbpy.data.DataClass` objects
are able to identify alternative field names as listed in this
document, as well as to perform transformations between a few field
names - see below for more details.

A `~sbpy.data.DataClass` object can hold as many data rows as you
want. All rows can refer to a single object, or each row can refer to
a separate object - this is usually up to the user, restrictions exist
only in a few cases as detailed in this documentation.


Why are there different Data Container Types?
=============================================

Although the data container classes `~sbpy.data.Ephem`,
`~sbpy.data.Orbit`, `~sbpy.data.Obs` and `~sbpy.data.Phys` look very
similar and are based on the same base class (`~sbpy.data.DataClass`),
there are subtle differences. The container classes have been
introduced to minimize confusion between properties of different
natures, i.e., to avoid mixing apples with oranges. For instance,
`~sbpy.data.Ephem` objects deal with ephemerides (e.g., RA and DEC and
other properties that change as a function of time); the diameter of
an object, however, does usually not vary with time. Adding a
"diameter" field to an `~sbpy.data.Ephem` object would hence be
inefficient: this object can contain hundreds of rows describing a
target's ephemerides, while the diameter would be the same in each of
these rows. Furthermore, providing separate data containers for
different properties enables the implementation of specifically
designed methods to query and modify the data held by these classes.


Data Container Type Overview
============================

Ephem
-----

`~sbpy.data.Ephem` has been designed to hold
**ephemerides**, i.e., properties that vary with time. 

`~sbpy.data.Ephem` currently provides convenience functions to query
ephemerides from the JPL Horizons system
(`~sbpy.data.Ephem.from_horizons`) the Minor Planet Center
(`~sbpy.data.Ephem.from_mpc`), IMCCE's Miriade system
(`~sbpy.data.Ephem.from_miriade`) as well as a convenience function to
derive ephemerides from an `~sbpy.data.Orbit` object using `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_.

Obs
---

`~sbpy.data.Obs` is tailored to holding **observational data**, e.g.,
magnitudes as a function of time. The `~sbpy.data.Obs` class is the
only data container that is not directly derived from
`~sbpy.data.DataClass`, but from `~sbpy.data.Ephem`, providing the
same functionality as the latter. `~sbpy.data.Obs.from_mpc` enables
you to query observations reported to the MPC for a specific target
and `~sbpy.data.Obs.supplement` queries one of the ephemerides service
to supplement your observation data.


Orbit
-----

`~sbpy.data.Orbit` should be used to hold **orbital elements** of one
or several bodies. Elements can be retrieved using the convenience
function `~sbpy.data.Orbit.from_horizons`, propagated using
`~sbpy.data.Orbit.oo_propagate`, and transformed into other frames
using `~sbpy.data.Orbit.oo_transform`.

Phys
----

`~sbpy.data.Phys` objects are meant to hold **physical properties**
that do not change over time. Known physical properties can currently
be queried from the JPL Small-Body Database Browser system using
`~sbpy.data.Phys.from_sbdb`.


Names
-----

`~sbpy.data.Names` objects are somewhat different from the other data
containers, as they don't hold properties but only object
**names**. These names can be used to identify object nature
(`~sbpy.data.Names.asteroid_or_comet`) and they can be parsed to
extract individual identifier components
(`~sbpy.data.Names.parse_asteroid` and
`~sbpy.data.Names.parse_comet`).

.. _How to use Data Containers:

How to use Data Containers
==========================

All of the data objects dealt with in `sbpy.data` share the same
common base class: `sbpy.data.DataClass`. `~sbpy.data.DataClass`
defines the basic functionality and makes sure that all `sbpy.data`
objects can used in the exact same way.

In plain words, this means that in the following examples you can
replace `~sbpy.data.DataClass`, `~sbpy.data.Ephem`,
`~sbpy.data.Orbit`, `~sbpy.data.Obs`, and `~sbpy.data.Phys` object
with each other. In order to show some useful use cases, we will
iterate between these types, but keep in mind: they all work the exact
same way.


Building Data Containers
------------------------

While `~sbpy.data.Ephem`, `~sbpy.data.Orbit`, `~sbpy.data.Obs`, and
`~sbpy.data.Phys` provide a range of convenience functions to build
objects containing data, for instance from online data archives, it is
easily possible to build these objects from scratch. This can be done
for input data stored in dictionaries
(`~sbpy.data.DataClass.from_dict`), lists or arrays
(`~sbpy.data.DataClass.from_columns` and
`~sbpy.data.DataClass.from_rows`), `~astropy.table.Table` objects
(`~sbpy.data.DataClass.from_table`), or from data files
(`~sbpy.data.DataClass.from_file`).

Depending on how your input data are organized, you can use different
options in different cases:

Building a Data Container from a Dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assume that you want to build an `~sbpy.data.Orbit` object to
propagate this orbit and obtain ephemerides. Since you are dealing
with a single orbit, the most convenient solution might be to use a
dictionary to build your object:

    >>> from sbpy.data import Orbit
    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> elements = {'a':1.234*u.au, 'e':0.1234, 'i':12.34*u.deg,
    ...             'argper': 123.4*u.deg, 'node': 45.2*u.deg,
    ...             'epoch': Time(2451200.5, format='jd'), 'true_anom':23.1*u.deg}
    >>> orb = Orbit.from_dict(elements)
    >>> orb
    <QTable length=1>
       a       e       i     argper   node    epoch   true_anom
       AU             deg     deg     deg                deg
    float64 float64 float64 float64 float64    Time    float64
    ------- ------- ------- ------- ------- --------- ---------
      1.234  0.1234   12.34   123.4    45.2 2451200.5      23.1

One quick note on building `~sbpy.data.DataClass` objects from
dictionaries: dictionaries have no intrinsic order. In dictionary
``elements`` as defined here, there is no guarantee that ``'a'`` will
always be located before ``'e'`` when reading out the dictionary item
by item, which happens when the data table is built in the
background. Hence, the order of the resulting data table columns has
to be considered random. If you want to force a specific order on the
columns in your data table, you can use an `~collections.OrderedDict`
instead of a simple dictionary. The order of elements in an
`~collections.OrderedDict` will be the same as the order of the data
table columns.

For details on how to build objects from dictionaries, see
`~sbpy.data.DataClass.from_dict`.

Building a Data Container from Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now assume that you want to build an `~sbpy.data.Obs` object holding
RA, Dec, and observation midtime for some target that you observed. In
this case, you can use `~sbpy.data.DataClass.from_columns` as shown
here:

    >>> from sbpy.data import Obs
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> from numpy import array
    >>> ra = [10.223423, 10.233453, 10.243452]*u.deg
    >>> dec = [-12.42123, -12.41562, -12.40435]*u.deg
    >>> epoch = Time(2451523.5 + array([0.1234, 0.2345, 0.3525]), format='jd')
    >>> obs = Obs.from_columns([ra, dec, epoch], names=['ra', 'dec', 't'])
    >>> obs
    <QTable length=3>
        ra       dec         t      
       deg       deg                
     float64   float64      Time    
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234
    10.233453 -12.41562 2451523.7345
    10.243452 -12.40435 2451523.8525

Note how ``epoch`` is handled differently: it is provided to
``Obs.from_column`` as a `~astropy.time.Time` object (see
:ref:`user_zen` for a discussion). 
    
For details on how to build objects from lists or arrays, see
`~sbpy.data.DataClass.from_columns` and also
`~sbpy.data.DataClass.from_rows`, depending on whether your data is
represented as rows or columns. Note that you could also use
`~sbpy.data.DataClass.from_dict` by providing column data to the
different fields.

Building a Data Container from a Table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your data are already available as a `~astropy.table.Table` or
`~astropy.table.QTable`, you can simply convert it into a
`~sbpy.data.DataClass` object using `~sbpy.data.DataClass.from_table`.

Building a Data Container from a File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also read in the data from a file that should be properly
formatted using `~sbpy.data.DataClass.from_file`. This function merely
serves as a wrapper for `astropy.table.Table.read` and uses the same
parameters as the latter function; please refer to `this document
<https://docs.astropy.org/en/stable/table/io.html>`_ for a review.

As an example, you can read in a properly formatted ASCII file using
the following lines:

   >>> from sbpy.data import Ephem
   >>> data = Ephem.from_file('data.txt', format='ascii') # doctest: +SKIP

Please note that the file formats available (see `here
<https://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers>`_
for a list of available formats) provide varying support for units and
meta data. For instance, ``basic``, ``csv``, ``html``, and ``latex``
do not provide unit or meta data information. However, ``fits``,
``cds``, ``daophot``, ``ecsv``, and ``ipac`` do support units and meta
data.


Building a Data Container from an Online Query
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most `~sbpy.data.DataClass` data containers offer convenience
functions to query data from online service. Please refer to the
corresponding classes for information and examples for querying data.


A Note on Field Names
---------------------

In order for `sbpy` to properly identify the fields that might be
necessary for calculations, default column names should be used to
name these fields. For instance, a column of Right Ascensions should
be named ``'RA'`` or ``'ra'``. For a list of acceptable field names,
please refer to the list of :ref:`field name list`.

Also note that `sbpy` is able to use :ref:`alternative field names
<fieldnames>`, but only those that are listed in the list of
:ref:`field name list`.


Accessing data
--------------

In order to obtain a list of field names in a `~sbpy.data.DataClass`
object, you can use `~sbpy.data.DataClass.field_names`:

.. testsetup::

    >>> # in case requirements for tests above are not met
    >>> from sbpy.data import Obs
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> from numpy import array
    >>> ra = [10.223423, 10.233453, 10.243452]*u.deg
    >>> dec = [-12.42123, -12.41562, -12.40435]*u.deg
    >>> epoch = Time(2451523.5 + array([0.1234, 0.2345, 0.3525]), format='jd')
    >>> obs = Obs.from_columns([ra, dec, epoch], names=['ra', 'dec', 't'])

.. doctest::

    >>> obs.field_names
    ['ra', 'dec', 't']

You can also use the `in` operator to check if a field is contained in
a `~sbpy.data.DataClass` object.  Alternative field names can be used
for the `in` test:

    >>> 'ra' in obs
    True
    >>> 'RA' in obs
    True

Each of these columns can be accessed easily, for instance:

    >>> obs['ra']
    <Quantity [10.223423, 10.233453, 10.243452] deg>

which will return an `~astropy.units.quantity.Quantity` object if that
column has a `~astropy.units.Unit` attached to it or a `~astropy.table.Column`
otherwise.

Similarly, if you are interested in the first set of observations in
``obs``, you can use:

    >>> obs[0]
    <QTable length=1>
        ra       dec         t
       deg       deg
     float64   float64      Time
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234

which returns you a new instance of the same class as your original
objet with only the requested subset of the
data. In order to retrieve RA from the second observation, you can
combine both examples and do:

    >>> obs[1]['ra']
    <Quantity [10.233453] deg>


Just like in any `~astropy.table.Table` or `~astropy.table.QTable`
object, you can use slicing to obtain subset tables from your data,
for instance:

    >>> obs['ra', 'dec']
    <QTable length=3>
        ra       dec
       deg       deg
     float64   float64
    --------- ---------
    10.223423 -12.42123
    10.233453 -12.41562
    10.243452 -12.40435
    <BLANKLINE>
    >>> obs[:2]
    <QTable length=2>
        ra       dec         t
       deg       deg
     float64   float64     Time
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234
    10.233453 -12.41562 2451523.7345
    <BLANKLINE>
    >>> obs[obs['ra'] <= 10.233453 * u.deg]
    <QTable length=2>
        ra       dec         t
       deg       deg
     float64   float64     Time
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234
    10.233453 -12.41562 2451523.7345

The results of these examples will be of the same data type as ``obs``
(or really just any type derived from `~sbpy.data.DataClass`, e.g.,
`~sbpy.data.Ephem`, `~sbpy.data.Orbit`, ...)  The latter example shown
here uses a condition to filter data (only those observations with RA
less than or equal to 10.233453 degrees; note that it is necessary
here to apply ``u.deg`` to the value that all the RAs are compared
against) but selects all the columns in the original table.

If you ever need to access the actual `~astropy.table.QTable` object
that is inside each `~sbpy.data.DataClass` object, you can access it
as ``obs.table``.

Modifying an object
-------------------

Individual elements, entire rows, and columns can be modified by
directly addressing them:

    >>> obs['ra']
    <Quantity [10.223423, 10.233453, 10.243452] deg>
    >>> obs['ra'] = obs['ra'] + 0.1 * u.deg
    >>> obs['ra']
    <Quantity [10.323423, 10.333453, 10.343452] deg>

The basic functionalities to modify the data table are implemented in
`~sbpy.data.DataClass`, including adding rows and columns and stack a
DataClass with another DataClass object or an `~astropy.table.Table`
object.

Let's assume you want to add some more observations to your ``obs``
object:

    >>> obs.add_row([10.255460 * u.deg, -12.39460 * u.deg, 2451523.94653 * u.d])
    >>> obs
    <QTable length=4>
        ra       dec          t      
       deg       deg      
     float64   float64      Time
    --------- --------- -------------
    10.323423 -12.42123  2451523.6234
    10.333453 -12.41562  2451523.7345
    10.343452 -12.40435  2451523.8525
     10.25546  -12.3946 2451523.94653
  

or if you want to add a column to your object:

    >>> obs.apply(['V', 'V', 'R', 'i'], name='filter')
    >>> obs
    <QTable length=4>
        ra       dec          t       filter
       deg       deg                        
     float64   float64       Time     str32
    --------- --------- ------------- ------
    10.323423 -12.42123  2451523.6234      V
    10.333453 -12.41562  2451523.7345      V
    10.343452 -12.40435  2451523.8525      R
     10.25546  -12.3946 2451523.94653      i

The same result can be achieved using the following syntax:

    >>> obs['filter2'] = ['V', 'V', 'R', 'i']
    >>> obs
    <QTable length=4>
        ra       dec          t       filter filter2
       deg       deg                                
     float64   float64       Time     str32    str1
    --------- --------- ------------- ------ -------
    10.323423 -12.42123  2451523.6234      V       V
    10.333453 -12.41562  2451523.7345      V       V
    10.343452 -12.40435  2451523.8525      R       R
     10.25546  -12.3946 2451523.94653      i       i

Similarly, existing columns can be modified using:

    >>> obs['filter'] = ['g', 'i', 'R', 'V']

If you want to stack two observations into a single object:

    >>> ra = [20.223423, 20.233453, 20.243452] * u.deg
    >>> dec = [12.42123, 12.41562, 12.40435] * u.deg
    >>> phase = [10.1, 12.3, 15.6] * u.deg
    >>> epoch = Time(2451623.5 + array([0.1234, 0.2345, 0.3525]), format='jd')
    >>> obs2 = Obs.from_columns([ra, dec, epoch, phase],
    ...     names=['ra', 'dec', 't', 'phase'])
    >>>
    >>> obs.vstack(obs2)
    >>> obs
    <QTable length=7>
        ra       dec          t       filter filter2  phase
       deg       deg                                   deg
     float64   float64       Time      str1    str1  float64
    --------- --------- ------------- ------ ------- -------
    10.323423 -12.42123  2451523.6234      g       V     ———
    10.333453 -12.41562  2451523.7345      i       V     ———
    10.343452 -12.40435  2451523.8525      R       R     ———
     10.25546  -12.3946 2451523.94653      V       i     ———
    20.223423  12.42123  2451623.6234     --      --    10.1
    20.233453  12.41562  2451623.7345     --      --    12.3
    20.243452  12.40435  2451623.8525     --      --    15.6

Note that the data table to be stacked doesn't have to have the same
columns as the original data table.  A keyword `join_type` is used to
decide how to process the different sets of columns.  See
`~astropy.table.Table.vstack()` for more detail.

Because the underlying `~astropy.table.QTable` can be exposed by the
`~sbpy.data.DataClass.table` property, it is possible to modify the data
table by directly accessing the underlying `~astropy.table.QTable` object.
However, this is not generally advised.  You should use the mechanisms provided
by `~sbpy.data.DataClass` to manipulate the data table as much as possible
to maintain the integrity of the data table.

Additional Data Container Concepts
==================================

.. _fieldnames:

Alternative field names
-----------------------

If you ask 3 different planetary astronomers which field name or
variable name they use for the orbital inclination, you will receive 3
different answers. Good candidates might be ``'i'``, ``'inc'``, or
``'incl'`` - it's a matter of personal taste. The `sbpy` developers
are aware of this ambiguity and hence `~sbpy.data.DataClass` provides
some flexibility in the use of field name. This functionality is
established through a list of acceptable field names that are
recognized by `sbpy`, which is provided in the
:ref:`field name list`.

As an example, if your `~sbpy.data.Orbit` object has a column named
``'incl'`` but you try to get column ``'i'``, the object will
internally check if ``'i'`` is a legitimate field name and what its
alternatives are, and it will find that a field name ``'incl'`` exists
in the object. The corresponding ``'incl'`` column is then
returned. If you try to get a field name that is not connected to any
existing field name, a ``KeyError`` will be raised.

    >>> from sbpy.data import Orbit
    >>> orb = Orbit.from_dict({'incl': [1, 2, 3]*u.deg})
    >>> orb['i']
    <Quantity [1., 2., 3.] deg>

The definition of alternative field names is done in the file
``sbpy/data/__init__.py``, using the list ``fieldnames``. This list is
automatically tested for potential naming conflicts, i.e., different
properties that share the same alternative field names, and a
human-readable list is compiled upon building `sbpy`.

The full list of field names is available here:
:ref:`field name list`.

Field conversions
-----------------

There are parameters and properties that can be used synonymously, a
good example for which are an object's radius and diameter. `sbpy`
acknowledges identities like this by providing internal conversions
for such properties. Consider the following example:

    >>> from sbpy.data import Phys
    >>> import astropy.units as u
    >>> data = Phys.from_dict({'d': 10*u.km})
    >>> print('{:.1f}'.format(data['d'][0]))
    10.0 km
    >>> print('{:.1f}'.format(data['radius'][0]))
    5.0 km

Note that the radius is not explicitly defined in ``data``, but
derived internally upon querying it and added to the internal data table:

    >>> print(data.field_names)
    ['d', 'radius']
    

.. _epochs:
    
Epochs and the use of astropy.time
----------------------------------

Epochs and data referring to specific points in time have to be
provided as `~astropy.time.Time` objects. The advantage of these
objects is their flexibility in terms of format and time
scale. `~astropy.time.Time` objects can be readily transformed into a
wide range of formats; for instance, ``Time('2019-07-23 10:49').jd``
can be used to convert an ISO epoch to a Julian date.

More importantly, `~astropy.time.Time` provides functionality to
transform epochs between different time scales. Hence, every
`~astropy.time.Time` object comes with a time scale (UTC is used
by default) and can be easily transformed into a different time
scale. The following example defines an epoch in UTC and as a Julian
date and transforms it to TDB:

    >>> from astropy.time import Time
    >>> epoch = Time(2451200, format='jd')
    >>> epoch
    <Time object: scale='utc' format='jd' value=2451200.0>
    >>> epoch.tdb
    <Time object: scale='tdb' format='jd' value=2451200.000742876>
    >>> epoch.tdb.iso
    '1999-01-21 12:01:04.184'

Using `~astropy.time.Time` in `~sbpy.data.DataClass` objects is
straightforward. The following example builds a simple
`~sbpy.data.Obs` object from a dictionary:

    >>> from sbpy.data import Obs
    >>> times = ['2018-10-01', '2018-11-01', '2018-12-01']
    >>> obs = Obs.from_dict({'epoch': Time(times), 'mag': [10, 12, 14]*u.mag})
    >>> obs
    <QTable length=3>
             epoch            mag  
                              mag  
              Time          float64
    ----------------------- -------
    2018-10-01 00:00:00.000    10.0
    2018-11-01 00:00:00.000    12.0
    2018-12-01 00:00:00.000    14.0

The ``'epoch'`` column in ``obs`` can be used like any other field or
`~astropy.time.Time` object. The following example converts the epoch
to TAI and Julian date:

    >>> obs['epoch'].tai.jd
    array([2458392.50042824, 2458423.50042824, 2458453.50042824])

Note that different functions in `sbpy` have different requirements on
the time scale of `~astropy.time.Time` objects. Fortunately,
`~astropy.time.Time` objects are able to convert most time scales
seamlessly. However, that requires that some user-defined time scale
might have to be converted to other time scale for compatibility
reasons internally, which also means that outpu t epochs usually
follow this forced time scale. In order to notify the user that the
time scale has been changed, a `~sbpy.data.TimeScaleWarning` will be
issued.


Writing object data to a file
-----------------------------

`~sbpy.data.DataClass` objects can be written to files using
`~sbpy.data.DataClass.to_file`:

.. testsetup::

    >>> import os
    >>> assert not os.path.exists('observations.dat')

.. doctest::

    >>> obs.to_file('observations.dat')

.. testcleanup::

    >>> os.unlink('observations.dat')

By default, the data are written in ASCII format, but other formats
are available, too (`list of file formats
<https://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers>`_). Please
note that not all file formats support units and meta data. For
instance, ``basic``, ``csv``, ``html``, and ``latex`` do not provide
unit or meta data information. However, ``fits``, ``cds``,
``daophot``, ``ecsv``, and ``ipac`` do support units and meta data.
