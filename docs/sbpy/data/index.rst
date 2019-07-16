Data Module (`sbpy.data`)
=========================

Introduction
------------

`sbpy.data` provides classes for dealing with all kinds of data that
planetary astronomers encounter: orbital elements
(`~sbpy.data.Orbit`), ephemerides (`~sbpy.data.Ephem`), observational
data (`~sbpy.data.Obs`), and physical properties
(`~sbpy.data.Phys`). These objects act as containers that bundle
similar data and have to be used in `sbpy` to maintain a transparent
and clean data flow between functions. Each of these classes is based
on the `~sbpy.data.DataClass` base class, which internally uses an
`~astropy.table.QTable` object in combination with some other
features. Hence, each of these objects can actually be thought of as a
table with different *fields* (or columns) and *rows*.

However, in contrast to simple tables, `~sbpy.data` objects provide
additional functionality that is detailed below.


What are Ephem, Orbit, Obs, and Phys and what are they used for?
----------------------------------------------------------------

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

What are Data Containers?
^^^^^^^^^^^^^^^^^^^^^^^^^

`~sbpy.data.DataClass` - and hence all the data containers presented
here - uses a `~astropy.table.QTable` object under the hood. You can
think of those as **tables** - consisting of **fields** (or columns)
and **rows** - that have `~astropy.units` attached to them, allowing
you to propagate these units through your programs. **We strongly urge
the user to make use of** `~astropy.units` in the definition of data
containers in order to minimize confusion and tap the full potential
of `sbpy`.

The user is free to add any fields they want to a
`~sbpy.data.DataClass` object. However, in order to enable the
seemless use of `sbpy` functions and methods, we require the user to
pick among a few common field names for different properties as listed
:ref:`here <field name list>`. `~sbpy.data.DataClass` objects
are able to identify alternative field names as listed in this
document, as well as to perform transformations between a few field
names - see below for more details.

A `~sbpy.data.DataClass` object can hold as many data rows as you
want. All rows can refer to a single object, or each row can refer to
a separate object - this is usually up to the user, restrictions exist
only in a few cases as detailed in this documentation.

So what are the data containers supposed to be used for?

Ephem
^^^^^

`~sbpy.data.Ephem` has been designed to hold
**ephemerides**, i.e., properties that vary with time. 

`~sbpy.data.Ephem` currently provides convenience functions to query
ephemerides from the JPL Horizons system
(`~sbpy.data.Ephem.from_horizons`) and the Minor Planet Center
(`~sbpy.data.Ephem.from_mpc`), as well as a convenience function to
derive ephemerides from an `~sbpy.data.Orbit` object using `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_.

Obs
^^^

`~sbpy.data.Obs` is tailored to holding **observational data**, e.g.,
magnitudes as a function of time. The `~sbpy.data.Obs` class is the
only data container that is not directly derived from
`~sbpy.data.DataClass`, but from `~sbpy.data.Ephem`, providing the
same functionality as the latter.



Orbit
^^^^^

`~sbpy.data.Orbit` should be used to hold **orbital elements** of one
or several bodies. Elements can be retrieved using the convenience
function `~sbpy.data.Orbit.from_horizons`, propagated using
`~sbpy.data.Orbit.oo_propagate`, and transformed into other frames
using `~sbpy.data.Orbit.oo_transform`.

Phys
^^^^

`~sbpy.data.Phys` objects are meant to hold **physical properties**
that do not change over time. Known physical properties can currently
be queried from the JPL Small-Body Database Browser system using
`~sbpy.data.Phys.from_sbdb`.


Names
^^^^^

`~sbpy.data.Names` objects are somewhat different from the other data
containers, as they don't hold properties but only object
**names**. These names can be used to identify object nature
(`~sbpy.data.Names.asteroid_or_comet`) and they can be parsed to
extract individual identifier components
(`~sbpy.data.Names.parse_asteroid` and
`~sbpy.data.Names.parse_comet`).



How to use Ephem, Orbit, Obs, and Phys objects
----------------------------------------------

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


Building an object
^^^^^^^^^^^^^^^^^^

While `~sbpy.data.Ephem`, `~sbpy.data.Orbit`, `~sbpy.data.Obs`, and
`~sbpy.data.Phys` provide a range of convience functions to build
objects containing data, for instance from online data archives, it is
easily possible to build these objects from scratch. This can be done
for input data stored in dictionaries
(`~sbpy.data.DataClass.from_dict`), lists or arrays
(`~sbpy.data.DataClass.from_columns` and
`~sbpy.data.DataClass.from_rows`), `~astropy.table.Table` objects
(`~sbpy.data.DataClass.from_table`), or from data files
(`~sbpy.data.DataClass.from_file`).

Depending on how your input data are organized, you cean use different
options in different cases:

1. Assume that you want to build an `~sbpy.data.Orbit` object to
   propagate this orbit and obtain ephemerides. Since you are dealing
   with a single orbit, the most convenient solution might be to use a
   dictionary to build your object:

    >>> from sbpy.data import Orbit
    >>> import astropy.units as u
    >>> elements = {'a':1.234*u.au, 'e':0.1234, 'i':12.34*u.deg,
    ...             'argper': 123.4*u.deg, 'node': 45.2*u.deg,
    ...             'epoch': 2451200.5*u.d, 'true_anom':23.1*u.deg}
    >>> orb = Orbit.from_dict(elements)
    >>> orb  # doctest: +SKIP
    <QTable length=1>
       a       e       i     argper   node    epoch   true_anom
       AU             deg     deg     deg       d        deg
    float64 float64 float64 float64 float64  float64   float64
    ------- ------- ------- ------- ------- --------- ---------
      1.234  0.1234   12.34   123.4    45.2 2451200.5      23.1

   One quick note on building `~sbpy.data.DataClass` objects from
   dictionaries: dictionaries have no intrinsic order. In dictionary
   ``elements`` as defined here, there is no guarantee that ``'a'``
   will always be located before ``'e'`` when reading out the
   dictionary item by item, which happens when the data table is built
   in the background. Hence, the order of the resulting data table
   columns has to be considered random. If you want to force a
   specific order on the columns in your data table, you can use and
   `~collections.OrderedDict` instead of a simple dictionary. The
   order of elements in an `~collections.OrderedDict` will be the same
   as the order of the data table columns.

   For details on how to build objects from dictionaries, see
   `~sbpy.data.DataClass.from_dict`.

   
2. Now assume that you want to build an `~sbpy.data.Obs` object
   holding RA, Dec, and observation midtime for some target that you
   observed. In this case, you can use
   `~sbpy.data.DataClass.from_columns` as shown here:

    >>> from sbpy.data import Obs
    >>> import astropy.units as u
    >>> from numpy import array
    >>> ra = [10.223423, 10.233453, 10.243452]*u.deg
    >>> dec = [-12.42123, -12.41562, -12.40435]*u.deg
    >>> epoch = (2451523.5 + array([0.1234, 0.2345, 0.3525]))*u.d
    >>> obs = Obs.from_columns([ra, dec, epoch], names=['ra', 'dec', 't'])
    >>> obs
    <QTable length=3>
        ra       dec         t
       deg       deg         d
     float64   float64    float64
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234
    10.233453 -12.41562 2451523.7345
    10.243452 -12.40435 2451523.8525

   For details on how to build objects from lists or arrays, see
   `~sbpy.data.DataClass.from_columns` and also
   `~sbpy.data.DataClass.from_rows`, depending on whether your data is
   represented as rows or columns. Note that you could also use
   `~sbpy.data.DataClass.from_dict` by providing column data to the
   different fields.

    
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

   Please note that some formats used by
   `~sbpy.data.DataClass.from_file` are not able to identify units
   automatically (see `here <https://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers>`_ for a list of available formats). 

Accessing data
^^^^^^^^^^^^^^

In order to obtain a list of field names in a `~sbpy.data.DataClass`
object, you can use `~sbpy.data.DataClass.field_names`:

    >>> obs.field_names
    <TableColumns names=('ra','dec','t')>

Each of these columns can be accessed easily, for instance:

    >>> obs['ra']  # doctest: +SKIP
    [10.223423 10.233453 10.243452] deg

which will return an `~astropy.units.quantity.Quantity` object if that
column has a `~astropy.units.Unit` attached to it.

Similarly, if you are interested in the first set of observations in
``obs``, you can use:

    >>> obs[0]  # doctest: +SKIP
        ra       dec         t
       deg       deg         d
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234

which returns you a table with only the requested subset of the
data. In order to retrieve RA from the second observation, you can
combine both examples and do:

    >>> obs[1]['ra'] # doctest: +SKIP
    10.233453 deg


Just like in any `~astropy.table.Table` or `~astropy.table.QTable`
object, you can use slicing to obtain subset tables from your data,
for instance:

    >>> obs['ra', 'dec']  # doctest: +SKIP
    <QTable length=3>
	ra       dec
       deg       deg
    --------- ---------
    10.223423 -12.42123
    10.233453 -12.41562
    10.243452 -12.40435

    >>> obs[obs['ra'] <= 10.233453*u.deg] # doctest: +SKIP
        ra       dec         t
       deg       deg         d
    --------- --------- ------------
    10.223423 -12.42123 2451523.6234
    10.233453 -12.41562 2451523.7345

The results of these examples will be of the same data type as ``obs``
(any type derived from `~sbpy.data.DataClass`, e.g.,
`~sbpy.data.Ephem`, `~sbpy.data.Orbit`, ...)  The latter example shown
here uses a condition to filter data (only those observations with RA
less than or equal to 10.233453 degrees; note that it is necessary
here to apply ``u.deg`` to the value that all the RAs are compared
against) but selects all the columns in the original table.

If you ever need to access the actual `~astropy.table.QTable` object
that is inside each `~sbpy.data.DataClass` object, you can access it
as ``obs.table``.

Alternative field names
^^^^^^^^^^^^^^^^^^^^^^^

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

The definition of alternative field names is done in the file
``sbpy/data/__init__.py``, using the list ``fieldnames``. This list is
automatically tested for potential naming conflicts, i.e., different
properties that share the same alternative field names, and a
human-readable list is compiled upon building `sbpy`.

The full list of field names is available here:
:ref:`field name list`.

Field conversions
^^^^^^^^^^^^^^^^^

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
    <TableColumns names=('d','radius')>
    

Modifying an object
^^^^^^^^^^^^^^^^^^^

Individual elements, entire rows, and columns can be modified by
directly addressing them:

    >>> obs['ra'] # doctest: +SKIP
    [10.223423 10.233453 10.243452 10.25546  10.265425 10.25546  10.4545
     10.5656  ] deg
    >>> obs['ra'][:] = obs['ra'] + 0.1*u.deg
    >>> obs['ra'] # doctest: +SKIP
    [10.323423 10.333453 10.343452 10.35546  10.365425 10.35546  10.5545
     10.6656  ] deg

Note the specific syntax in this case (``obs['ra'][:] = ...``) that
is required by `~astropy.table.Table` if you want to replace
an entire column.

More complex data table modifications are possible by directly
accessing the underlying `~astropy.table.QTable` object as shown below.

`~sbpy.data.DataClass` provides a direct interface to the table
modification functions provided by `astropy.table.Table`:
`~astropy.table.Table.add_row`, `~astropy.table.Table.add_column`,
`~astropy.table.Table.add_columns`, etc. For instance, it is trivial to add
additional rows and columns to these objects.

Let's assume you want to add some more observations to your ``obs``
object:

    >>> obs.table.add_row([10.255460*u.deg, -12.39460*u.deg, 2451523.94653*u.d])
    >>> obs
    <QTable length=4>
	ra       dec          t      
       deg       deg          d      
     float64   float64     float64
    --------- --------- -------------
    10.323423 -12.42123  2451523.6234
    10.333453 -12.41562  2451523.7345
    10.343452 -12.40435  2451523.8525
     10.25546  -12.3946 2451523.94653
  

or if you want to add a column to your object:

    >>> from astropy.table import Column
    >>> obs.table.add_column(Column(['V', 'V', 'R', 'i'], name='filter'))
    >>> obs  # doctest: +SKIP
    <QTable length=4>
	ra       dec          t       filter
       deg       deg          d             
     float64   float64     float64     str1 
    --------- --------- ------------- ------
    10.223423 -12.42123  2451523.6234      V
    10.233453 -12.41562  2451523.7345      V
    10.243452 -12.40435  2451523.8525      R
     10.25546  -12.3946 2451523.94653      i

The same result can be achieved using the following syntax:

    >>> obs['filter2'] = ['V', 'V', 'R', 'i']  # doctest: +SKIP
    >>> obs  # doctest: +SKIP
    <QTable length=4>
	ra       dec          t       filter filter2
       deg       deg          d                     
     float64   float64     float64     str1    str1 
    --------- --------- ------------- ------ -------
    10.223423 -12.42123  2451523.6234      V       V
    10.233453 -12.41562  2451523.7345      V       V
    10.243452 -12.40435  2451523.8525      R       R
     10.25546  -12.3946 2451523.94653      i       i

Similarly, exisiting columns can be modified using:

    >>> obs['filter'] = ['g', 'i', 'R', 'V']  # doctest: +SKIP
    
Note how the `~astropy.table.Table.add_column` and
`~astropy.table.Table.add_row` functions are called from
``obs.table``. `~sbpy.data.DataClass.table` is a property that exposes
the underlying `~astropy.table.QTable` object so that the user can
directly interact with it. Please refer to the `~astropy.table.Table`
reference and
[documentation](https://docs.astropy.org/en/stable/table/index.html)
for more information on how to modify `~astropy.table.QTable` objects.


Writing object data to a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`~sbpy.data.DataClass` objects can be written to files using
`~sbpy.data.DataClass.to_file`:

    >>> obs.to_file('observations.dat')

By default, the data are written in ASCII format, but other formats
are available, too (cf. `~astropy.table.Table.write`). In order to
preserve units and meta data, we suggest to use the ``'FITS'`` format.


How to use Ephem and Obs
------------------------

As shown above (`How to use Ephem, Orbit, Obs, and Phys objects`_),
`~sbpy.data.Ephem` objects can be created on the fly. However,
`~sbpy.data.Ephem` can also be used to access ephemerides information
from remote services. For instance, the following few lines will query
ephemerides for asteroid Ceres on a given date and for the position of
Mauna Kea Observatory (IAU observatory code ``568``) from the `JPL Horizons service <https://ssd.jpl.nasa.gov/horizons.cgi>`_:

    >>> from sbpy.data import Ephem
    >>> from astropy.time import Time
    >>> epoch = Time('2018-08-03 14:20', scale='utc') # time in UT
    >>> eph = Ephem.from_horizons('Ceres',
    ...                           location='568',
    ...                           epochs=epoch)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable masked=True length=1>
    targetname       datetime_str          datetime_jd    ...  PABLat timescale
						d         ...   deg
       str7             str24                float64      ... float64    str3
    ---------- ------------------------ ----------------- ... ------- ---------
       1 Ceres 2018-Aug-03 14:20:00.000 2458334.097222222 ...  9.3473       UTC

    >>> eph.field_names  # doctest: +REMOTE_DATA
    <TableColumns names=('targetname','datetime_str','datetime_jd','H','G','solar_presence','flags','RA','DEC','RA_app','DEC_app','RA*cos(Dec)_rate','DEC_rate','AZ','EL','AZ_rate','EL_rate','sat_X','sat_Y','sat_PANG','siderealtime','airmass','magextinct','V','surfbright','illumination','illum_defect','sat_sep','sat_vis','ang_width','PDObsLon','PDObsLat','PDSunLon','PDSunLat','SubSol_ang','SubSol_dist','NPole_ang','NPole_dist','EclLon','EclLat','r','r_rate','delta','delta_rate','lighttime','vel_sun','vel_obs','elong','elongFlag','alpha','lunar_elong','lunar_illum','sat_alpha','sunTargetPA','velocityPA','OrbPlaneAng','constellation','TDB-UT','ObsEclLon','ObsEclLat','NPole_RA','NPole_DEC','GlxLon','GlxLat','solartime','earth_lighttime','RA_3sigma','DEC_3sigma','SMAA_3sigma','SMIA_3sigma','Theta_3sigma','Area_3sigma','RSS_3sigma','r_3sigma','r_rate_3sigma','SBand_3sigma','XBand_3sigma','DoppDelay_3sigma','true_anom','hour_angle','alpha_true','PABLon','PABLat','timescale')>

`~sbpy.data.Ephem.from_horizons` uses one or more target names, an
observer location in the form of an IAU observatory code, and a list
of discrete epochs or a range of epochs defined in a dictionary (see
`~sbpy.data.Ephem.from_horizons`) to query the JPL Horizons
service. Due to different requirements of the JPL Horizons service for
the epoch format, we recommend to use `~astropy.time.Time`
objects. The column names in the data table can be inquired using
`~sbpy.data.DataClass.field_names`.

`~sbpy.data.Ephem.from_horizons` is actually a wrapper around
`~astroquery.jplhorizons.HorizonsClass.ephemerides`. This function
conveniently combines the creation of a
`~astroquery.jplhorizons.HorizonsClass` query and the actual
ephemerides information retrieval into a single function. Additional
optional parameters provided to `~sbpy.data.Ephem.from_horizons` are
directly passed on to
`~astroquery.jplhorizons.HorizonsClass.ephemerides`, maintaining the
full flexibility of the latter function:

    >>> epoch1 = Time('2018-08-03 14:20', scale='utc')
    >>> epoch2 = Time('2018-08-04 07:30', scale='utc')
    >>> eph = Ephem.from_horizons('Ceres',
    ...                           location='568',
    ...                           epochs={'start': epoch1,
    ...                                   'stop': epoch2,
    ...                                   'step': '10m'},
    ...                           skip_daylight=True)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable masked=True length=26>
    targetname    datetime_str      datetime_jd    ...  PABLon   PABLat timescale
					 d         ...   deg      deg
       str7          str17            float64      ... float64  float64    str3
    ---------- ----------------- ----------------- ... -------- ------- ---------
       1 Ceres 2018-Aug-03 14:20 2458334.097222222 ...  171.275  9.3473       UTC
       1 Ceres 2018-Aug-03 14:30 2458334.104166667 ... 171.2774  9.3472       UTC
       1 Ceres 2018-Aug-03 14:40 2458334.111111111 ... 171.2798  9.3471       UTC
       1 Ceres 2018-Aug-03 14:50 2458334.118055556 ... 171.2822   9.347       UTC
       1 Ceres 2018-Aug-03 15:00       2458334.125 ... 171.2846  9.3469       UTC
       1 Ceres 2018-Aug-03 15:10 2458334.131944444 ... 171.2869  9.3468       UTC
	   ...               ...               ... ...      ...     ...       ...
       1 Ceres 2018-Aug-04 06:40 2458334.777777778 ... 171.5076  9.3369       UTC
       1 Ceres 2018-Aug-04 06:50 2458334.784722222 ... 171.5099  9.3368       UTC
       1 Ceres 2018-Aug-04 07:00 2458334.791666667 ... 171.5123  9.3367       UTC
       1 Ceres 2018-Aug-04 07:10 2458334.798611111 ... 171.5147  9.3366       UTC
       1 Ceres 2018-Aug-04 07:20 2458334.805555556 ... 171.5171  9.3365       UTC
       1 Ceres 2018-Aug-04 07:30      2458334.8125 ... 171.5195  9.3364       UTC

Note that ``skip_daylight`` is an optional parameter of
`~astroquery.jplhorizons.HorizonsClass.ephemerides` and it can be used
here as well. An additional feature of
`~sbpy.data.Ephem.from_horizons` is that you can automatically
concatenate queries for a number of objects:

    >>> eph = Ephem.from_horizons(['Ceres', 'Pallas', 12893, '1983 SA'],
    ...                           location='568',
    ...                           epochs=epoch1)  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable masked=True length=4>
	    targetname               datetime_str       ...  PABLat  timescale
							...   deg
	      str26                     str24           ... float64     str3
    -------------------------- ------------------------ ... -------- ---------
		       1 Ceres 2018-Aug-03 14:20:00.000 ...   9.3473       UTC
		      2 Pallas 2018-Aug-03 14:20:00.000 ... -20.1396       UTC
     12893 Mommert (1998 QS55) 2018-Aug-03 14:20:00.000 ...  -2.0567       UTC
    3552 Don Quixote (1983 SA) 2018-Aug-03 14:20:00.000 ...  13.3365       UTC

Please be aware that these queries are not simultaneous. The more
targets you query, the longer the query will take. Furthermore, keep
in mind that asteroids and comets have slightly different table
layouts (e.g., different magnitude systems: ``T-mag`` and ``N-mag``
instead of ``V-mag``), which will complicate the interpretation of the
data. It might be safest to query asteroids and comets separately.

Also note that the two examples shown above use different ways to
define epochs. The first example uses a dictionary that defines a
``start`` and ``stop`` epoch, as well as a ``step`` size (see
`~astroquery.jplhorizons.HorizonsClass.ephemerides` for
details). Instead of a ``step`` size, the user can also provide a
``number`` of steps as an integer; ``step`` is then automatically
derived from the interval and ``number`` in units of full minutes. The
second example uses a specific epoch as input. The
`~astropy.time.Time` object provided to ``epochs`` can also be
initialized with a list of epochs, querying multiple epochs at the
same time. Note that the total number of epochs queried using this
option should be less than a few hundred to prevent corruption of the
query (see `~astroquery.jplhorizons.HorizonsClass.ephemerides` for
details).

Observer locations can be defined as strings using offical `IAU
observatory codes
<https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html>`__
or using `~astropy.coordinates.EarthLocation` as shown in the following example:

    >>> from astropy.coordinates import EarthLocation
    >>> lowell = EarthLocation.of_site('Lowell Observatory')
    >>> eph = Ephem.from_horizons(1, epochs=Time('2018-01-01', format='iso'),
    ... 			  location=lowell) # doctest: +REMOTE_DATA
    >>> eph # doctest: +REMOTE_DATA
    <QTable masked=True length=1>
    targetname       datetime_str       datetime_jd ...  PABLon   PABLat timescale
					     d      ...   deg      deg            
       str7             str24             float64   ... float64  float64    str3  
    ---------- ------------------------ ----------- ... -------- ------- ---------
    1 Ceres 2018-Jan-01 00:00:00.000   2458119.5 ... 130.4303  9.2004       UTC

Offering almost identical functionality, the
`~sbpy.data.Ephem.from_mpc` method will retrieve ephemerides from the
Minor Planet Center:

    >>> eph = Ephem.from_mpc('2P', location='568',
    ...                      epochs={'start': '2018-10-22',
    ...                              'stop': '2018-10-26',
    ...                              'step': '1d'})  # doctest: +REMOTE_DATA
    >>> eph  # doctest: +REMOTE_DATA
    <QTable length=5>
	      Date          timescale ... Moon distance Moon altitude
				      ...      deg           deg
	     object            str3   ...    float64       float64
    ----------------------- --------- ... ------------- -------------
    2018-10-22 00:00:00.000       UTC ...          28.0         -33.0
    2018-10-23 00:00:00.000       UTC ...          41.0         -41.0
    2018-10-24 00:00:00.000       UTC ...          54.0         -48.0
    2018-10-25 00:00:00.000       UTC ...          67.0         -53.0
    2018-10-26 00:00:00.000       UTC ...          81.0         -56.0

Ephemerides can also be derived from `~Orbit` objects using `sbpy`'s
interface to `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_ with the function
`~sbpy.data.Ephem.from_oorb`. The following example computes
ephemerides for the next ten days in steps of 1 hr for Ceres as seen
from the Discovery Channel Telescope:

     >>> import numpy as np
     >>> from sbpy.data import Orbit, Ephem
     >>> from astropy.time import Time
     >>> epochs = Time.now().jd + np.arange(0, 10, 1/24)
     >>> ceres = Orbit.from_horizons('1')  # doctest: +REMOTE_DATA
     >>> eph = Ephem.from_oo(ceres, epochs, 'G37') # doctest: +SKIP 
     >>> print(eph) # doctest: +SKIP 
     <QTable length=240>
     targetname       epoch        ...           obsz               trueanom
			d          ...            AU                  deg
	str7         float64       ...         float64              float64
     ---------- ------------------ ... ----------------------- -----------------
	1 Ceres 2458519.2878717002 ...   4.886414464166933e-06 68.07980642088688
	1 Ceres 2458519.3295383668 ...  2.3814767035612583e-06  68.0893160393968
	1 Ceres 2458519.3712050337 ...  -7.136200919632962e-07 68.09882544202566
	1 Ceres 2458519.4128717002 ...   -4.18340743346679e-06 68.10833462855386
	1 Ceres 2458519.4545383668 ...  -7.786747377891423e-06 68.11784359908062
	1 Ceres 2458519.4962050337 ... -1.1273355301266719e-05 68.12735235370518
	    ...                ... ...                     ...               ...
	1 Ceres 2458529.0378717002 ...   1.093565783852335e-05 70.29915515170745
	1 Ceres 2458529.0795383668 ...  1.3089531693877277e-05  70.3086140523456
	1 Ceres 2458529.1212050337 ...  1.4402894355114437e-05 70.31807273565124
	1 Ceres 2458529.1628717002 ...  1.4786143903738891e-05 70.32753120140761
	1 Ceres 2458529.2045383668 ...  1.4213398342149963e-05 70.33698944971509
	1 Ceres 2458529.2462050337 ...  1.2724269065650384e-05 70.34644748067402

The properties computed by pyoorb and listed in the resulting table
are defined in the `pyoorb documentation
<https://github.com/oorb/oorb/tree/master/python>`_. Note that this function requires pyoorb to be installed, which is not a requirement for `sbpy`.

`~sbpy.data.Obs` works exactly like `~sbpy.data.Ephem`, but this class
will feature in the future some convenience functions to be able to
better deal with observational data.


How to use Orbit
----------------

`~sbpy.data.Orbit.from_horizons` enables the query of Solar System
body osculating elements from the `JPL Horizons service
<https://ssd.jpl.nasa.gov/horizons.cgi>`_:

    >>> from sbpy.data import Orbit
    >>> from astropy.time import Time
    >>> epoch = Time('2018-05-14', scale='utc')
    >>> elem = Orbit.from_horizons('Ceres', epochs=epoch)  # doctest: +REMOTE_DATA
    >>> elem  # doctest: +SKIP
    <QTable masked=True length=1>
    targetname datetime_jd ...         P         timescale
		    d      ...         d
       str7      float64   ...      float64         str2
    ---------- ----------- ... ----------------- ---------
       1 Ceres   2458252.5 ... 1681.218128428134        TT
    >>> elem.field_names  # doctest: +REMOTE_DATA
    <TableColumns names=('targetname','datetime_jd','datetime_str','H','G','e','q','incl','Omega','w','Tp_jd','n','M','nu','a','Q','P','timescale')>

If ``epochs`` is not set, the osculating elements for the current
epoch (current time) are queried. Similar to
`~sbpy.data.Ephem.from_horizons`, this function is a wrapper for
`~astroquery.jplhorizons.HorizonsClass.elements` and passes optional
parameter on to that function. Furthermore, it is possible to query
orbital elements for a number of targets:

    >>> epoch = Time('2018-08-03 14:20', scale='utc')
    >>> elem = Orbit.from_horizons(['3749', '2009 BR60'],
    ...                            epochs=epoch,
    ...                            refplane='earth')  # doctest: +REMOTE_DATA
    >>> elem # doctest: +SKIP
    <QTable length=2>
	  targetname         datetime_jd    ...         P         timescale
				  d         ...         d
	    str21              float64      ...      float64         str2
    --------------------- ----------------- ... ----------------- ---------
    3749 Balam (1982 BG1) 2458334.097222222 ... 1221.865723414031        TT
       312497 (2009 BR60) 2458334.097222222 ... 1221.776912893334        TT

An existing `~Orbit` instance can be transformed to a different
orbital element definition system (e.g., Keplerian, cometary,
cartesian) using `~sbpy.data.Orbit.oo_transform` or it can be
propagated into the future or past using
`~sbpy.data.Orbit.oo_propagate`. Both functions are implemented in
`sbpy` to provide an interface to `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_, a Python module
using `OpenOrb <https://github.com/oorb/oorb>`_.

In order to transform some current orbits to a state vector in
cartesian coordinates, one could use the following code:

    >>> elem = Orbit.from_horizons(['Ceres', 'Pallas', 'Vesta'])  # doctest: +REMOTE_DATA
    >>> statevec = elem.oo_transform('CART') # doctest: +SKIP 
    >>> statevec # doctest: +SKIP
    <QTable length=3>
       id             x                   y           ...    H       G    timescale
		      AU                  AU          ...   mag
      str8         float64             float64        ... float64 float64    str2
    -------- ------------------- -------------------- ... ------- ------- ---------
     1 Ceres -1.9673670927605356   -1.788869179608663 ...    3.34    0.12        TT
    2 Pallas  -2.354147777522819 -0.20413910825654025 ...    4.13    0.11        TT
     4 Vesta   2.142974769357926  -0.8590480100896669 ...     3.2    0.32        TT

Orbits can currently be transformed to the following definitions:
cartesian (``'CART'``), Keplerian (``'KEP'``), and cometary
(``'COM'``).

Orbit propagation requires the epoch to which the orbit should be
propagated to either as `~astropy.time.Time` object, or as float in
terms of Julian date. The following example propagates the current
orbit of Ceres back to year 2000:

    >>> elem = Orbit.from_horizons('Ceres')  # doctest: +REMOTE_DATA
    >>> epoch = Time('2000-01-01', format='iso')
    >>> newelem = elem.oo_propagate(epoch) # doctest: +SKIP 
    >>> newelem # doctest: +SKIP
    <QTable length=1>
       id           a                   e          ...    H       G    timescale
		    AU                             ...   mag
      str7       float64             float64       ... float64 float64    str3
    ------- ------------------ ------------------- ... ------- ------- ---------
    1 Ceres 2.7664942134894703 0.07837504303420217 ...    3.34    0.12       UTC

Note that both functions require pyoorb to be installed, which is
not a requirement for `sbpy`.

How to use Phys
---------------

`~sbpy.data.Phys` is designed to contain query physical properties for
small bodies; functions to query these properties are
available. `~sbpy.data.Phys.from_sbdb` queries the `JPL Small-body
Database Browser (SBDB) <https://ssd.jpl.nasa.gov/sbdb.cgi>`_ for physical
properties and stores the data in a `~sbpy.data.Phys` object, offering
the same functionality as all the other `~sbpy.data` functions,
including the use of `~astropy.units`.

As an example, the following code will query the properties for a
small number of asteroids:

    >>> from sbpy.data import Phys
    >>> phys = Phys.from_sbdb(['Ceres', '12893', '3552'])  # doctest: +REMOTE_DATA
    >>> phys['targetname', 'H', 'diameter'] # doctest: +SKIP
    <QTable length=3>
	    targetname            H    diameter
					  km
	      str26            float64 float64
    -------------------------- ------- --------
		       1 Ceres    3.34    939.4
     12893 Mommert (1998 QS55)    13.9    5.214
    3552 Don Quixote (1983 SA)    12.9     19.0


Please note that the SBDB database is not complete with respect to
physical properties and should be considered as a sparse dataset.

`~sbpy.data.Phys` also contains a function to query molecular data that
might be useful for various calculations such as production rate calculations.
`~sbpy.data.Phys.from_jplspec` queries the `JPL Molecular Spectroscopy Catalog
<https://spec.jpl.nasa.gov/home.html>`_ molecular properties, and stores the
data in a `~sbpy.data.Phys` object, offering the same functionality as all the
other `~sbpy.data` functions, including the use of `~astropy.units`. The results
from `~sbpy.data.Phys.from_jplspec` include the following data:

    | Transition frequency
    | Temperature
    | Integrated line intensity at 300 K
    | Partition function at 300 K
    | Partition function at designated temperature
    | Upper state degeneracy
    | Upper level energy in Joules
    | Lower level energy in Joules
    | Degrees of freedom

.. doctest-skip::

    >>> from sbpy.data.phys import Phys
    >>> import astropy.units as u
    >>> temp_estimate = 47. * u.K
    >>> transition_freq = (230.53799 * u.GHz).to('MHz')
    >>> mol_tag = '^CO$'
    >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    >>> mol_data
    <QTable length=1>
    Transition frequency Temperature ... Degrees of freedom Molecule Identifier
            MHz               K      ...
          float64          float64   ...       int64               int64
    -------------------- ----------- ... ------------------ -------------------
                230538.0        47.0 ...                  2               28001


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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Reference/API
-------------
.. automodapi:: sbpy.data
    :no-heading:
