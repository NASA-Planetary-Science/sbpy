========================================
Physical Data Objects (`sbpy.data.Phys`)
========================================

`~sbpy.data.Phys` is designed to contain and query physical properties for
small bodies; functions to query these properties are
available. `~sbpy.data.Phys.from_sbdb` queries the `JPL Small-body
Database Browser (SBDB) <https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html>`_ for physical
properties and stores the data in a `~sbpy.data.Phys` object, offering
the same functionality as all the other `~sbpy.data` functions,
including the use of `~astropy.units`.

As an example, the following code will query the properties for a
small number of asteroids:

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> from sbpy.data import Phys
    >>> phys = Phys.from_sbdb(['Ceres', '12893', '3552'])
    >>> phys['targetname', 'H', 'diameter'] # doctest: +SKIP
    <QTable length=3>
            targetname            H    diameter
                                 mag      km   
              str26            float64 float64 
    -------------------------- ------- --------
             1 Ceres (A801 AA)    3.56    939.4
     12893 Mommert (1998 QS55)   13.98    5.214
    3552 Don Quixote (1983 SA)   12.96     19.0


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

.. .. doctest-requires:: astroquery
.. doctest-remote-data:: 

    >>> from sbpy.data.phys import Phys
    >>> import astropy.units as u
    >>> temp_estimate = 47. * u.K
    >>> transition_freq = (230.53799 * u.GHz).to('MHz')
    >>> mol_tag = '^CO$'
    >>> mol_data = Phys.from_jplspec(temp_estimate, transition_freq, mol_tag)
    >>> mol_data  # doctest: +SKIP
    <QTable length=1>
    Transition frequency Temperature ... Degrees of freedom Molecule Identifier
            MHz               K      ...
          float64          float64   ...       int64               int64
    -------------------- ----------- ... ------------------ -------------------
                230538.0        47.0 ...                  2               28001

