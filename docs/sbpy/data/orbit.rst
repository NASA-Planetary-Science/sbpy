=============
 Using Orbit
=============

Orbit Queries
=============

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

       Alternatively, orbital elements can also be queried from the `Minor Planet Center <https://minorplanetcenter.net/iau/MPEph/MPEph.html>`_, although in this case only the most recent elements are accessible:

    >>> elem = Orbit.from_mpc(['3552', '12893']) # doctest: +SKIP
    >>> elem # doctest: +SKIP
    <QTable length=2>
     absmag    Q      arc       w     ...     a        Tj   moid_uranus moid_venus
      mag      AU      d       deg    ...     AU                 AU         AU
    float64 float64 float64  float64  ...  float64  float64   float64    float64
    ------- ------- ------- --------- ... --------- ------- ----------- ----------
       12.9   7.278 12955.0 316.44802 ... 4.2589272     2.3     11.7518    0.56105
       13.9   3.028 12990.0 184.31369 ... 2.8281991     3.3     15.1738    1.90535


Orbit Transformations
=====================
       
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

Orbit Propagations
==================

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
