======================================
Orbit Data Objects (`sbpy.data.Orbit`)
======================================

Orbit Queries
=============

`~sbpy.data.Orbit.from_horizons` enables the query of Solar System
body osculating elements from the `JPL Horizons service
<https://ssd.jpl.nasa.gov/horizons/>`_:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> from sbpy.data import Orbit
    >>> from astropy.time import Time
    >>> epoch = Time('2018-05-14', scale='utc')
    >>> elem = Orbit.from_horizons('Ceres', epochs=epoch)
    >>> elem  # doctest: +SKIP
    <...>/sbpy/data/orbit.py:113: TimeScaleWarning: converting utc epochs to tdb for use in astroquery.jplhorizons
      TimeScaleWarning)
    <QTable masked=True length=1>
    targetname    H       G    ...         P               epoch      
		 mag           ...         d                          
       str7    float64 float64 ...      float64            object     
    ---------- ------- ------- ... ----------------- -----------------
       1 Ceres    3.34    0.12 ... 1681.218136185659 2458252.500800755

Please note the ``TimeScaleWarning``, which is raised when the time
scale of the desired epoch is not supported by the query function and
hence converted to a scale that is supported (``tdb`` in this case).
The following fields are available in this object:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> elem.field_names
    ['targetname', 'H', 'G', 'e', 'q', 'incl', 'Omega', 'w', 'n', 'M', 'nu', 'a', 'Q', 'P', 'epoch', 'Tp']

If ``epochs`` is not set, the osculating elements for the current
epoch (current time) are queried. Similar to
`~sbpy.data.Ephem.from_horizons`, this function is a wrapper for
`~astroquery.jplhorizons.HorizonsClass.elements` and passes optional
parameter on to that function. Furthermore, it is possible to query
orbital elements for a number of targets:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> epoch = Time('2018-08-03 14:20', scale='tdb')
    >>> elem = Orbit.from_horizons(['3749', '2009 BR60'],
    ...                            epochs=epoch,
    ...                            refplane='earth')
    >>> elem # doctest: +SKIP
    <QTable masked=True length=2>
	  targetname         H       G    ...         P               epoch      
			    mag           ...         d                          
	    str21         float64 float64 ...      float64            object     
    --------------------- ------- ------- ... ----------------- -----------------
    3749 Balam (1982 BG1)    13.3    0.15 ... 1221.865723743203 2458334.097222222
       312497 (2009 BR60)    17.7    0.15 ... 1221.776912893334 2458334.097222222

Alternatively, orbital elements can also be queried from the `Minor
Planet Center <https://minorplanetcenter.net/iau/MPEph/MPEph.html>`_,
although in this case only the most recent elements are accessible:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> elem = Orbit.from_mpc(['3552', '12893'])
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

.. .. doctest-requires:: astroquery, oorb
.. doctest-remote-data::

    >>> elem = Orbit.from_horizons(['Ceres', 'Pallas', 'Vesta'])
    >>> statevec = elem.oo_transform('CART') # doctest: +SKIP 
    >>> statevec # doctest: +SKIP
    <QTable length=3>
       id             x                   y          ...    H       G   
		      AU                  AU         ...   mag          
      str8         float64             float64       ... float64 float64
    -------- ------------------- ------------------- ... ------- -------
     1 Ceres -0.4867631007775121 -2.7702346649193696 ...    3.34    0.12
    2 Pallas -1.7745931352186222 -1.7169356664520194 ...    4.13    0.11
     4 Vesta    2.24552918427612  1.0169886872736296 ...     3.2    0.32

Orbits can currently be transformed to the following definitions:
cartesian (``'CART'``), Keplerian (``'KEP'``), and cometary
(``'COM'``).

Orbit Propagations
==================

Orbit propagation requires the epoch to which the orbit should be
propagated to either as `~astropy.time.Time` object, or as float in
terms of Julian date. The following example propagates the current
orbit of Ceres back to year 2000:

.. .. doctest-requires:: astroquery, oorb
.. doctest-remote-data::

    >>> elem = Orbit.from_horizons('Ceres')
    >>> epoch = Time('2000-01-01', scale='tdb')
    >>> newelem = elem.oo_propagate(epoch) 
    >>> newelem # doctest: +SKIP
    <QTable length=1>
       id           a                  e          ...   epoch      H       G   
		    AU                            ...             mag          
      str7       float64            float64       ...   object  float64 float64
    ------- ----------------- ------------------- ... --------- ------- -------
    1 Ceres 2.766494220549446 0.07837504411299284 ... 2451544.5    3.34    0.12

Note that both functions require `pyoorb
<https://github.com/oorb/oorb/tree/master/python>`_ to be installed.


Calculate dynamical parameters
==============================

The Tisserand parameter is a commonly used dynamic parameter to characterize
the orbit of a small body, especially a comet, when its orbital evolution is
dominated by the gravitational effect of a particular planet.  The Tisserand
parameter with respect to Jupiter is used in the dynamical classification of
comets.  The Tisserand parameter can be calculated by `~sbpy.Orbit.tisserand`
as follows:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> epoch = Time(2449400.5, format='jd', scale='tdb')
    >>> halley = Orbit.from_horizons('1P', id_type='designation',
    ...     closest_apparition=True, epochs=epoch)
    >>> T = halley.tisserand()
    >>> print('{:.4f}'.format(T)) # doctest: +SKIP
    -0.6050

One can also specify the planet with respect to which the Tisserand parameter
is calculated with optional parameter `planet`.  It also allows multiple
planet to be specified simultaneously:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> chariklo = Orbit.from_dict({
    ...     "e": 0.16778,
    ...     "a": 15.78694 * u.au,
    ...     "incl": 23.39153 * u.deg,
    ...     "Omega": 300.44770 * u.deg,
    ...     "w": 242.01787 * u.deg,
    ...     "n": 0.015713 * u.deg / u.day,
    ...     "M":  113.36375 * u.deg,
    ... })
    >>> T = chariklo.tisserand(planet=['599', '699', '799', '899'], epoch=Time("2023-08-25"))
    >>> with np.printoptions(precision=3):
    ...     print(T)  # doctest: +FLOAT_CMP
    [3.482 2.93  2.859 3.225]

`~sbpy.Orbit` also provides a method to compare the orbits of two objects
in terms of the "D-criterion" (`Jopek 1993 <https://ui.adsabs.harvard.edu/abs/1993Icar..106..603J/abstract>`_).  The `~sbpy.Orbit.D_criterion` method
implements all three versions of the D-criterion, including
Southworth & Hawkins function (`Southworth and Hawkins 1963 <https://ui.adsabs.harvard.edu/abs/1963SCoA....7..261S/abstract>`_),
Drummond function (`Drummond 1991 <https://ui.adsabs.harvard.edu/abs/1981Icar...45..545D/abstract>`_), and the hybrid function (`Jopek 1993 <https://ui.adsabs.harvard.edu/abs/1993Icar..106..603J/abstract>`_).
The code example below demonstrates the calculation of three versions of
D_criterion:

.. .. doctest-requires:: astroquery
.. doctest-remote-data::

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>>
    >>> comet_252P = Orbit.from_dict({  # P/2016 BA14
    ...     "e": 0.67309,
    ...     "q": 0.99605 * u.au,
    ...     "incl": 10.42220 * u.deg,
    ...     "Omega": 190.94850 * u.deg,
    ...     "w": 343.31047 * u.deg,
    ...     "n": 0.18533 * u.deg / u.day,
    ...     "Tp": Time("2016-03-15 06:19:30", scale="tdb")
    ... })
    >>>
    >>> comet_460P = Orbit.from_dict({  # P/2016 BA14
    ...     "e": 0.66625,
    ...     "q": 1.00858 * u.au,
    ...     "incl": 18.91867 * u.deg,
    ...     "Omega": 180.53368 * u.deg,
    ...     "w": 351.89672 * u.deg,
    ...     "n": 0.18761 * u.deg / u.day,
    ...     "Tp": Time("2016-03-15 12:24:19", scale="tdb")
    ... })
    >>>
    >>> # Southworth & Hawkins function
    >>> D_SH = comet_252P.D_criterion(comet_460P)
    >>> # Drummond function
    >>> D_D = comet_252P.D_criterion(comet_460P, version='d')
    >>> # hybrid function
    >>> D_H = comet_252P.D_criterion(comet_460P, version='h')
    >>> print('D_SH = {:.4f}, D_D = {:.4f}, D_H = {:.4f}'.
    ...    format(D_SH, D_D, D_H))
    D_SH = 0.1562, D_D = 0.0503, D_H = 0.1558

