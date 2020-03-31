==============
Using DASTCOM5
==============


For using the DASTCOM5 Module, you have to first download the databse locally.
That can be done by:

    >>> from sbpy.data.utils import dastcom5
    >>> dastcom5.download_dastcom5()  # doctest: +SKIP

After the database is downloaded, all the queries can be done easily.

DASTCOM5 is a subset of Small Body Database provided by JPL, NASA.
For querying the database, either name or record number for the object
can be used.

    >>> dastcom5.orbit_from_name('atira')  # doctest: +SKIP
    >>> dastcom5.orbit_from_record(900001)  # doctest: +SKIP

More information about the DASTCOM5 Database can be taken from it's README file.
