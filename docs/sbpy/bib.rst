Bibliography Tracking Module (`sbpy.bib`)
=========================================

Introduction
------------

`sbpy` classes and functions can automatically report citations to the
user through a bibliography registry. The idea behind this service is
to make it easier to properly acknowledge and reference those who
designed methods and tools used.

How to use `~sbpy.bib`
----------------------

Use `~sbpy.bib.track` to enable citation tracking. Every method called
after activating the tracking will register citations relevant to
it. Each citation is associated with a tag that enables the user to
identify which aspect of the method the citation is relevant to:

    >>> from sbpy import bib, data
    >>> bib.track()
    >>> eph = data.Ephem.from_horizons('433', epochs=None, location='500')
    >>> print(bib.to_text())  # doctest: +REMOTE_DATA
    sbpy.data.Ephem:
      data service: Giorgini et al. 1996, 1996DPS....28.2504G

In this case, ``Giorgini et al. 1996, 1996DPS....28.2504G`` is
relevant to the implementation of the JPL Horizons system that is
queried by the `~sbpy.data.Ephem.from_horizons`
function. `~sbpy.bib.to_text` outputs the current citation registry in
simple text form.


Bibliography tracking can also be used in a context manager:

    >>> from sbpy import bib
    >>> from sbpy.data import Ephem
    >>> with bib.Tracking():
    ...     eph = Ephem.from_horizons('Ceres', epochs=None, location='500')
    >>> print(bib.to_text())  # doctest: +REMOTE_DATA
    sbpy.data.Ephem:
      data service:
        Giorgini et al. 1996, 1996DPS....28.2504G


Output formats
--------------

Bibliographies can be generated in different output formats:

Simple text (`~sbpy.bib.to_text`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    >>> bib.to_text()  # doctest: +REMOTE_DATA
    sbpy.data.Ephem:
      data service:
        Giorgini et al. 1996, 1996DPS....28.2504G


BibTeX (`~sbpy.bib.to_bibtex`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**not yet implemented**


Other bibliography styles, for instance in accordance to specific
journal rules, can be readily implemented


Reference/API
-------------
.. automodapi:: sbpy.bib
    :no-heading:
