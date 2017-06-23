Bibliography Module (`sbpy.bib`)
================================

Introduction
------------

`sbpy.bib` is the reference tracking system of `SBPy`. The idea
is to provide the user with a list of bibliographic references that
should be cited if `SBPy` functionality has been used in their
research. These references can be compiled as simple text, or as bibtex
elements for use in LATeX.

Every top-level function of `SBPy` has an optional `bib` argument to
which a bibliography instance can be provided. By passing a
bibliography instance to each function used, a complete list of
references can be generated.

Reference/API
-------------
.. automodapi:: sbpy.bib
    :no-heading:
