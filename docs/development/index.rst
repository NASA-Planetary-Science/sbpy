.. _contributing:

Contributing to sbpy
====================

If you are interested in contributing to `sbpy`, you can do that
through testing existing code (somewhat easy), contributing Python
code (somewhat hard), or simply by letting us know what you think. We
are also interested in code donations: if you have code written in
Python or a different programming language and you think this code
might be useful to others, please let us know.

Communicating with the developers
---------------------------------

We would prefer any communication with the development team to go
through the `github issue system
<https://docs.github.com/issues>`_. This way, all
communications are centralized in one place and accessible to
everyone. `sbpy` issues can be posted `here
<https://github.com/NASA-Planetary-Science/sbpy/issues>`__.

Another way to communicate in a less official way (if you don't think
your question or suggestion is worth being remembered in the far
future) is the sbpy slack channel, which is part of the Astropy slack
(see `here <https://www.astropy.org/help.html>`__ for information on
how to join).

Testing
-------

We are looking for volunteers to test `sbpy` functionality. In order
to make this task worthwhile for you, we only ask you to test routines
that are stable and which you can already use for your research. We
will indicate which functions are ready to be tested on this dedicated
:doc:`/status`. Keep in mind that there might be issues, so feel free to
compare the function results with your own results.

`sbpy` tests are run with `pytest`.  To install `pytest` and all
requirements for testing:

.. code-block:: bash

  pip install sbpy[test]

or install `sbpy` from the source tree in editable mode:

.. code-block:: bash

  pip install -e .[test]

Then the tests may be run

.. code-block:: bash

  pytest sbpy
  # or pytest sbpy --remote-data to enable tests requiring an internet connection

For more testing options, including testing multiple dependency versions,
see `astropy`'s `testing guidelines <https://docs.astropy.org/en/latest/development/testguide.html>`__.

Reporting Problems
------------------

If you find a problem, please create an `issue report
<https://github.com/NASA-Planetary-Science/sbpy/issues>`__. An issue
report template is available there.
 

Code Contributions
------------------

If you would like to implement (or modify) some functionality
yourself, you are welcome to do so, but please make sure that your
contribution meets the following rules and requirements.

Contribution Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~

Topical requirements
^^^^^^^^^^^^^^^^^^^^

* Contributions must have general relevance to astronomers studying
  small solar System bodies, especially comets and asteroids.
* Small contributions are welcome, especially if they can be extended
  in the future.
* Highly specialized code with limited use cases is not within the scope of `sbpy`.
* Functionality that has a wider interest in astronomy should instead be
  considered by other Astropy affiliated packages or Astropy itself.

.. _technical requirements:
  
Technical requirements
^^^^^^^^^^^^^^^^^^^^^^

* code must adhere to `astropy's contributing guidelines
  <https://www.astropy.org/contribute.html>`__, the guidelines
  described in this document and `PEP8
  <https://peps.python.org/pep-0008/>`_
* code must be accompanied by corresponding tests; 100% of the
  implemented tests must pass, a test coverage >= 90% is required; if
  possible, results should be checked against results from the
  literature
* code must be accompanied by docstrings that describe the input and
  output parameters and includes example code, documentation, and at
  least one science task-oriented notebook that goes into the `sbpy
  tutorial repository
  <https://github.com/NASA-Planetary-Science/sbpy-tutorial>`_
* the API used in the code must follow sbpy's :doc:`design-principles`:

    1. physical parameters must be `~astropy.units.Quantity`;
    2. epochs are `~astropy.time.Time` objects;
    3. use sbpy `~sbpy.data.DataClass` objects: `~sbpy.data.Orbit`, `~sbpy.data.Phys`, `~sbpy.data.Ephem`, and `~sbpy.data.Obs`;
    4. append fields to ``DataClass`` at the user's request via ``append_results``;
    5. cite relevant works with :func:`sbpy.bib.register` or `~sbpy.bib.cite`;
    6. except for private functions or speed.

* consider class method names following the pattern ``.to_XXX`` and ``.from_XXX``
* add references to docstrings, documentation, and tests where
  applicable; use `~sbpy.bib`!
* customized exceptions and warnings are encouraged, and should
  ultimately be derived from the base classes in `sbpy.exceptions`
* if you use `~sbpy.data.DataClass` objects, extend the :ref:`field
  name list` list where it makes sense
* consider using the sbpy function and decorator helpers to test for the
  presence of optional dependencies:
  * `~sbpy.utils.required_packages` and `~@sbpy.utils.decorators.requires` raise an exception if a package cannot be imported.
  * `~sbpy.utils.optional_packages` and `~@sbpy.utils.decorators.optionally_uses` warn the user if a package cannot be imported.
* a CHANGELOG entry is required; also update the :doc:`/status` where applicable



Contribution Workflow
~~~~~~~~~~~~~~~~~~~~~

This is the proposed workflow for code contributions:

* Before you write any code, please issue a `Feature Request
  <https://github.com/NASA-Planetary-Science/sbpy/issues/new?assignees=&labels=feature+request&template=feature_request.md&title=feature+request>`_,
  fill out the form and submit it.
* The `sbpy` team and any interested parties will discuss the proposal
  in the issue comments to (1) determine if it is within the scope of
  `sbpy`, (2) to avoid duplication of effort and functionality,
  and (3) to determine what the best location within `sbpy` is.
* A final decision on the proposal will be made by the `sbpy` core
  team and, if accepted, coding may begin, and a pull request made.
* To make your code compatible with the `sbpy` API, please follow the
  :ref:`technical requirements` for new code.
* The pull request will be merged after it successfully passed a
  review process conducted by at least one `sbpy` core developer team
  member.

Please also check out `astropy's contributing guidelines
<https://www.astropy.org/contribute.html>`__ for a general introduction
on coding techniques and additional hints.

Please follow the `astropy code of conduct`_ at any time.

.. _astropy code of conduct: https://docs.astropy.org/en/latest/development/codeguide.html
