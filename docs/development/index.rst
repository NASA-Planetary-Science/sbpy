Contributing to `sbpy`
======================

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
<https://guides.github.com/features/issues/>`_. This way, all
communications are centralized in one place and accessible to
everyone. `sbpy` issues can be posted `here
<https://github.com/NASA-Planetary-Science/sbpy/issues>`_.

Another way to communicate in a less official way (if you don't think
your question or suggestion is worth being remembered in the far
future) is the `sbpy slack channel <http://sbpy.slack.com>`_ (requires
sign-in). In order to request access to the slack channel, please
email sbpy.dev(at)gmail.com.

Testing
-------

We are looking for volunteers to test `sbpy` functionality. In order
to make this task worthwhile for you, we only ask you to test routines
that are stable and which you can already use for your research. We
will indicate which functions are ready to be tested on this dedicated
`status page <status.rst>`_. Keep in mind that there might be issues,
so feel free to compare the function results with your own results. If
you find a problem, please create an `issue report
<https://github.com/NASA-Planetary-Science/sbpy/issues>`__.

Code Contributions
------------------

If you would like to implement (or modify) some functionality
yourself, you are welcome to do so, but please adhere to `astropy's
contributing guildelines <http://www.astropy.org/contribute.html>`__;
this link also provides a general guide to the workflow involving git
and other hints.

This is the proposed workflow for code contributions:

* Before you write any code, please inform the developer via an `issue
  report <https://github.com/NASA-Planetary-Science/sbpy/issues>`__ of
  your plans, so that duplication of work can be prevented.
* Given that `sbpy` is supported by a NASA grant, the developers have
  to maintain the package's stability and integrity; hence, the
  developers' decisions have to be followed at any time.
* If your plans have been approved by the developers, they will point
  you to a place in the code where your contribution can be
  implemented.
* Fork the latest version of the `sbpy` repository from
  `https://github.com/NASA-Planetary-Science/sbpy
  <https://github.com/NASA-Planetary-Science/sbpy>`__, clone it to
  your machine, and create a new branch with a descriptive name
  relevant to your project. Implement your code into this new branch.
* To make your code compatible with the `sbpy` API, please follow the
  `function/method design guidelines
  <https://github.com/NASA-Planetary-Science/sbpy/wiki/function-method-design>`_.
* All `sbpy` code contributions require docstrings in agreement with
  the `sbpy` `docstring guidelines
  <https://github.com/NASA-Planetary-Science/sbpy/wiki/docstring-guidelines>`__,
  unit tests, and useful `documentation
  <http://docs.astropy.org/en/latest/development/docguide.html>`__.
  Code should be written following the `astropy coding guidelines
  <http://docs.astropy.org/en/latest/development/codeguide.html>`__.
* `sbpy` follows the `astropy testing guidelines
  <http://docs.astropy.org/en/latest/development/testguide.html>`__.
* Once you are done working on your code, push your branch to your
  forked repository on github and create a pull request describing
  your code and linking to your original issue report.
  
Please follow the `astropy code of conduct`_ at any time.

.. _astropy code of conduct: http://docs.astropy.org/en/latest/development/codeguide.html
