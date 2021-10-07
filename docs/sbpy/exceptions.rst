Exceptions Module (`sbpy.exceptions`)
=====================================

Introduction
------------

`~sbpy.exceptions` provides the ``SbpyExecption`` and ``SbpyWarning`` base classes, from which all sbpy-specific exceptions and warnings should be derived.

To suppress all sbpy warnings:

```
import warnings
from sbpy.exceptions import SbpyWarning
warnings.simplefilter('ignore', SbpyWarning)
```

See the `Astropy documentation <https://docs.astropy.org/en/stable/warnings.html>`_ for more ways to suppress warnings.


Reference/API
-------------
.. automodapi:: sbpy.exceptions
    :no-heading:
