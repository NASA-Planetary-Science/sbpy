from warnings import catch_warnings
import pytest
from ..decorators import requires, optional
from ...exceptions import RequiredPackageUnavailable, OptionalPackageUnavailable

def test_requires():
    @requires("unavailable_package")
    def f():
        pass

    with pytest.raises(RequiredPackageUnavailable):
        f()

def test_optional():
    @optional("unavailable_package")
    def f():
        pass

    with pytest.warns(OptionalPackageUnavailable):
        f()
