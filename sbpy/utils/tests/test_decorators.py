from warnings import catch_warnings
import pytest
from ..decorators import requires, optionally_uses
from ...exceptions import RequiredPackageUnavailable, OptionalPackageUnavailable


@requires("unavailable_package")
def f_required():
    pass


@optionally_uses("unavailable_package")
def f_optional():
    pass


def test_requires():
    with pytest.raises(RequiredPackageUnavailable):
        f_required()


def test_optional():
    with pytest.warns(OptionalPackageUnavailable):
        f_optional()
