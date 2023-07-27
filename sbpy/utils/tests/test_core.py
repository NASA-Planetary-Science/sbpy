import pytest
from ..core import requires, optional
from ...exceptions import RequiredPackageUnavailable, OptionalPackageUnavailable


def test_requires():
    with pytest.raises(RequiredPackageUnavailable):
        requires("unavailable_package")


def test_requires_message():
    try:
        message = "Because these integrations are tricky."
        requires("unavailable_package", message=message)
    except RequiredPackageUnavailable as exc:
        assert message in str(exc)


def test_optional():
    with pytest.warns(OptionalPackageUnavailable):
        assert not optional("unavailable_package")

    assert optional("astropy")


def test_optional_message():
    message = "Using linear interpolation."
    with pytest.warns(OptionalPackageUnavailable, match=message) as record:
        optional("unavailable_package", message=message)
