import pytest
from ..core import required_packages, optional_packages
from ...exceptions import RequiredPackageUnavailable, OptionalPackageUnavailable


def test_requires():
    with pytest.raises(RequiredPackageUnavailable):
        required_packages("unavailable_package")


def test_requires_message():
    try:
        message = "Because these integrations are tricky."
        required_packages("unavailable_package", message=message)
    except RequiredPackageUnavailable as exc:
        assert message in str(exc)


def test_optional():
    with pytest.warns(OptionalPackageUnavailable):
        assert not optional_packages("unavailable_package")

    assert optional_packages("astropy")


def test_optional_message():
    message = "Using linear interpolation."
    with pytest.warns(OptionalPackageUnavailable, match=message) as record:
        optional_packages("unavailable_package", message=message)
