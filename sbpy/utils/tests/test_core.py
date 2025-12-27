import pytest
import numpy as np
from astropy.utils.masked import Masked
from ..core import required_packages, optional_packages, unmasked
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


def test_unmasked():
    array = [1, 2, 3]
    mask = [False, False, True]
    ma1 = Masked([1.0, 2.0, 3.0], mask=[False, False, True])
    ma2 = np.ma.MaskedArray(array, mask=mask)

    assert unmasked(array) is array
    assert (unmasked(ma1) == array).all()
    assert (unmasked(ma2) == array).all()
