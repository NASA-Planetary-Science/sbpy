def test_rarray():
    import numpy as np
    from ..utils import rarray

    r = rarray((5, 5))
    assert r[0, 0] == 2 * np.sqrt(2)
    assert r[2, 2] == 0.0

    r = rarray((5, 5), yx=(0, 0))
    assert r[2, 2] == 2 * np.sqrt(2)
    assert r[0, 0] == 0.0

    r = rarray((5, 5), subsample=10)
    assert np.isclose(r[2, 2], 0.3826, rtol=0.01, atol=0.01)


def test_rebin():
    import numpy as np
    from ..utils import rebin

    a = np.ones((10, 10))
    assert rebin(a, -10) == 1.0
    assert rebin(a, -10, flux=True) == 100

    b = rebin(a, 2)
    assert b.shape == (20, 20)
    assert b[0, 0] == 1.0

    b = rebin(a, 2, flux=True)
    assert b[0, 0] == 0.25

    b = rebin(a, -3, trim=True)
    assert b.shape == (3, 3)


def test_refine_pixel():
    import numpy as np
    from ..utils import rarray, refine_pixel

    yx = (2, 2)
    f = refine_pixel(rarray, 10, (2, 2), yx)
    assert np.isclose(f, 0.03826, rtol=0.01, atol=0.01)

    f = refine_pixel(rarray, 10, (2, 1), yx)
    assert np.isclose(f, 0.10434, rtol=0.01, atol=0.01)

    f = refine_pixel(rarray, 10, (1, 1), yx)
    assert np.isclose(f, 0.14435, rtol=0.01, atol=0.01)


def test_xarray():
    import numpy as np
    import astropy.units as u
    from ..utils import xarray, yarray

    x = xarray((3, 3))
    xx = np.array([[0, 1, 2], [0, 1, 2], [0, 1, 2]])
    assert np.all(x == xx)

    x = xarray((3, 3), yx=(1, 0))
    assert np.all(x == xx)

    x = xarray((3, 3), yx=(0, 0.5))
    xx = xx - 0.5
    assert np.all(x == xx)

    x = xarray((3, 3), th=90 * u.deg)
    y = yarray((3, 3))
    assert np.allclose(x, y)


def test_yarray():
    import numpy as np
    import astropy.units as u
    from ..utils import yarray, xarray

    y = yarray((3, 3))
    yy = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    assert np.all(y == yy)

    y = yarray((3, 3), yx=(0, 1))
    assert np.all(y == yy)

    y = yarray((3, 3), yx=(0.5, 0))
    yy = yy - 0.5
    assert np.all(y == yy)

    y = yarray((3, 3), th=-90 * u.deg)
    x = xarray((3, 3))
    assert np.allclose(y, x)
