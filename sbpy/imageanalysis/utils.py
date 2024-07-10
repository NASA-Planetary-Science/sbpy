# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
===================================
SBPy Imageanalysis Module Utilities
===================================

created on July 02, 2017
"""


def rarray(shape, yx=None, subsample=0):
    """2D array of distances from a point.

    Parameters
    ----------
    shape : array
      The shape of the resulting array, (y, x).
    yx : array, optional
      The center of the array, (y, x).  If set to `None`, then the
      center is `(shape - 1.0) / 2.0`.  Integer values refer to the
      center of the pixel.
    subsample : int, optional
      Set to `>1` to sub-pixel sample the 11x11 pixel core of the
      array using `refine_pixel`.

    Returns
    -------
    r : ndarray
      The array of radial values.

    Example
    -------
    >>> from sbpy.imageanalysis.utils import rarray
    >>> r = rarray((5, 5))
    >>> r[2, 2]  # docest: +FLOAT_CMP
    0.0
    >>> r = rarray((5, 5), yx=(0, 0))
    >>> r[2, 2]  # docest: +FLOAT_CMP
    2.8284271247461903

    """

    import numpy as np

    if yx is None:
        yx = (np.array(shape) - 1) / 2
    else:
        yx = np.array(yx)

    y = yarray(shape, yx=yx)
    x = xarray(shape, yx=yx)
    r = np.sqrt(x**2 + y**2)

    if subsample > 1:
        Ny = min(shape[0], 11)
        Nx = min(shape[1], 11)
        y0, x0 = np.floor((yx + 0.5)).astype(int)
        for y in range(-Ny // 2 + y0, Ny // 2 + y0 + 1):
            for x in range(-Nx // 2 + x0, Nx // 2 + x0 + 1):
                f = refine_pixel(rarray, subsample, (y, x), yx) * subsample
                r[y, x] = f
    return r


def rebin(a, factor, flux=False, trim=False):
    """Rebin a 1, 2, or 3 dimensional array by integer amounts.

    Parameters
    ----------
    a : ndarray
      Image to rebin.
    factor : int
      Rebin factor.  For `factor < 0`, shrink the image.  When
      shrinking, all axes must be an integer multiple of factor or use
      `trim`.
    flux : bool
      If `True`, preserve flux, otherwise, surface brightness is
      preserved.
    trim : bool
      If `True`, trim the shape of the input array to an integer
      multiple of factor.

    Returns
    -------
    b : ndarray
      The rebinned array.

    Example
    -------
    >>> from sbpy.imageanalysis.utils import xarray, rebin
    >>> x = xarray((10, 10))
    >>> x[0, :2].mean()  # docest: +FLOAT_CMP
    0.5
    >>> x2 = rebin(x, -2)
    >>> x2[0, 0]  # docest: +FLOAT_CMP
    0.5

    """

    import numpy as np

    def mini(a, factor):
        b = a[::-factor]
        for i in range(-factor - 1):
            b += a[(i + 1)::-factor]
        if not flux:
            b = b / -factor
        return b

    def magni(a, factor):
        s = np.array(a.shape)
        s[0] *= factor
        b = np.zeros(s)
        for i in range(factor):
            b[i::factor] = a
        if flux:
            b = b / factor
        return b

    if factor == 1:
        # done!
        return a

    _a = a.copy()
    if factor < 0:
        for i in range(_a.ndim):
            if trim:
                r = _a.shape[i] % abs(factor)
                if r != 0:
                    _a = np.rollaxis(np.rollaxis(_a, i)[:-r], 0, i + 1)

            if (_a.shape[i] % factor) != 0:
                raise ValueError(
                    "Axis {0} must be an integer multiple of "
                    "the minification factor.".format(i)
                )
        f = mini
    else:
        f = magni

    b = f(_a, factor)
    for i in range(len(_a.shape) - 1):
        c = f(np.rollaxis(b, i + 1), factor)
        b = np.rollaxis(c, 0, i + 2)

    return b


def refine_pixel(func, subsample, yx_pixel, yx, **kwargs):
    """Subsample an array-generating function over a pixel.

    The function is numerically averaged over the area.

    Parameters
    ----------
    func : function
      The array-generating function.  The first argument must be
      `shape`.  See `xarray` for an example.
    subsample : int
      The subsample factor.
    yx_pixel : array of int
      The coordinates of the pixel to consider.
    yx : array
      The function's origin.
    **kwargs
      Keyword arguments to pass on to `func`.

    Returns
    -------
    f : float
      The average pixel value over the subsampled area.

    Example
    -------
    >>> import numpy as np
    >>> from sbpy.imageanalysis.utils import rarray, refine_pixel
    >>> yx = (2, 2)  # the center of the radial array
    >>> r = rarray((5, 5), yx=yx)
    >>> r[2, 2]  # docest: +FLOAT_CMP
    0.0
    >>> f = refine_pixel(rarray, 10, (2, 2), yx)
    >>> np.isclose(f, 0.03826, rtol=0.01, atol=0.01)
    True
    """

    import numpy as np

    yx_s = ((np.array(yx) - np.array(yx_pixel)) + 0.5) * subsample - 0.5
    f = func((subsample, subsample), yx=yx_s, **kwargs) / subsample**2
    return f.mean()


def xarray(shape, yx=[0, 0], th=None):
    """2D array of x values.

    Parameters
    ----------
    shape : tuple of int
      The shape of the resulting array, (y, x).
    yx : tuple of float
      Offset the array to align with this y, x center.
    th : Quantity, optional
      Place the x-axis along this position angle, measured
      counterclockwise from the original x-axis.

    Returns
    -------
    x : ndarray
      An array of x values.

    Examples
    --------
    >>> from sbpy.imageanalysis.utils import xarray
    >>> x = xarray((10, 10))
    >>> x[0, 3]
    3

    """

    import numpy as np
    import astropy.units as u

    y, x = np.indices(shape)[-2:]
    y = y - yx[0]
    x = x - yx[1]

    if th is not None:
        x = x * np.cos(th.to(u.rad).value) + y * np.sin(th.to(u.rad).value)

    return x


def yarray(shape, yx=[0, 0], th=None):
    """2D array of y values.

    Parameters
    ----------
    shape : tuple of int
      The shape of the resulting array, (y, x).
    yx : tuple of float
      Offset the array to align with this y, x center.
    th : Quantity, optional
      Place the y-axis along this position angle, measured
      counterclockwise from the original y-axis.

    Returns
    -------
    y : ndarray
      An array of y values.

    >>> from sbpy.imageanalysis.utils import yarray
    >>> y = yarray((10, 10))
    >>> y[3, 0]
    3

    """

    import numpy as np
    import astropy.units as u

    y, x = np.indices(shape)[-2:]
    y = y - yx[0]
    x = x - yx[1]

    if th is not None:
        y = -x * np.sin(th.to(u.rad).value) + y * np.cos(th.to(u.rad).value)

    return y
