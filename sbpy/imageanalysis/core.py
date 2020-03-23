# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy Imageanalysis Module

created on June 23, 2017
"""

__all__ = ['CometaryEnhancement', 'PSFSubtraction', 'centroid']


def centroid(im, gyx, method='weighted'):
    """Method to centroid on a target used different methods

    Parameters
    ----------
    im : array, mandatory
        image
    gyx : tuple, mandatory
        initial guess coordinates (y, x)
    method : str, optional, default: `weighted`, choices: [`weighted`,
        `peak`, `comet`]
        method to use for finding centroid

    Returns
    -------
    tuple, (y,x) image coordinates

    Examples
    --------
    >>> from astropy.io import fits
    >>> from sbpy.imageanalysis import centroid # doctest: +SKIP
    >>> hdu = fits.open('test.fits') # doctest: +SKIP
    >>> cyx = centroid(hdu[0].data, (21, 45)) # doctest: +SKIP

    not yet implemented

    """


class CometaryEnhancement():
    def __init__(self, im, cyx):
        self.im = im
        self.cyx = cyx

    def azavg_norm(self):
        """Normalize image using azimuthal average

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometaryEnhancement # doctest: +SKIP
        >>> import matplotlib.pyplot as plt # doctest: +SKIP
        >>> hdu = fits.open('test.fits') # doctest: +SKIP
        >>> cyx = centroid(hdu[0].data, (21, 45)) # doctest: +SKIP
        >>> enhance = CometaryEnhancement(hdu[0].data, cyx) # doctest: +SKIP
        >>> plt.imshow(enhance.azavg_norm()) # doctest: +SKIP

        not yet implemented

        """

    def rho_norm(self):
        """Normalize image by a 1/rho profile.

        Returns
        -------
        enhanced : ndarray
          The enhanced image.

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometaryEnhancement # doctest: +SKIP
        >>> import matplotlib.pyplot as plt # doctest: +SKIP
        >>> hdu = fits.open('test.fits') # doctest: +SKIP
        >>> cyx = centroid(hdu[0].data, (21, 45)) # doctest: +SKIP
        >>> enhance = CometaryEnhancement(hdu[0].data, cyx) # doctest: +SKIP
        >>> plt.imshow(enhance.rho_norm()) # doctest: +SKIP

        """

        from .utils import rarray

        r = rarray(self.im.shape, yx=self.cyx, subsample=10)

        return self.im * r

    def rvsf(self, **parameters):
        """Normalize image using a radially variable spatial filter.

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometaryEnhancement # doctest: +SKIP
        >>> import matplotlib.pyplot as plt # doctest: +SKIP
        >>> hdu = fits.open('test.fits') # doctest: +SKIP
        >>> cyx = centroid(hdu[0].data, (21, 45)) # doctest: +SKIP
        >>> enhance = CometaryEnhancement(hdu[0].data, cyx) # doctest: +SKIP
        >>> plt.imshow(enhance.rvsf_norm(a=1, b=1, n=0.1)) # doctest: +SKIP

        not yet implemented

        """

    # how to implement Ginga interface?


class PSFSubtraction():

    def create_psfmodel(im, cyx, **parameters):
        """Create a PSF model, ignoring the actual target

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, PSFSubtraction # doctest: +SKIP
        >>> import matplotlib.pyplot as plt # doctest: +SKIP
        >>> hdu = fits.open('test.fits') # doctest: +SKIP
        >>> cyx = centroid(hdu[0].data, (21, 45)) # doctest: +SKIP
        >>> psfmodel = PSFSubtraction(hdu[0].data, cyx) # doctest: +SKIP
        >>> plt.imshow(psfmodel) # doctest: +SKIP

        not yet implemented

        """

    def subtract_psfmodel(im, psfmodel, cyx, scaling='peak'):
        """subtract a PSF model from the image

        Parameters
        ----------
        im : array, mandatory
            image
        psfmodel : array, mandatory
            PSF model image
        cyx : tuple, mandatory
            target centroid coordinates `(y, x)`
        scaling : str, optional, default: `peak`
            select scaling method (`peak` fits the height of PSF model to
            the target's image)


        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, PSFSubtraction # doctest: +SKIP
        >>> import matplotlib.pyplot as plt # doctest: +SKIP
        >>> hdu = fits.open('test.fits') # doctest: +SKIP
        >>> cyx = centroid(hdu[0].data, (21, 45)) # doctest: +SKIP
        >>> psfmodel = PSFSubtraction.create_psfmodel(hdu[0].data, cyx) # doctest: +SKIP
        >>> diff = PSFSubtraction.subtract_psfmodel(hdu[0].data, psfmodel, cyx) # doctest: +SKIP
        >>> plt.imshow(diff) # doctest: +SKIP

        not yet implemented

        """
