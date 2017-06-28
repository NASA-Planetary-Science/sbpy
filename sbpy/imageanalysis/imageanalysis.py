# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
=========================
SBPy Imageanalysis Module
=========================

created on June 23, 2017
"""

__all__ = ['CometEnhancement', 'PSFSubtraction', 'centroid']



def centroid(im, gyx, method='weighted'):
    """Method to centroid on a target used different methods
   
    Parameters
    ----------
    im : array, mandatory
        image 
    gyx : tuple, mandatory
        initual guess coordinates (y, x)
    method : str, optional, default: `weighted`, choices: [`weighted`, `peak`, `comet`]
        method to use for finding centroid
    
    Returns
    -------
    tuple, (y,x) image coordinates

    Examples
    --------
    >>> from astropy.io import fits
    >>> from sbpy.imageanalysis import centroid
    >>> hdu = fits.open('test.fits')
    >>> cyx = centroid(hdu[0].data, (21, 45))
    
    not yet implemented

    """

    
class CometEnhancement():
    def __init__(self, im, cyx):
        self.im = im
        self.cyx = cyx

    def azavg_norm(self):
        """Normalize image using azimuthal average 

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometEnhancement
        >>> import matplotlib.pyplot as plt
        >>> hdu = fits.open('test.fits')
        >>> cyx = centroid(hdu[0].data, (21, 45))
        >>> enhance = CometEnhancement(hdu[0].data, cyx)
        >>> plt.imshow(enhance.azavg_norm())

        not yet implemented

        """

    def rho_norm(self):
        """Normalize image using 1/rho fit

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometEnhancement
        >>> import matplotlib.pyplot as plt
        >>> hdu = fits.open('test.fits')
        >>> cyx = centroid(hdu[0].data, (21, 45))
        >>> enhance = CometEnhancement(hdu[0].data, cyx)
        >>> plt.imshow(enhance.rho_norm())

        not yet implemented

        """

    def rho_norm(self, **parameters):
        """Normalize image using a radially variable spatial filter

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometEnhancement
        >>> import matplotlib.pyplot as plt
        >>> hdu = fits.open('test.fits')
        >>> cyx = centroid(hdu[0].data, (21, 45))
        >>> enhance = CometEnhancement(hdu[0].data, cyx)
        >>> plt.imshow(enhance.rvsf_norm(a=1, b=1, n=0.1))

        not yet implemented

        """
        
    # how to implement Ginga interface?

    
class PSFSubtraction():

    def create_psfmodel(im, cyx, **parameters):
        """Create a PSF model, ignoring the actual target

        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometEnhancement
        >>> import matplotlib.pyplot as plt
        >>> hdu = fits.open('test.fits')
        >>> cyx = centroid(hdu[0].data, (21, 45))
        >>> psfmodel = PSFSubtraction(hdu[0].data, cyx)
        >>> plt.imshow(psfmodel)

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
            select scaling method (`peak` fits the height of PSF model to the target's image)


        Examples
        --------
        >>> from astropy.io import fits
        >>> from sbpy.imageanalysis import centroid, CometEnhancement
        >>> import matplotlib.pyplot as plt
        >>> hdu = fits.open('test.fits')
        >>> cyx = centroid(hdu[0].data, (21, 45))
        >>> psfmodel = PSFSubtraction.create_psfmodel(hdu[0].data, cyx)
        >>> diff = PSFSubtraction.subtract_psfmodel(hdu[0].data, psfmodel, cyx) 
        >>> plt.imshow(diff)

        not yet implemented 

        """
