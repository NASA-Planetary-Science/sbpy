"""
This is the toy version of the toy model
MM
"""

import numpy as np
import astropy.units as u


class Photometry():

    def apparent_magnitudes(self):
        if self.eph == None:
            raise ValueError('No ephemerides available')

        mag = []
        # use simple H-G model from Bowell, Ast2, for now
        for epoch in self.eph:
            alpha = np.deg2rad(epoch['alpha'])
            phi_1 = np.exp(-3.33*pow(np.tan(0.5*alpha), 0.63))
            phi_2 = np.exp(-1.87*pow(np.tan(0.5*alpha), 1.22))

            mag.append(self.absmag-2.5*np.log10((1-self.slopepar)*phi_1
                                                +self.slopepar*phi_2)
                        +5.0*np.log10(epoch['r']*epoch['delta']))
        return np.array(mag)


class Data():
    
    # use a decorator to extract phys properties and populate Body attributes
    def get_ephemerides(self, obscode, date):
        
        # use callhorizons for now, will use new method
        import callhorizons

        eph = callhorizons.query(self.ident)
        eph.set_discreteepochs(date)
        eph.get_ephemerides(obscode)
        self.eph = eph

        # workaround for missing decorator
        self.absmag = eph[0]['H'] 
        self.slopepar = eph[0]['G']

        return len(eph)
        

class Activity():

    def __init__(self):
        self.afrho = None
        self.dustrate = None
        self.dustmodel = None
    
    def model_dust_coma(self, dust_rate, dust_size):
        self.dustmodel = dust_rate*dust_size # bogus

    def surface_brightness(self, aprad):
        if self.dustmodel is None:
            raise ValueError('No dust coma model available')
        
        return self.dustmodel*np.pi*aprad**2
    
class Body(Photometry, Data, Activity):

    def __init__(self, ident):
        # initialize basic physical and orbital properties
        self.ident = ident
        self.diam = None
        self.pv = None
        self.eph = None
        self.absmag = None
        self.slopepar = None
        
    def __repr__(self):
        return(self.ident)



    
ceres = Body('ceres')
ceres.get_ephemerides('568', ['2451234.123456', '2451254.123456', '2451274.123456'])
print('derived magnitudes and Horizons magnitudes:')
print(ceres.apparent_magnitudes())
print(ceres.eph['V'])

# add dust coma
ceres.model_dust_coma(1, 2)
print('coma surface brightness (bogus):', ceres.surface_brightness(2))
totalmag = (ceres.apparent_magnitudes()
            +ceres.surface_brightness(2)
            -2.5*np.log10(np.pi*2**2))
print('total magnitude in 2 arcsec aperture (more bogus):', totalmag)


