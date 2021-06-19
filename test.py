import astropy.units as u
from sbpy.activity import gas

Q = 1e25 / u.s
v = 1 * u.km / u.s
parent = gas.photo_lengthscale('H2O')
daughter = gas.photo_lengthscale('OH')
coma = gas.Haser(Q, v, parent, daughter)

from sbpy.activity import CircularAperture
ap = CircularAperture(1 * u.arcsec).as_length(1 * u.au)
N = coma.total_number(ap)