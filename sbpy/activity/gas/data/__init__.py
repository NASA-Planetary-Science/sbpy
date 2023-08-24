# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Gas Data"""

__all__ = [
    'photo_lengthscale',
    'photo_timescale',
    'fluorescence_band_strength',
    'OHFluorescenceSA88'
]

import numpy as np

import astropy.units as u
from astropy.io import ascii
from astropy.utils.data import get_pkg_data_path

from .... import data as sbd
from .... import bib
from ....utils import optional_packages

photo_lengthscale = {   # (value, {key feature: ADS bibcode})
    'H2O': {
        'CS93': (2.4e4 * u.km,
                 {'H2O photodissociation lengthscale':
                  '1993Icar..105..235C'})
    },
    'OH': {
        'CS93': (1.6e5 * u.km,
                 {'OH photodissociation lengthscale':
                  '1993Icar..105..235C'})
    },
}

photo_timescale = {   # (value, {key feature: ADS bibcode})
    'H2O': {
        'CS93': (5.2e4 * u.s,
                 {'H2O photodissociation timescale':
                  '1993Icar..105..235C'})
    },
    'OH': {
        'CS93': (1.6e5 * u.s,
                 {'OH photodissociation timescale':
                  '1993Icar..105..235C'})
    },
    'HCN': {
        'C94': (6.7e4 * u.s,
                {'HCN photodissociation timescale':
                 '1994JGR....99.3777C'})
    },
    'CH3OH': {
        'C94': (7.7e4 * u.s,
                {'CH3OH photodissociation timescale':
                 '1994JGR....99.3777C'})
    },
    'H2CO': {
        'C94': (5.0e3 * u.s,
                {'H2CO photodissociation timescale':
                 '1994JGR....99.3777C'})
    },
    'CO': {
        'CE83': (1.5e6 * u.s,
                 {'CO photodissociation timescale':
                  '1983A&A...126..170C'})
    },
    'CO2': {
        'CE83': (5.0e5 * u.s,
                 {'CO2 photodissociation timescale':
                  '1983A&A...126..170C'})
    },
    'CN': {
        'H92': ([3.15e5, 1.35e5] * u.s,
                {'CN photodissociation timescale':
                 '1992Ap&SS.195....1H'})
    },
}


class OHFluorescenceSA88:
    """OH fluorescence.

    Valid for heliocentric radial velocities between -60 and 60 km/s,
    and heliocentric distances greater than 0.5 au.

    Based on Table 5 of Schleicher and A'Hearn 1988, The Fluorescence
    of Cometary OH, ApJ 331, 1058-1077.


    Parameters
    ----------
    band : string
        Initialize for this OH band.  Valid bands may be found via
        `OHFluorescenceSchleicher88.BANDS`.


    Examples
    --------

    >>> from sbpy.activity.gas.data import OHFluorescenceSA88
    >>> LN = OHFluorescenceSA88('0-0')
    >>> print(LN(-1 * u.km / u.s))    # doctest: +FLOAT_CMP
    [1.54e-15] erg / s

    """

    BANDS = ['0-0', '1-0', '1-1', '2-2',
             '0-1', '0-2', '1-2', '2-0', '2-1']

    def __init__(self, band):
        fn = get_pkg_data_path('schleicher88.txt')
        self.table5 = ascii.read(fn)
        self._rdot = self.table5['rdot'].data * u.km / u.s
        self._inversion = self.table5['I']
        self._tau = self.table5['tau'] * u.s

        self.basis = ['0-0', '1-0', '1-1', '2-2',
                      '0-0', '0-0', '1-1', '2-2', '2-2']
        self.scales = [1.0, 1.0, 1.0, 1.0,
                       0.00356, 0.00021, 0.00610, 0.274, 1.921]

        i = self.BANDS.index(band)
        k = self.basis[i]
        self._LN = u.Quantity(self.table5[k].data * self.scales[i],
                              'erg / s')

        if optional_packages("scipy"):
            from scipy.interpolate import splrep
            self._tck = splrep(self.rdot.value, self.LN.value)
            self._interp = self._spline
        else:
            self._interp = self._linear

    def _spline(self, rdot, rh):
        from scipy.interpolate import splev
        return splev(rdot, self._tck, ext=2) / rh**2

    def _linear(self, rdot, rh):
        return np.interp(rdot, self.rdot.value, self.LN.value) / rh**2

    @bib.cite({'OH fluorescence band efficiency': '1988ApJ...331.1058S'})
    @sbd.dataclass_input(eph=sbd.Ephem)
    @sbd.quantity_to_dataclass(eph=(sbd.Ephem, 'rdot'))
    def __call__(self, eph):
        """Fluorescence band efficiency.

        Evaluated at the given heliocentric radial velocity, and,
        optionally, heliocentric distance.


        Parameters
        ----------
        eph : `~astropy.units.Quantity`, `~sbpy.data.Ephem`
            Heliocentric radial velocity as a
            `~astropy.units.Quantity` or a column in an
            `sbpy.data.Ephem` object.  If heliocentric distance is
            present, the value will be scaled by ``rh**-2``,
            otherwise, the value at 1 au is returned.


        Returns
        -------
        LN : `~astropy.units.Quantity`
            Fluorescence band efficiency or luminosity per molecule.


        Raises
        ------
        `ValueError` when the heliocentric distance is < 0.5 au or
        r-dot is outside the tabulated range.

        """

        rdot = eph['rdot'].to(self.rdot.unit).value
        if np.any(np.abs(rdot) > 60):
            raise ValueError('r-dot must be between -60 and 60 km/s')

        try:
            rh = eph['rh'].to(u.au).value
        except KeyError:
            rh = 1

        if np.any(rh < 0.5):
            raise ValueError(
                'At rh < 0.5 au the pumping rate is not small compared '
                'to the rotational decay rate.  See Schleicher & A\'Hearn '
                '1988 for details.')

        return self._interp(rdot, rh) * self.LN.unit

    @property
    def rdot(self):
        """Heliocentric radial velocity for tabulated data."""
        return self._rdot

    @property
    def inversion(self):
        """Inversion (n_u - n_l) / (n_u + n_l)."""
        return self._inversion

    @property
    def tau(self):
        """Lifetime via A^2 Σ^+(ν=2,3)."""
        return self._tau

    @property
    def LN(self):
        """Tabulated fluorescence band efficiency (L/N)."""
        return self._LN


fluorescence_band_strength = {
    # (function, note, citation or None)
    # for OHFluorescenceSA88, the citation is handled by the class
    'OH 0-0': {
        'SA88': (OHFluorescenceSA88('0-0'), 'Requires r-dot', None)
    },
    'OH 1-0': {
        'SA88': (OHFluorescenceSA88('1-0'), 'Requires r-dot', None)
    },
    'OH 1-1': {
        'SA88': (OHFluorescenceSA88('1-1'), 'Requires r-dot', None)
    },
    'OH 2-2': {
        'SA88': (OHFluorescenceSA88('2-2'), 'Requires r-dot', None)
    },
    'OH 0-1': {
        'SA88': (OHFluorescenceSA88('0-1'), 'Requires r-dot', None)
    },
    'OH 0-2': {
        'SA88': (OHFluorescenceSA88('0-2'), 'Requires r-dot', None)
    },
    'OH 1-2': {
        'SA88': (OHFluorescenceSA88('1-2'), 'Requires r-dot', None)
    },
    'OH 2-0': {
        'SA88': (OHFluorescenceSA88('2-0'), 'Requires r-dot', None)
    },
    'OH 2-1': {
        'SA88': (OHFluorescenceSA88('2-1'), 'Requires r-dot', None)
    },
}
