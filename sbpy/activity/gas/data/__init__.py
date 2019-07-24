# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Gas Data"""

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
                  '1983A%26A...126..170C'})
    },
    'CO2': {
        'CE83': (5.0e5 * u.s,
                 {'CO2 photodissociation timescale':
                  '1983A%26A...126..170C'})
    },
    'CN': {
        'H92': ([3.15e5, 1.35e5] * u.s,
                {'CN photodissociation timescale':
                 '1992Ap%26SS.195....1H'})
    },
}

data = {   # (value, {key feature: bibcode})
    'OH 0-0': {
        'SA88': (func0_0,
                 {'OH 0-0 fluorescence band efficiency':
                  '1988ApJ...331.1058S'})
    },
    'OH 1-0': {
        'SA88': (func1_0,
                 {'OH 1-0 fluorescence band efficiency':
                  '1988ApJ...331.1058S'})
    },
    'OH 1-1': {
        'SA88': (func1_1,
                 {'OH 1-1 fluorescence band efficiency':
                  '1988ApJ...331.1058S'})
    },
    'OH 2-2': {
        'SA88': (func2_2,
                 {'OH 2-2 fluorescence band efficiency':
                  '1988ApJ...331.1058S'})
    },
}
