# Licensed under a 3-clause BSD style license - see LICENSE.rst

# name: expected result from parse_comet()
comets = {
    '1P/Halley': {'type': 'P', 'number': 1, 'name': 'Halley'},
    '3D/Biela': {'type': 'D', 'number': 3, 'name': 'Biela'},
    '6P/d\'Arrest': {'type': 'P', 'number': 6, 'name': 'd\'Arrest'},
    '9P/Tempel 1': {'type': 'P', 'number': 9, 'name': 'Tempel 1'},
    '73P/Schwassmann Wachmann 3 C': {'type': 'P', 'number': 73,
                                     'name': 'Schwassmann Wachmann 3 C'},
    '73P-C/Schwassmann Wachmann 3 C': {'type': 'P', 'number': 73,
                                       'fragment': 'C',
                                       'name': 'Schwassmann Wachmann 3 C'},
    '73P-BB': {'type': 'P', 'number': 73, 'fragment': 'BB'},
    '122P/de Vico': {'type': 'P', 'number': 122, 'name': 'de Vico'},
    '322P': {'type': 'P', 'number': 322},
    'X/1106 C1': {'type': 'X', 'desig': '1106 C1'},
    'P/1994 N2 (McNaught-Hartley)': {'type': 'P', 'desig': '1994 N2',
                                     'name': 'McNaught-Hartley'},
    'P/2001 YX127 (LINEAR)': {'type': 'P', 'desig': '2001 YX127',
                              'name': 'LINEAR'},
    'C/-146 P1': {'type': 'C', 'desig': '-146 P1'},
    'C/2001 A2-A (LINEAR)': {'type': 'C', 'desig': '2001 A2',
                             'fragm': 'A', 'name': 'LINEAR'},
    'C/2013 US10': {'type': 'C', 'desig': '2013 US10'},
    'C/2015 V2 (Johnson)': {'type': 'C', 'desig': '2015 V2', 'name': 'Johnson'}
}

# name: expected result from parse_asteroid()
asteroids = {
    '1': {'number': 1},
    '(2) Pallas': {'number': 2, 'name': 'Pallas'},
    '(2001) Einstein': {'number': 2001, 'name': 'Einstein'},
    '2001 AT1': {'desig': '2001 AT1'},
    '(1714) Sy': {'number': 1714, 'name': 'Sy'},
    '1714 SY': {'desig': '1714 SY'},  # not real, just for testing
    '2014 MU69': {'desig': '2014 MU69'},
    '(228195) 6675 P-L': {'number': 228195, 'desig': '6675 P-L'},
    '4101 T-3': {'desig': '4101 T-3'},
    '4015 Wilson-Harrington (1979 VA)': {'number': 4015, 'desig': '1979 VA',
                                         'name': 'Wilson-Harrington'},
    'J95X00A': {'desig': '1995 XA'},
    'K07Tf8A': {'desig': '2007 TA418'},
    'G3693': {'number': 163693}
}


def test_asteroid_or_comet():
    """Test target name identification."""
    from ..names import Names
    print(dir())
    for comet in comets:
        assert Names.asteroid_or_comet(comet) == 'comet', \
            'failed for {}'.format(comet)

    for asteroid in asteroids:
        if asteroid != '2017 U1':
            assert Names.asteroid_or_comet(asteroid) == 'asteroid', \
                'failed for {}'.format(asteroid)


def test_parse_comet():
    """Test comet name parsing."""

    from ..names import Names, TargetNameParseError
    import pytest

    for comet, result in comets.items():
        r = Names.parse_comet(comet)
        assert r == result, 'Parsed {}: {} != {}'.format(comet, r, result)

    # bad names
    with pytest.raises(TargetNameParseError):
        Names.parse_comet('73p')

    with pytest.raises(TargetNameParseError):
        Names.parse_comet('c/2001 A2')

    with pytest.raises(TargetNameParseError):
        Names.parse_comet('2001 a2')


def test_parse_asteroid():
    """Test asteroid name parsing."""

    from ..names import Names, TargetNameParseError
    import pytest

    for asteroid, result in asteroids.items():
        r = Names.parse_asteroid(asteroid)
        assert r == result, 'Parsed {}: {} != {}'.format(asteroid, r, result)

    # bad names
    with pytest.raises(TargetNameParseError):
        Names.parse_asteroid('C/2015 V2 (Johnson)')

    with pytest.raises(TargetNameParseError):
        Names.parse_asteroid('1P/Halley')

    with pytest.raises(TargetNameParseError):
        Names.parse_asteroid('1P')

    with pytest.raises(TargetNameParseError):
        Names.parse_asteroid('P/1994 N2')

    with pytest.raises(TargetNameParseError):
        Names.parse_asteroid('2001 at1')
