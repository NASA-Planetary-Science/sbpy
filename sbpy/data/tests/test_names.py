# Licensed under a 3-clause BSD style license - see LICENSE.rst

# name: expected result from parse_comet()
comets = {
    '1P/Halley': ('1P', 'Halley'),
    '3D/Biela': ('3D', 'Biela'),
    '9P/Tempel 1': ('9P', 'Tempel 1'),
    '73P/Schwassmann-Wachmann 3 C': ('73P', 'Schwassmann-Wachmann 3 C'),
    '73P-C/Schwassmann-Wachmann 3 C': ('73P-C', 'Schwassmann-Wachmann 3 C'),
    '73P-BB': ('73P-BB', ''),
    '322P': ('322P', ''),
    'X/1106 C1': ('X/1106 C1', ''),
    'P/1994 N2 (McNaught-Hartley)': ('P/1994 N2', 'McNaught-Hartley'),
    'P/2001 YX127 (LINEAR)': ('P/2001 YX127', 'LINEAR'),
    'C/-146 P1': ('C/-146 P1', ''),
    'C/2001 A2-A (LINEAR)': ('C/2001 A2-A', 'LINEAR'),
    'C/2013 US10': ('C/2013 US10', ''),
    'C/2015 V2 (Johnson)': ('C/2015 V2', 'Johnson'),
}

# name: expected result from parse_asteroid()
asteroids = {
    '1': ('1', ''),
    '(2) Pallas': ('2', 'Pallas'),
    '(2001) Einstein': ('2001', 'Einstein'),
    '2001 AT1': ('2001 AT1', ''),
    '(1714) Sy': ('1714', 'Sy'),
    '1714 SY': ('1714 SY', ''), # not real; note confusion with prev. item
    '2014 MU69': ('2014 MU69', ''),
    '2017 AA': ('2017 AA', ''),
    '(20231) 1997 YK': ('20231', '1997 YK'),
    '2040 P-L': ('2040 P-L', ''),
    '3138 T-1': ('3138 T-1', ''),
    '1010 T-2': ('1010 T-2', ''),
    '4101 T-3': ('4101 T-3', ''),
}

def test_asteroid_or_comet():
    """Test target name idenfication."""

    from ..names import asteroid_or_comet
    
    for comet in comets:
        assert asteroid_or_comet(comet) == 'comet'

    for asteroid in asteroids:
        assert asteroid_or_comet(asteroid) == 'asteroid'

    for name in ['Fred', 'S/2005 P 1']:
        assert asteroid_or_comet(name) is None

def test_parse_comet():
    """Test comet name parsing."""

    from ..names import parse_comet, TargetNameParseError
    import pytest

    for comet, result in comets.items():
        r = parse_comet(comet)
        assert r == result, 'Parsed {}: {} != {}'.format(comet, r, result)

    # bad names
    with pytest.raises(TargetNameParseError):
        parse_comet('(1) Ceres')
        
    with pytest.raises(TargetNameParseError):
        parse_comet('73p')
        
    with pytest.raises(TargetNameParseError):
        parse_comet('c/2001 A2')

    with pytest.raises(TargetNameParseError):
        parse_comet('C/2001 a2')

    with pytest.raises(TargetNameParseError):
        parse_comet('C/2015 V2 Johnson')

def test_parse_asteroid():
    """Test asteroid name parsing."""

    from ..names import parse_asteroid, TargetNameParseError
    import pytest

    for asteroid, result in asteroids.items():
        r = parse_asteroid(asteroid) 
        assert r == result, 'Parsed {}: {} != {}'.format(asteroid, r, result)

    with pytest.raises(TargetNameParseError):
        parse_asteroid('C/2015 V2 (Johnson)')

    with pytest.raises(TargetNameParseError):
        parse_asteroid('1P/Halley')

    with pytest.raises(TargetNameParseError):
        parse_asteroid('1P')

    with pytest.raises(TargetNameParseError):
        parse_asteroid('P/1994 N2')

    with pytest.raises(TargetNameParseError):
        parse_asteroid('2001 at1')
