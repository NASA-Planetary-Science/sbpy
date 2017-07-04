# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
SBPy data.names Module
======================

Functions concerning small body names and designations.

created on June 04, 2017
"""

__all__ = ['altident', 'asteroid_or_comet', 'from_packed', 'to_packed', 'parse_comet', 'parse_asteroid', 'TargetNameParseError']

class TargetNameParseError(Exception):
    pass

def altident(targetid, bib=None):
    """Query Lowell database to obtain alternative target names for `targetitd`.

    Examples
    --------
    >>> from sbpy.data import names
    >>> names.altident('3552')
    ['3552', 'Don Quixote', '1983 SA']

    not yet implemented

    """

def asteroid_or_comet(targetid):
    """Checks if an object is an asteroid, or a comet, based on its targetid.

    Parameters
    ----------
    targetid : string
      The name of the target.
    
    Returns
    -------
    target_type : string
      The target identification: 'comet', 'asteroid', or `None`.

    Examples
    --------
    >>> from sbpy.data import names
    >>> names.asteroid_or_comet('2P')
    'comet'
    >>> names.asteroid_or_comet('(1) Ceres')
    'asteroid'
    >>> names.asteroid_or_comet('Fred')
    None

    """

    try:
        parse_comet(targetid)
        return 'comet'
    except TargetNameParseError:
        pass

    try:
        parse_asteroid(targetid)
        return 'asteroid'
    except TargetNameParseError:
        pass

    return None
        
def from_packed(pkd):
    """Convert MPC packed designation/number to unpacked identifier.

    Parameters
    ----------
    pkd : string
      The packed designation.
    
    Returns
    -------
    desig : string
      Unpacked designation.

    Examples
    --------
    >>> from sbpy.data import from_packed
    >>> names.from_packed('J95A010')
    '1995 A1'

    not yet implemented

    """

def to_packed(targetid):
    """Convert target designation/number to packed identifier.

    Parameters
    ----------
    pkd : string
      The packed designation.
    
    Returns
    -------
    desig : string
      Unpacked designation.

    Examples
    --------
    >>> from sbpy.data import from_packed
    >>> names.from_packed('J95A010')
    '1995 A1'

    not yet implemented

    """

def parse_comet(s):
    """Parse a string as if it were a comet name.

    Only considers IAU-formatted permanent and new-style designations.
    Note that letter case is important.

    Parameters
    ----------
    s : string
      The string to parse.

    Returns
    -------
    des : string
      The designation of the comet, e.g., '1P', or 'C/1995 O1'.
    name : string
      The name of the comet, if provided in `s`, e.g., 'Halley', or
      'Hale-Bopp'.

    Raises
    ------
    TargetNameParseError : Exception
      If the string does not appear to be a comet name.

    Examples
    --------
    >>> from sbpy.data import names
    >>> names.parse_comet('9P/Tempel 1')
    ('9P', 'Tempel 1')
    >>> names.parse_comet('C/2001 A2-A (LINEAR)')
    ('C/2001 A2-A', 'LINEAR')

    The following table shows results of the designation parsing:

      +-------------------------------+-------------+-------------------------+
      |targetname                     |des          |name                     |
      +===============================+=============+=========================+
      |1P/Halley                      |1P           |Halley                   |
      +-------------------------------+-------------+-------------------------+
      |3D/Biela                       |3D           |Biela                    |
      +-------------------------------+-------------+-------------------------+
      |9P/Tempel 1                    |9P           |Tempel 1                 |
      +-------------------------------+-------------+-------------------------+
      |73P/Schwassmann-Wachmann 3 C   |73P          |Schwassmann-Wachmann 3 C |
      +-------------------------------+-------------+-------------------------+
      |73P-C/Schwassmann-Wachmann 3 C |73P-C        |Schwassmann-Wachmann 3 C |
      +-------------------------------+-------------+-------------------------+
      |73P-BB                         |73P-BB       |                         |
      +-------------------------------+-------------+-------------------------+
      |322P                           |322P         |                         |
      +-------------------------------+-------------+-------------------------+
      |X/1106 C1                      |X/1106 C1    |                         |
      +-------------------------------+-------------+-------------------------+
      |P/1994 N2 (McNaught-Hartley)   |P/1994 N2    |McNaught-Hartley         |
      +-------------------------------+-------------+-------------------------+
      |P/2001 YX127 (LINEAR)          |P/2001 YX127 |LINEAR                   |
      +-------------------------------+-------------+-------------------------+
      |C/-146 P1                      |C/-146 P1    |                         |
      +-------------------------------+-------------+-------------------------+
      |C/2001 A2-A (LINEAR)           |C/2001 A2-A  |LINEAR                   |
      +-------------------------------+-------------+-------------------------+
      |C/2013 US10                    |C/2013 US10  |                         |
      +-------------------------------+-------------+-------------------------+
      |C/2015 V2 (Johnson)            |C/2015 V2    |Johnson                  |
      +-------------------------------+-------------+-------------------------+

    """

    import re

    pat = ('((^([1-9][0-9]*[PD](-[A-Z]{1,2})?)(/(.+))?$)'
           '|(^([CPX]/-?[0-9]{1,4} [A-Z]{1,2}[1-9][0-9]{0,2}(-[A-Z]{1,2})?)\s*(\((.+)\))?$))')

    m = re.findall(pat, s.strip())
    if len(m) > 0:
        m = m[0]
        if len(m[2]) > 0:    # Permament
            return m[2], m[5]
        elif len(m[7]) > 0:  # Provisional
            return m[7], m[10]

    raise TargetNameParseError('{} does not appear to be a comet name'.format(s))

def parse_asteroid(s):
    """Parse a string as if it were an asteroid name.

    Only considers IAU-formatted permanent and new-style designations.
    Note that letter case is important.

    Parameters
    ----------
    s : string
      The string to parse.

    Returns
    -------
    des : string
      The designation of the asteroid, e.g., '1', or '2014 MU69'.
    name : string
      The name of the asteroid, if provided in `s`, e.g., 'Ceres'.

    >>> from sbpy.data import names
    >>> names.parse_asteroid('(1) Ceres')
    ('1', 'Ceres')
    >>> names.parse_asteroid('2014 MU69')
    ('2014 MU69', '')

    Examples
    --------
    The following table shows the result of the parsing:

      +-----------------+-----------+----------+
      |targetname       |des        |name      |
      +=================+===========+==========+
      |1                |1          |          |
      +-----------------+-----------+----------+
      |(2001) Einstein  |2001       |Einstein  |
      +-----------------+-----------+----------+
      |2001 AT1         |2001 AT1   |          |
      +-----------------+-----------+----------+
      |(1714) Sy        |1714       |Sy        |
      +-----------------+-----------+----------+
      |1714 SY          |1714 SY    |          |  # compare with previous
      +-----------------+-----------+----------+
      |2014 MU69        |2014 MU69  |          |
      +-----------------+-----------+----------+
      |2017 AA          |2017 AA    |          |
      +-----------------+-----------+----------+
      |(20231) 1997 YK  |20231      |1997 YK   |
      +-----------------+-----------+----------+
      |2040 P-L         |2040 P-L   |          |
      +-----------------+-----------+----------+
      |3138 T-1         |3138 T-1   |          |
      +-----------------+-----------+----------+
      |1010 T-2         |1010 T-2   |          |
      +-----------------+-----------+----------+
      |4101 T-3         |4101 T-3   |          |
      +-----------------+-----------+----------+

    """

    import re

    # Provisional, permanent, then survey-specific
    pat = ('((^([1-9]|A)[0-9]{3,3}( [A-Z]{2,2}([1-9][0-9]{0,2})?)$)'
           '|(^\(([1-9][0-9]*)\) (.+)?$)|(^[1-9][0-9]*$)'
           '|(^[1-9][0-9]* ((P-L)|(T-[123]))$))')

    m = re.findall(pat, s.strip())
    if len(m) > 0:
        m = m[0]
        if len(m[1]) > 0:    # Provisional
            return m[1], ''
        elif len(m[6]) > 0:  # Permanent with name
            return m[6], m[7]
        elif len(m[8]) > 0:  # Permenent without name
            return m[8], ''
        elif len(m[9]) > 0:  # Survey-specific
            return m[9], ''

    raise TargetNameParseError('{} does not appear to be an asteroid name'.format(s))
