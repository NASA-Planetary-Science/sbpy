# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Names Module
======================

Class for dealing with object naming conventions.

created on August 28, 2017

"""

from ..exceptions import SbpyException
from .core import DataClass
from numpy import ndarray

__all__ = ['Names', 'TargetNameParseError', 'natural_sort_key']


def natural_sort_key(s):
    """List sort keys considering strings of numbers as integers.

    Intended to be used with `list.sort` or the `sorted` built-in
    function.


    Parameters
    ----------
    s : string
        String to parse into keys.


    Returns
    -------
    keys : tuple
        Keys for sorting.


    Examples
    --------
    >>> from sbpy.data.names import natural_sort_key
    >>> comets = ['9P/Tempel 1',
    ...           '101P/Chernykh',
    ...           '10P/Tempel 2',
    ...           '2P/Encke']
    >>> sorted(comets)
    ['101P/Chernykh', '10P/Tempel 2', '2P/Encke', '9P/Tempel 1']
    >>> sorted(comets, key=natural_sort_key)
    ['2P/Encke', '9P/Tempel 1', '10P/Tempel 2', '101P/Chernykh']

    """
    import re
    keys = tuple()
    for k in re.split('([0-9]+)', str(s)):
        keys += (int(k) if k.isdigit() else k,)
    return keys


class TargetNameParseError(SbpyException):
    pass


class Names():
    """Class for parsing target identifiers. The functions in this class will
    identify designation, name strings, and number for both comets and
    asteroids. It also includes functionality to distinguish between comet and
    asteroid identifiers."""

    # packed numbers translation string
    pkd = ('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
           'abcdefghifklmnopqrstuvwxyz')

    @staticmethod
    def to_packed(s):
        """Convert asteroid designation/number to packed identifier.

        Parameters
        ----------
        s : str
           Target identifier.

        Returns
        -------
        p : str
           Packed designation/number.

        Examples
        --------
        >>> from sbpy.data import Names
        >>> Names.to_packed('1995 AA1')
        'J95A01A'
        """

        if s.isdigit() and not s.isalpha():
            # number
            s = int(s)
            if s < 100000:
                return ('{:05d}'.format(s))
            elif s > 619999:
                raise TargetNameParseError(
                    ('{} cannot be turned into a '
                     'packed number').format(s))
            else:
                mod = (s % 10000)
                return ('{}{:04d}'.format(
                    Names.pkd[int((s-mod)/10000)],
                    mod))
        elif s.isalnum and not s.isdigit() and not s.isalpha():
            # designation
            yr = s.strip()[:4]
            yr = Names.pkd[int(float(yr[:2]))]+yr[2:]
            let = s.strip()[4:7].strip()
            num = s.strip()[7:].strip()
            if num == '':
                num = '00'
            elif len(num) == 1:
                num = '0' + num
            elif len(num) > 2:
                try:
                    num = Names.pkd[int(float(num[:-1]))]+num[-1]
                except (IndexError, ValueError):
                    raise TargetNameParseError(
                        ('{} cannot be turned into a '
                         'packed designation').format(s))
            return (yr + let[0] + num + let[1])

        else:
            raise TargetNameParseError(
                ('{} cannot be turned into a '
                 'packed number or designation').format(s))

    @staticmethod
    def from_packed(p):
        """Unpack asteroid designation/number.

        Parameters
        ----------
        p : str
           Packed target identifier.

        Returns
        -------
        s : str
           Unpacked designation/number.

        Examples
        --------
        >>> from sbpy.data import Names
        >>> Names.from_packed('J95A01A')
        '1995 AA1'
        """
        # packed number
        if p.isdigit():
            return int(p)
        elif p[0].isalpha() and p[1:].isdigit():
            return int(str(Names.pkd.find(p[0])) + p[1:])

        # old designation style, e.g.: 1989AB
        if (len(p.strip()) < 7 and p[:4].isdigit() and
                p[4:6].isalpha()):
            return p[:4]+' '+p[4:6]
        # Palomar Survey
        elif p.find("PLS") == 0:
            return p[3:] + " P-L"
        # Trojan Surveys
        elif p.find("T1S") == 0:
            return p[3:] + " T-1"
        elif p.find("T2S") == 0:
            return p[3:] + " T-2"
        elif p.find("T3S") == 0:
            return p[3:] + " T-3"
        # insert blank in designations
        elif (p[0:4].isdigit() and p[4:6].isalpha() and
              p[4] != ' '):
            return p[:4]+" "+p[4:]
        # MPC packed 7-digit designation
        elif (p[0].isalpha() and p[1:3].isdigit() and
              p[-1].isalpha() and p[-2].isdigit()):
            yr = str(Names.pkd.find(p[0]))+p[1:3]
            let = p[3]+p[-1]
            num = str(Names.pkd.find(p[4]))+p[5]
            num = num.lstrip("0")
            return yr+' '+let+num
        # nothing to do
        else:
            return p

    @staticmethod
    def parse_comet(s):
        """Parse a string as if it were a comet name.

        Considers IAU-formatted permanent and new-style
        designations. Note that comet types (P, D, C etc) are required
        and letter case is important.

        Parameters
        ----------
        s : str or list/array of str
           String, or a list/array of strings, to parse.

        Returns
        -------
        r : dict
           The dictionary contains the components identified from ``s``:
           number, orbit type, designation, name, and/or fragment. If
           none of these components are identified, a
           `TargetNameParseError` is raised.

        Raises
        ------
        TargetNameParseError : Exception
           If the string does not appear to be a comet name.

        Notes
        -----
        This function has absolutely no knowledge whether the Solar
        System small body ``s`` is an asteroid or a comet. It simply
        searches for common patterns in string ``s`` that are common for
        comet names and designations. For instance, if ``s`` contains an
        asteroid name, this function will identify that part as a
        comet name. Hence, the user is advised to take that into
        account when interpreting the parsing results.


        Examples
        --------
        >>> from sbpy.data import Names
        >>> tempel = Names.parse_comet('9P/Tempel 1')
        >>> tempel['type'], tempel['name']
        ('P', 'Tempel 1')
        >>> linear = Names.parse_comet('C/2001 A2-A (LINEAR)')
        >>> linear['desig'], linear['fragm'], linear['name']
        ('2001 A2', 'A', 'LINEAR')

        The following table shows results of the parsing:

        +------------------------------+------+----+-----+----------+------------------------+
        |targetname                    |number|type|fragm|desig     |name                    |
        +==============================+======+====+=====+==========+========================+
        |1P/Halley                     | 1    | P  |     |          |Halley                  |
        +------------------------------+------+----+-----+----------+------------------------+
        |3D/Biela                      | 3    | D  |     |          |Biela                   |
        +------------------------------+------+----+-----+----------+------------------------+
        |P/Encke                       |      | P  |     |          |Encke                   |
        +------------------------------+------+----+-----+----------+------------------------+
        |9P/Tempel 1                   | 9    | P  |     |          |Tempel 1                |
        +------------------------------+------+----+-----+----------+------------------------+
        |73P/Schwassmann Wachmann 3 C  | 73   | P  |     |          |Schwassmann Wachmann 3 C|
        +------------------------------+------+----+-----+----------+------------------------+
        |73P-C/Schwassmann Wachmann 3 C| 73   | P  | C   |          |Schwassmann Wachmann 3 C|
        +------------------------------+------+----+-----+----------+------------------------+
        |73P-BB                        | 73   | P  | BB  |          |                        |
        +------------------------------+------+----+-----+----------+------------------------+
        |322P                          | 322  | P  |     |          |                        |
        +------------------------------+------+----+-----+----------+------------------------+
        |X/1106 C1                     |      | X  |     | 1066 C1  |                        |
        +------------------------------+------+----+-----+----------+------------------------+
        |P/1994 N2 (McNaught-Hartley)  |      | P  |     |  1994 N2 |McNaught-Hartley        |
        +------------------------------+------+----+-----+----------+------------------------+
        |P/2001 YX127 (LINEAR)         |      | P  |     |2001 YX127|LINEAR                  |
        +------------------------------+------+----+-----+----------+------------------------+
        |P/2010 WK (LINEAR)            |      | P  |     | 2010 WK  |LINEAR                  |
        +------------------------------+------+----+-----+----------+------------------------+
        |C/-146 P1                     |      | C  |     | -146 P1  |                        |
        +------------------------------+------+----+-----+----------+------------------------+
        |C/2001 A2-A (LINEAR)          |      | C  | A   | 2001 A2  |LINEAR                  |
        +------------------------------+------+----+-----+----------+------------------------+
        |C/2013 US10                   |      | C  |     | 2013 US10|                        |
        +------------------------------+------+----+-----+----------+------------------------+
        |C/2015 V2 (Johnson)           |      | C  |     | 2015 V2  |Johnson                 |
        +------------------------------+------+----+-----+----------+------------------------+

        """

        import re

        # define comet matching pattern
        pat = ('^(([1-9][0-9]*[PDCX]'
               '(-[A-Z]{1,2})?)|[PDCX]/)'  # type/number/fragm [0,1,2]
               '|([-]?[0-9]{3,4}[ _][A-Z]{1,2}[0-9]{0,3}(-[1-9A-Z]{0,2})?)'
               # designation [3,4]
               '|(([dvA-Z][a-z\']? ?[A-Z]*[a-z]*[ -]?[A-Z]?[1-9]*[a-z]*)'
               '( [1-9A-Z]{1,2})*)'  # name [5,6]
               )

        # regex patterns that will be rejected
        rej_pat = ('(([1-9][0-9]*[pdcxai]\b)'  # small-caps comet number
                   '|([pdcxai]/))'  # small-caps comet type
                   )

        raw = s.translate(str.maketrans('()', '  ')).strip()

        # reject rej_pat patterns
        rej = re.findall(rej_pat, raw)

        if len(rej) > 0:
            raise TargetNameParseError('{} does not appear to be a '
                                       'comet identifier'.format(s))

        m = re.findall(pat, s)

        r = {}

        if len(m) > 0:
            for el in m:
                # type & number & fragment
                if len(el[0]) > 0:
                    typnumber = el[0].replace('/', '')
                    try:
                        r['type'] = re.findall('[PDCXAI]', typnumber)[0]
                    except IndexError:
                        pass
                    try:
                        r['number'] = int(re.findall('[0-9]*', typnumber)[0])
                    except (IndexError, ValueError):
                        pass
                    try:
                        r['fragment'] = re.findall('-[A-Z]{1,2}',
                                                   typnumber)[0][1:]
                    except IndexError:
                        pass
                # designation & fragment
                if len(el[3]) > 0:
                    r['desig'] = el[3].replace('_', ' ')
                    try:
                        r['fragm'] = re.findall('-[A-Z]{1,2}',
                                                r['desig'])[0][1:]
                        r['desig'] = r['desig'][:r['desig'].find('-' +
                                                                 r['fragm'])]
                    except IndexError:
                        pass
                # name
                if len(el[5]) > 0:
                    if len(el[5]) > 1:
                        r['name'] = el[5]

        if len(r) == 0 or 'type' not in r:
            raise TargetNameParseError(('{} does not appear to be a '
                                        'comet name').format(s))
        else:
            return r

    @staticmethod
    def parse_asteroid(s):
        """Parse a string as if it were an asteroid name.

        Considers IAU-formatted permanent and new-style designations,
        as well as MPC packed designations and numbers. Note that
        letter case is important. Parentheses are ignored in the parsing.

        Parameters
        ----------
        s : str or list of str
           The string, or a list/array of strings, to parse.

        Returns
        -------
        r : dict
           The dictionary contains the components identified from ``s``:
           IAU number, designation, and/or name. If none of these
           components are identified, a `TargetNameParseError` is raised

        Raises
        ------
        TargetNameParseError : Exception
           If the string does not appear to be an asteroid name.

        Notes
        -----
        This function has absolutely no knowledge whether the Solar
        System small body ``s`` is an asteroid or a comet. It simply
        searches for common patterns in string ``s`` that are common
        for asteroid names, numbers, or designations. For instance, if
        ``s`` contains a comet name, this function will identify that
        part as an asteroid name. Hence, the user is advised to take
        that into account when interpreting the parsing results.

        Examples
        --------
        >>> from sbpy.data import Names
        >>> ceres = Names.parse_asteroid('(1) Ceres')
        >>> ceres['number'], ceres['name']
        (1, 'Ceres')
        >>> mu = Names.parse_asteroid('2014 MU69')
        >>> mu['desig']
        '2014 MU69'

        The following table shows results of the parsing:

        +----------------------------------+----------+------+-----------------+
        |targetname                        |desig     |number|name             |
        +==================================+==========+======+=================+
        |1                                 |          |1     |                 |
        +----------------------------------+----------+------+-----------------+
        |2 Pallas                          |          |2     |Pallas           |
        +----------------------------------+----------+------+-----------------+
        |\(2001\) Einstein                 |          |2001  |Einstein         |
        +----------------------------------+----------+------+-----------------+
        |1714 Sy                           |          |1714  |Sy               |
        +----------------------------------+----------+------+-----------------+
        |2014 MU69                         |2014 MU69 |      |                 |
        +----------------------------------+----------+------+-----------------+
        |\(228195\) 6675 P-L               |6675 P-L  |228195|                 |
        +----------------------------------+----------+------+-----------------+
        |4101 T-3                          |4101 T-3  |      |                 |
        +----------------------------------+----------+------+-----------------+
        |4015 Wilson-Harrington \(1979 VA\)|1979 VA   |4015  |Wilson-Harrington|
        +----------------------------------+----------+------+-----------------+
        |J95X00A                           |1995 XA   |      |                 |
        +----------------------------------+----------+------+-----------------+
        |K07Tf8A                           |2007 TA418|      |                 |
        +----------------------------------+----------+------+-----------------+
        |G3693                             |          |163693|                 |
        +----------------------------------+----------+------+-----------------+
        |1A                                |1A        |      |                 |
        +----------------------------------+----------+------+-----------------+

        """

        import re

        pat = ('(([1A][8-9][0-9]{2}[ _][A-Z]{2}[0-9]{0,3}|'
               '20[0-9]{2}[ _][A-Z]{2}[0-9]{0,3})'
               # designation [0,1]
               '|([1-9][0-9]{3}[ _](P-L|T-[1-3])))'
               # Palomar-Leiden  [0,2,3]
               '|([IJKL][0-9]{2}[A-Z][0-9a-z][0-9][A-Z]'
               '|PLS[1-9][0-9]{3}|T1S[1-9][0-9]{3}|T2S[1-9][0-9]{3}'
               '|T3S[1-9][0-9]{3})'
               # packed desig [4]
               '|(^[A-Za-z][0-9]{4}| [A-Za-z][0-9]{4})'
               # packed number [5]
               '|([A-Z]{3,} |[A-Z]{3,}$'  # capitalized acronymns
               '|van de [A-Z][a-z]*[ ^ 0-9]*[-]?[A-Z]?[a-z]*[^0-9] *'
               '|de [A-Z][a-z]*[ ^ 0-9]*[-]?[A-Z]?[a-z]*[^0-9] *'
               "|['`]?[A-Z][A-Z]*['`]?[a-z][a-z]*['`]?[^0-9]*"
               "[ -]?[A-Z]?[a-z]*[^0-9]*)"
               # name [6]
               '|((^|\b)[1-9][0-9]*(\b|$| |_))'
               # number [7,8]
               '|^(([1-9][0-9]*A))'
               # comet-style designations: 1A [10]
               )

        # regex patterns that will be rejected
        rej_pat = ('([1-2][0-9]{0,3}[ _][A-Z][0-9]*(\b|$))'
                   # comet desig
                   '|([1-9][0-9]*[PDCXAI]\b)'
                   # comet number
                   '|([PDCXAI]/)'
                   # comet type
                   '|([1-2][0-9]{0,3}[ _][a-z]{2}[0-9]{0,3})'
                   )

        raw = s.translate(str.maketrans('()_', '   ')).strip()

        # reject rej_pat patterns
        rej = re.findall(rej_pat, raw)

        if len(rej) > 0:
            raise TargetNameParseError('{} does not appear to be an '
                                       'asteroid identifier'.format(s))

        # match target patterns
        m = re.findall(pat, raw)

        r = {}

        if len(m) > 0:
            for el in m:
                # designation
                if len(el[0]) > 0:
                    if el[0][0] == 'A':
                        r['desig'] = '1'+el[0][1:]
                    else:
                        r['desig'] = el[0]
                # packed designation
                elif len(el[4]) > 0:
                    ident = el[4]
                    r['desig'] = Names.from_packed(ident)
                # packed number
                elif len(el[5]) > 0:
                    ident = el[5]
                    r['number'] = Names.from_packed(ident)
                # number
                elif len(el[7]) > 0:
                    r['number'] = int(float(el[7]))
                # name
                elif len(el[6]) > 0:
                    if len(el[6].strip()) > 1:
                        r['name'] = el[6].strip()
                # comet-style designation
                elif len(el[10]) > 0:
                    r['desig'] = el[10].strip()

        if len(r) == 0:
            raise TargetNameParseError(('{} does not appear to be an '
                                        'asteroid name'.format(s)))
        else:
            return r

    @staticmethod
    def asteroid_or_comet(s):
        """Checks if an object identifier is more likely to belong to an
        asteroid or a comet.

        Parameters
        ----------
        s : str
           Target identifier.

        Returns
        -------
        target_type : str
           The target identification: ``'comet'`` or ``'asteroid'``.

        Notes
        -----
        This function uses the results of
        `~sbpy.data.Names.parse_asteroid` and
        `~sbpy.data.Names.parse_comet`. Hence, it is affected by
        ambiguities in the name/number/designation identification. If
        the name is ambiguous, a `~sbpy.data.names.TargetNameParseError`
        will be
        raised. Note that
        for any identifier that does not contain a comet type (P, D, C
        etc.), it is likely that the object gets identified as an
        asteroid.

        Examples
        --------
        >>> from sbpy.data import Names
        >>> print(Names.asteroid_or_comet('2P'))
        comet
        >>> print(Names.asteroid_or_comet('(1) Ceres'))
        asteroid

        """

        # compare lengths of dictionaries from parse_asteroid and
        # parse_comet; the longer one is more likely to describe the
        # nature of the target, if both dictionaries have the same
        # length, the nature is ambiguous
        ast = {}
        com = {}

        try:
            com = Names.parse_comet(s)
        except TargetNameParseError:
            pass

        try:
            ast = Names.parse_asteroid(s)
        except TargetNameParseError:
            pass

        if len(ast) > 0 and len(com) == 0:
            return 'asteroid'
        elif len(com) > 0 and len(ast) == 0:
            return 'comet'
        else:
            raise TargetNameParseError('Target nature unclear.')
