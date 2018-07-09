# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
======================
sbpy data.Names Module
======================

Class for dealing with object naming conventions.

created on August 28, 2017

"""

from .core import DataClass
from numpy import ndarray

__all__ = ['Names', 'TargetNameParseError']


class TargetNameParseError(Exception):
    pass


class Names():
    """Class for dealing with object naming conventions"""

    @staticmethod
    def altident(identifier, bib=None):
        """Query Lowell database to obtain alternative target names for
        `identifier`.

        Examples
        --------
        >>> from sbpy.data import Names
        >>> Names.altident('3552')  # doctest: +SKIP
        ['3552', 'Don Quixote', '1983 SA']

        not yet implemented

        """

        pass

    @staticmethod
    def to_packed(s):
        """Convert asteroid designation/number to packed identifier.

        Parameters
        ----------
        s : string
        The long target identifier.

        Returns
        -------
        p : string
        The packed designation/number.

        Examples
        --------
        >>> from sbpy.data import Names
        >>> Names.to_packed('1995 AA1')
        'J95A01A'
        """

        if isinstance(s, dict):
            ident = s
        else:
            ident = Names.parse_asteroid(s)

        # packed numbers translation string
        pkd = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghifklmnopqrstuvwxyz'

        if 'number' in ident:
            if ident['number'] < 100000:
                return ('{:05d}'.format(ident['number']))
            elif ident['number'] > 619999:
                raise TargetNameParseError(('{} cannot be turned into a '
                                            'packed number').format(ident['number']))
            else:
                mod = (ident['number'] % 10000)
                return ('{}{:04d}'.format(pkd[int((ident['number']-mod)/10000)],
                                          mod))
        elif 'desig' in ident:
            yr = ident['desig'].strip()[:4]
            yr = pkd[int(float(yr[:2]))]+yr[2:]
            let = ident['desig'].strip()[4:7].strip()
            num = ident['desig'].strip()[7:].strip()
            if len(num) == 1:
                num = '0' + num
            elif len(num) > 2:
                try:
                    num = pkd[int(float(num[:-1]))]+num[-1]
                except IndexError:
                    raise TargetNameParseError(('{} cannot be turned into a '
                                                'packed designation').format(ident['desig']))
            return (yr + let[0] + num + let[1])
        else:
            raise TargetNameParseError(('{} cannot be turned into a '
                                        'packed number or designation').format(s))

    @staticmethod
    def parse_comet(s):
        """Parse a string as if it were a comet name.

        Considers IAU-formatted permanent and new-style
        designations. Note that letter case is important.

        Parameters
        ----------
        s : string or list/array of strings
        The string, or a list/array of strings, to parse.

        Returns
        -------
        r : dict
        The dictionary contains the components identified from `s`:
        number, orbit type, designation, name, and/or fragment. If
        none of these components are identified, a
        `TargetNameParseError` is raised

        Raises
        ------
        TargetNameParseError : Exception
        If the string does not appear to be a comet name.

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

        +--------------------------------+-------+------+-------+------------+----------------------------+
        |targetname                      |number | type | fragm | desig      |  name                      |
        +================================+=======+======+=======+============+============================+
        |1P/Halley                       | 1     | 'P'  |       |            | 'Halley'                   |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |3D/Biela                        | 3     | 'D'  |       |            | 'Biela'                    |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |9P/Tempel 1                     | 9     | 'P'  |       |            | 'Tempel 1'                 |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |73P/Schwassmann Wachmann 3 C    | 73    | 'P'  |       |            | 'Schwassmann Wachmann 3 C' |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |73P-C/Schwassmann Wachmann 3 C  | 73    | 'P'  | 'C'   |            | 'Schwassmann Wachmann 3 C' |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |73P-BB                          | 73    | 'P'  | 'BB'  |            |                            |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |322P                            | 322   | 'P'  |       |            |                            |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |X/1106 C1                       |       | 'X'  |       | '1066 C1'  |                            |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |P/1994 N2 (McNaught-Hartley)    |       | 'P'  |       |  '1994 N2' | 'McNaught-Hartley'         |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |P/2001 YX127 (LINEAR)           |       | 'P'  |       |'2001 YX127'| 'LINEAR'                   |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |C/-146 P1                       |       | 'C'  |       | '-146 P1'  |                            |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |C/2001 A2-A (LINEAR)            |       | 'C'  | 'A'   | '2001 A2'  | 'LINEAR'                   |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |C/2013 US10                     |       | 'C'  |       | '2013 US10'|                            |
        +--------------------------------+-------+------+-------+------------+----------------------------+
        |C/2015 V2 (Johnson)             |       | 'C'  |       | '2015 V2'  | 'Johnson'                  |
        +--------------------------------+-------+------+-------+------------+----------------------------+

        """

        import re

        # define comet matching pattern
        pat = ('^(([1-9][0-9]*[PDCXAI]'
               '(-[A-Z]{1,2})?)|[PDCXAI]/)'  # typ/number/fragm [0,1,2]
               '|([-]?[0-9]{3,4}[ _][A-Z]{1,2}[0-9]{1,3}(-[1-9A-Z]{0,2})?)'
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

        if len(r) == 0:
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
        s : string or list/array of strings
        The string, or a list/array of strings, to parse.

        Returns
        -------
        r : dict
        The dictionary contains the components identified from `s`:
        IAU number, designation, and/or name. If none of these
        components are identfied, a `TargetNameParseError` is raised

        >>> from sbpy.data import Names
        >>> ceres = Names.parse_asteroid('(1) Ceres')
        >>> ceres['number'], ceres['name']
        (1, 'Ceres')
        >>> mu = Names.parse_asteroid('2014 MU69')
        >>> mu['desig']
        '2014 MU69'

        Examples
        --------
        The following table shows results of the parsing:

        +--------------------------------+--------------+--------+---------------------+
        |targetname                      | desig        | number | name                |
        +================================+==============+========+=====================+
        |1                               |              | 1      |                     |
        +--------------------------------+--------------+--------+---------------------+
        |2 Pallas                        |              | 2      | 'Pallas'            |
        +--------------------------------+--------------+--------+---------------------+
        |\(2001\) Einstein               |              | 2001   | 'Einstein'          |
        +--------------------------------+--------------+--------+---------------------+
        |1714 Sy                         |              | 1714   | 'Sy'                |
        +--------------------------------+--------------+--------+---------------------+
        |2014 MU69                       | '2014 MU69'  |        |                     |
        +--------------------------------+--------------+--------+---------------------+
        |(228195) 6675 P-L               | '6675 P-L'   | 228195 | None                |
        +--------------------------------+--------------+--------+---------------------+
        |4101 T-3                        | '4101 T-3'   |        |                     |
        +--------------------------------+--------------+--------+---------------------+
        |4015 Wilson-Harrington (1979 VA)| '1979 VA'    | 4015   | 'Wilson-Harrington' |
        +--------------------------------+--------------+--------+---------------------+
        |J95X00A                         | '1995 XA'    |        |                     |
        +--------------------------------+--------------+--------+---------------------+
        |K07Tf8A                         | '2007 TA418' |        |                     |
        +--------------------------------+--------------+--------+---------------------+
        |G3693                           |              | 163693 |                     |
        +--------------------------------+--------------+--------+---------------------+
        """

        import re

        # packed numbers translation string
        pkd = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghifklmnopqrstuvwxyz'

        pat = ('(([1-2][0-9]{0,3}[ _][A-Z]{2}[0-9]{0,3})'  # designation [0,1]
               '|([1-9][0-9]{3}[ _](P-L|T-[1-3])))'  # Palomar-Leiden  [0,2,3]
               '|([IJKL][0-9]{2}[A-Z][0-9a-z][0-9][A-Z])'  # packed desig [4]
               '|([A-Za-z][0-9]{4})'  # packed number [5]
               '|([A-Z][A-Z]*[a-z][a-z]*[^0-9]*'
               '[ -]?[A-Z]?[a-z]*[^0-9]*)'  # name [6]
               '|([1-9][0-9]*(\b|$| |_))')  # number [7,8]

        # regex patterns that will be rejected
        rej_pat = ('([1-2][0-9]{0,3}[ _][A-Z][0-9]*(\b|$))'  # comet desig
                   '|([1-9][0-9]*[PDCXAI]\b)'  # comet number
                   '|([PDCXAI]/)'  # comet type
                   # small-caps desig
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
                    r['desig'] = el[0]
                # packed designation (unpack here)
                elif len(el[4]) > 0:
                    ident = el[4]
                    # old designation style, e.g.: 1989AB
                    if (len(ident.strip()) < 7 and ident[:4].isdigit() and
                            ident[4:6].isalpha()):
                        r['desig'] = ident[:4]+' '+ident[4:6]
                    # Palomar Survey
                    elif ident.find("PLS") == 0:
                        r['desig'] = ident[3:] + " P-L"
                    # Trojan Surveys
                    elif ident.find("T1S") == 0:
                        r['desig'] = ident[3:] + " T-1"
                    elif ident.find("T2S") == 0:
                        r['desig'] = ident[3:] + " T-2"
                    elif ident.find("T3S") == 0:
                        r['desig'] = ident[3:] + " T-3"
                    # insert blank in designations
                    elif (ident[0:4].isdigit() and ident[4:6].isalpha() and
                          ident[4] != ' '):
                        r['desig'] = ident[:4]+" "+ident[4:]
                    # MPC packed 7-digit designation
                    elif (ident[0].isalpha() and ident[1:3].isdigit() and
                          ident[-1].isalpha() and ident[-2].isdigit()):
                        yr = str(pkd.find(ident[0]))+ident[1:3]
                        let = ident[3]+ident[-1]
                        num = str(pkd.find(ident[4]))+ident[5]
                        num = num.lstrip("0")
                        r['desig'] = yr+' '+let+num
                    # nothing to do
                    else:
                        r['desig'] = ident
                # packed number (unpack here)
                elif len(el[5]) > 0:
                    ident = el[5]
                    r['number'] = int(float(str(pkd.find(ident[0]))+ident[1:]))
                # number
                elif len(el[7]) > 0:
                    r['number'] = int(float(el[7]))
                # name (strip here)
                elif len(el[6]) > 0:
                    if len(el[6].strip()) > 1:
                        r['name'] = el[6].strip()

        if len(r) == 0:
            raise TargetNameParseError(('{} does not appear to be an '
                                        'asteroid name'.format(s)))
        else:
            return r

    @staticmethod
    def asteroid_or_comet(s):
        """Checks if an object is an asteroid, or a comet, based on its
        identifier.

        Parameters
        ----------
        s : string
        target identifier

        Returns
        -------
        target_type : string
        The target identification: 'comet', 'asteroid', or `None`.

        Examples
        --------
        >>> from sbpy.data import Names
        >>> Names.asteroid_or_comet('2P')
        'comet'
        >>> Names.asteroid_or_comet('(1) Ceres')
        'asteroid'
        >>> Names.asteroid_or_comet('Fred')


        """

        # compare lengths of dictionaries from parse_asteroid and parse_comet;
        # the longer one is more likely to describe the nature of the target,
        # if both dictionaries have the same length, the nature is ambiguous
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

        if len(ast) > len(com):
            return 'asteroid'
        elif len(com) > len(ast):
            return 'comet'
        elif ('desig' in ast and 'desig' in com and
              (ast['desig'] == com['desig'] or ast['name'] == com['name'])):
            # in this case, it's most likely to be an asteroid
            return 'asteroid'
        elif (('desig' in ast or 'number' in ast) and not
              ('desig' in com or 'number' in com)):
            return 'asteroid'
        elif (('desig' in com or 'number' in com) and not
              ('desig' in ast or 'number' in ast)):
            return 'asteroid'
        else:
            return None
