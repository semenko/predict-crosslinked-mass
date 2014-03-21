"""
Code adapted from the awesome Pyteomics package.

Code:
    http://hg.theorchromo.ru/pyteomics/overview
    https://pypi.python.org/pypi/pyteomics
    http://pythonhosted.org/pyteomics/


Original Pyteomics copyright:

Copyright 2012 Anton Goloborodko, Lev Levitsky

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from itertools import chain
from collections import deque
import re

"""
This dict contains regular expressions for cleavage rules of the most
popular proteolytic enzymes. The rules were taken from the
`PeptideCutter tool
<http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
at Expasy.

These describe the C-terminal sites of cleavage.
"""
expasy_rules = {
    'arg-c':         'R',
    'asp-n':         '\w(?=D)',
    'bnps-skatole' : 'W',
    'caspase1':     '(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase2':     '(?<=DVA)D(?=[^PEDQKR])',
    'caspase3':     '(?<=DMQ)D(?=[^PEDQKR])',
    'caspase4':     '(?<=LEV)D(?=[^PEDQKR])',
    'caspase5':     '(?<=[LW]EH)D',
    'caspase6':     '(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase7':     '(?<=DEV)D(?=[^PEDQKR])',
    'caspase8':     '(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase9':     '(?<=LEH)D',
    'caspase10':    '(?<=IEA)D',
    'chymotrypsin-low-specificity' : '([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin-high-specificity':
        '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain':   'R',
    'cnbr':          'M',
    'enterokinase':  '(?<=[DN][DN][DN])K',
    'factor-xa':     '(?<=[AFGILTVM][DE]G)R',
    'formic-acid':   'D',
    'glutamyl-endopeptidase': 'E',
    'granzyme-b':    '(?<=IEP)D',
    'hydroxylamine': 'N(?=G)',
    'iodosobezoic-acid': 'W',
    'lysc':          'K',
    'ntcb':          '\w(?=C)',
    'pepsin-ph1.3':  '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                     '((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'pepsin-ph2.0':  '((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                     '((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'proline-endopeptidase': '(?<=[HKR])P(?=[^P])',
    'proteinase-k':  '[AEFILTVWY]',
    'staphylococcal-peptidase-i': '(?<=[^E])E',
    'thermolysin':   '[^DE](?=[AFILMV])',
    'thrombin':      '((?<=G)R(?=G))|'
                     '((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin':       '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
    }


def cleave(sequence, rule, missed_cleavages, overlap):
    """Cleaves a polypeptide sequence using a given rule.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    rule : str
        The name of a digestion enzyme in expasy_rules.
    missed_cleavages : int, optional
        The maximal number of allowed missed cleavages. Defaults to 0.
    overlap : bool, optional
        Set this to :py:const:`True` if the cleavage rule is complex and
        it is important to get all possible peptides when the matching
        subsequences overlap (e.g. 'XX' produces overlapping matches when
        the sequence contains 'XXX'). Default is :py:const:`False`.
        Use with caution: enabling this results in exponentially growing
        execution time.

    Returns
    -------
    out : set
        A set of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0, False)
    set(['AK', 'BK'])
    >>> cleave('AKAKBKCK', expasy_rules['trypsin'], 2, False)
    set(['CK', 'AKBK', 'BKCK', 'AKAK', 'AKBKCK', 'AK', 'AKAKBK', 'BK'])

    """
    rule_regex = expasy_rules[rule]
    peptides = set()
    cleavage_sites = deque([0], maxlen=missed_cleavages+2)
    for i in chain(map(lambda x: x.end(), re.finditer(rule_regex, sequence)),
                   [None]):
        cleavage_sites.append(i)
        for j in range(0, len(cleavage_sites)-1):
            peptides.add(sequence[cleavage_sites[j]:cleavage_sites[-1]])
        if overlap and i not in {0, None}:
            peptides.update(
                    cleave(sequence[i:], rule, missed_cleavages, overlap))

    if '' in peptides:
        peptides.remove('')
    return peptides
