#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Nick Semenkovich <semenko@alum.mit.edu> / WUSTL
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
# http://gordonlab.wustl.edu/
#
# This software is released under the MIT License:
#  <http://www.opensource.org/licenses/mit-license.php>
#
# Source: https://github.com/semenko/predict-crosslinked-mass

from __future__ import absolute_import, division, print_function, unicode_literals
from collections import deque, OrderedDict
from itertools import combinations
import argparse
from pyteomics_derived import cleave  # Cherry-picked & modified from the Apache-licensed Pyteomics package.

AA_SHORT_CODES = "ARNDCEQGHILKMFPSTWYVUOBZJX"

# Digestion options
# Definitions are: (["MOTIFS", "MOTIFS"], [CUT_INDEX, CUT_INDEX])
ENZYME_RULES = {
    "trypsin": (["K","R"], [1,1])
}

# Crosslinker options
CROSSLINKER_RULES = {
    "BF3": True
}

def parse_input_faa(in_faa):
    """
    Read & store an input FAA in memory.

    Input:
        in_fna = path_to_input.faa (string)
    Output:
        PASS
    """
    faa_dict = OrderedDict()
    working_key = None
    with open(in_faa) as f:
        for line in f:
            if line.startswith(">"):
                working_key = line[1:].strip()
                working_peptide = []
            else:
                working_peptide.append(line.strip())
            faa_dict[working_key] = ''.join(working_peptide)

    # Be paranoid and make sure we understand the AA codes we were given
    for peptide in faa_dict.itervalues():
        assert(all([aa in AA_SHORT_CODES for aa in peptide]))
    return faa_dict


def predict_digestion(peptide_sequence, enzyme):
    """
    Given an input peptide sequence & digestion enzyme, return digests.

    Input:
        peptide_sequence = "TWNTGIMLLLITMATAFMGYVLPWGQMSFWGA" (string)
        enzyme = ""
    Output:
        PASS (all possible digested fragments)
    """
    enzyme_rule = ENZYME_RULES[enzyme]
    for motif, cutindex in enzyme_rule:
        split_list = []
        # Because we can cut at different areas after/before motifs, we can't use .split()

        print(peptide_sequence)
        if str(motif) in peptide_sequence:
            pass
        print(motif)
        print(cutindex)
    return True



def compute_crosslinked_mass(peptide_sequence, crosslinker, mode):
    """
    Given an input peptide,
    """

    return True


def main():
    parser = argparse.ArgumentParser(description='Compute protein/peptide crosslinks and masses.')
    parser.add_argument('--input', dest="input", metavar='input.faa', type=str,
                        help='Input protein/peptide list (fasta format).', required=True)
    parser.add_argument('--linker', dest='linker', type=str,
                        choices=['BS3'], help='Linker to simulate.', required=True)
    parser.add_argument('--enzyme', dest='enzyme', type=str,
                        choices=['trypsin'], help='Digestion enzyme.', required=True)

    args = parser.parse_args()

    # First, read our input file
    faa_dict = parse_input_faa(args.input)

    # Cleave
    print(cleave('AKAKBKCK', 'proteinase-k', missed_cleavages=2))

    for peptide_id, peptide_sequence in faa_dict.iteritems():
        #print(peptide_id)
        # print(peptide_sequence)
       #  digest = predict_digestion(peptide_sequence, args.enzyme)

        break


    # Generate all possible digestion fragments



if __name__ == '__main__':
    main()
