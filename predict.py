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
from collections import OrderedDict
from itertools import combinations_with_replacement, product
from pyteomics_derived import cleave, expasy_rules  # Cherry-picked & modified from the Apache-licensed Pyteomics package.
import argparse

AA_SHORT_CODES = "ARNDCEQGHILKMFPSTWYVUOBZJX"

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


def compute_crosslinked_mass(peptide_sequence, crosslinker, mode):
    """
    Given an input peptide,
    """

    return True


def main():
    # Parse & interpret command line flags.
    parser = argparse.ArgumentParser(description='Compute protein/peptide crosslinks and masses.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', dest="input", metavar='input.faa', type=str,
                        help='Input protein list (fasta format).', required=True)
    parser.add_argument('--linker', dest='linker', type=str,
                        choices=['BS3'], help='Linker to simulate.', required=True)
    parser.add_argument('--enzyme', dest='enzyme', type=str,
                        choices=expasy_rules.keys(), help='Digestion enzyme.', required=True)

    # Options for the Pyteomics-derived cleavage step
    digest_group = parser.add_argument_group('digest options')
    digest_group.add_argument('--missed-cleavages', dest="cleavages", type=int, default=0, choices=range(0,10),
                              help='Max number of missed cleavages.', required=False)
    digest_group.add_argument('--find-overlaps', dest="overlap", action='store_true',
                              help='Find overlapping cleavages [Slow!].', required=False)

    args = parser.parse_args()


    # Read & parse our input file
    faa_dict = parse_input_faa(args.input)  # KEY, PEP

    ## We need to enumerate three different types of peptide linkages:
    # 1: Silly, orphan linkages (peptide1 + linker)
    # 2: Simple self+self internal linkages (peptide1 + linker + peptide1)
    # 3: Regular, self+other linkages (peptide1 + linker + peptide2)

    # Type 1 (peptide1 + linker)
    for peptide_id, peptide_sequence in faa_dict.iteritems():
        pass

    # Type 2 (peptide1 + linker + peptide1)
    for peptide_id, peptide_sequence in faa_dict.iteritems():
        pass

    # Type 3 (peptide1 + linker + peptide2)

    for protein1, protein2, in combinations_with_replacement(faa_dict.iterkeys(), 2):
        # For each protein combination, generate digestions & the subsequent permutations
        protein1_cleavages = cleave(faa_dict[protein1], args.enzyme,
                                    missed_cleavages=args.cleavages, overlap=args.overlap)
        protein2_cleavages = cleave(faa_dict[protein2], args.enzyme,
                                    missed_cleavages=args.cleavages, overlap=args.overlap)

        print(protein1 + " " + protein2)
        for peptide1, peptide2, in product(protein1_cleavages, protein2_cleavages):
            print("\t" + peptide1 + " " + peptide2)





    # Generate all possible digestion fragments



if __name__ == '__main__':
    main()
