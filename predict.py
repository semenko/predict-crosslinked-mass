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
from pyteomics_derived import cleave, expasy_rules  # Adapted from the Apache-licensed Pyteomics package.
from memoization import memoize_single, memoize_args
import argparse

# Wrap external cleave function in memoize decorator.
cleave = memoize_args(cleave)


# http://physics.nist.gov/cuu/Constants/
MASS_PROTON = 1.007276466812

# For sanity-checking .faa input proteins
AA_SHORT_CODES = "ACEDGFIHKMLNQPSRTWVY"

AA_AVERAGE_MASSES = {
    'A':  71.0788, 'C': 103.1388, 'E': 129.1155, 'D': 115.0886,
    'G':  57.0519, 'F': 147.1766, 'I': 113.1594, 'H': 137.1411,
    'K': 128.1741, 'M': 131.1926, 'L': 113.1594, 'N': 114.1038,
    'Q': 128.1307, 'P':  97.1167, 'S':  87.0782, 'R': 156.1875,
    'T': 101.1051, 'W': 186.2132, 'V':  99.1326, 'Y': 163.1760
}

# Crosslinker definitions
#   Rule: ((, mass_shift_when_singly_linked, mass_shift_when_doubly_linked])
#
CROSSLINKER_RULES = {
    "BS3": ((("K"), ("K")),  [156.08, 138.07])
}


@memoize_single
def mass(peptide_sequence):
    """ Compute a given peptide's mass. """
    return sum([AA_AVERAGE_MASSES[aa] for aa in peptide_sequence])


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
    for protein_id, protein_sequence in faa_dict.iteritems():
        assert("," not in protein_id)  # Unescaped commas would break our CSV output.
        assert(all([aa in AA_SHORT_CODES for aa in protein_sequence]))
    return faa_dict


@memoize_args
def can_crosslink(crosslinker, peptide1, peptide2):
    # Ignore peptides of length 2 or less
    if len(peptide1) > 2 and len(peptide2) > 2:
        if "K" in peptide1 and "K" in peptide2:  ## TODO: fix this hack
            return True
    return False


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
    digest_group.add_argument('--missed-cleavages', dest="cleavages", type=int, default=0, choices=range(0, 10),
                              help='Max number of missed cleavages.', required=False)
    digest_group.add_argument('--find-overlaps', dest="overlap", action='store_true',
                              help='Find overlapping cleavages [Slow!].', required=False)

    # Other various options
    other_group = parser.add_argument_group('other options')
    other_group.add_argument('--ox-met', dest="ox_met", action='store_true',
                             help='Oxidize all methionines.', required=False)

    args = parser.parse_args()

    if args.ox_met is True:
        AA_AVERAGE_MASSES['M'] += 15.9949  # Oxidized methionine

    # Read & parse our input file
    faa_sequence_dict = parse_input_faa(args.input)

    print("type,protein1_id,peptide1,protein2_id,peptide2,mass")

    ## We need to enumerate three different types of peptide linkages:
    # 1: Mono-links, e.g. protein1's peptide + linker (other end of crosslinker may be hydrolysis product, etc.)
    # 2: Interpeptide links, e.g. protein1's peptide + linker + protein1's peptide
    # 3: Intrapeptide links, e.g. protein1's peptide + linker + protein2's peptide

    # Types 1 & 2 (linker associated with one protein's peptides only)
    for protein_id, protein_sequence in faa_sequence_dict.iteritems():
        protein_cleavages = cleave(faa_sequence_dict[protein_id], args.enzyme, args.cleavages, args.overlap)
        # Type 1, mono-links
        for peptide in protein_cleavages:
            if can_crosslink(args.linker, peptide, peptide):
                print(1, protein_id, peptide, protein_id, peptide, mass(peptide) + 156.07864 + MASS_PROTON, sep=",")

        # Type 2, interpeptide links
        for peptide1, peptide2 in combinations_with_replacement(protein_cleavages, 2):
            if can_crosslink(args.linker, peptide1, peptide2):
                print(2, protein_id, peptide1, protein_id, peptide2, mass(peptide1) + mass(peptide2) + 138.06808 + MASS_PROTON, sep=",")


    # Type 3 (intrapeptide combinations)
    for protein1_id, protein2_id, in combinations_with_replacement(faa_sequence_dict.iterkeys(), 2):
        # For each protein combination, generate digestions & the subsequent permutations
        protein1_cleavages = cleave(faa_sequence_dict[protein1_id], args.enzyme, args.cleavages, args.overlap)
        protein2_cleavages = cleave(faa_sequence_dict[protein2_id], args.enzyme, args.cleavages, args.overlap)

        # print(protein1_id + " " + protein2_id)
        for peptide1, peptide2, in product(protein1_cleavages, protein2_cleavages):
            print(3, protein1_id, peptide1, protein2_id, peptide2, mass(peptide1) + mass(peptide2) + 138.06808 + MASS_PROTON, sep=",")


if __name__ == '__main__':
    main()
