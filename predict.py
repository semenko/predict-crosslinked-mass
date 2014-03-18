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
import argparse


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
    return True

def compute_crosslinked_mass(peptide_sequence, crosslinker, mode):
    """
    Given an input peptide,
    """

    return True


def main():
    parser = argparse.ArgumentParser(description='Compute protein/peptide crosslinks and masses.')
    parser.add_argument('--in_faa', metavar='input.faa', type=str,
                        help='Input protein/peptide list (fasta format).', required=True)
    parser.add_argument('--linker', dest='linker', type=str,
                        choices=['BS3'], help='Linker to simulate.', required=True)
    parser.add_argument('--enzyme', dest='enzyme', type=str,
                        choices=['trypsin'], help='Digestion enzyme.', required=True)

    args = parser.parse_args()
    print(args.accumulate(args.integers))

    pass


if __name__ == '__main__':
    main()
