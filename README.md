# Predict Crosslinked Mass

This script is designed to assist in determining protein-protein interactions, by predicting the mass of
crosslinked peptide fragments following an enzymatic digestion of a protein complex.


## Crosslink Types

There are three types of crosslinks predicted by the program:

1. Mono-links, e.g. protein1's peptide + linker (other end of crosslinker may be hydrolysis product, etc.)
2. Interpeptide links, e.g. protein1's peptide + linker + protein1's peptide
3. Intrapeptide links, e.g. protein1's peptide + linker + protein2's peptide


## Authors
**Nick Semenkovich** | https://nick.semenkovich.com/

## License
Copyright 2014, Nick Semenkovich \<semenko@alum.mit.edu\>

Released under the MIT License. See LICENSE for details.


(Contains code adapted from Pyteomics, under the Apache license.)
