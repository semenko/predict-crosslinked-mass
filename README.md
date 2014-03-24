# Predict Crosslinked Mass

This script is designed to assist in determining protein-protein interactions, by predicting the mass of
crosslinked peptide fragments following an enzymatic digestion of a protein complex.


## Crosslinkers

There are three types of crosslinks predicted by the program:

1. Mono-links, e.g. protein1's peptide + linker (other end of crosslinker may be hydrolysis product, etc.)
2. Interpeptide links (e.g. protein1's peptide + linker + protein1's peptide)
3. Intrapeptide links (e.g. protein1's peptide + linker + protein2's peptide)

### Crosslink Support

Currently this supports:

* BS3 bis(sulfosuccinimidyl)suberate

Other linkers can be easily defined as functions in ```crosslinkers.py```

## Usage

    usage: predict.py [options]
    
    Compute protein/peptide crosslinks and masses.
    
    optional arguments:
      -h, --help            show this help message and exit
      --input input.faa     Input protein list (fasta format). (default: None)
      --linker {bs3}        Linker to simulate. (default: None)
      --enzyme {clostripain,chymotrypsin-high-specificity,staphylococcal-peptidase-i,granzyme-b,iodosobezoic-acid,hydroxylamine,proline-endopeptidase,enterokinase,cnbr,arg-c,pepsin-ph1.3,pepsin-ph2.0,glutamyl-endopeptidase,bnps-skatole,caspase6,caspase7,caspase4,caspase5,caspase2,caspase3,caspase1,thermolysin,thrombin,proteinase-k,lysc,asp-n,ntcb,caspase8,caspase9,factor-xa,formic-acid,caspase10,trypsin,chymotrypsin-low-specificity}
                            Digestion enzyme. (default: None)
    
    digest options:
      --missed-cleavages {0,1,2,3,4,5,6,7,8,9}
                            Max number of missed cleavages. (default: 0)
      --find-overlaps       Find overlapping cleavages [SLOW!]. (default: False)
    
    other options:
      --ox-met              Oxidize all methionines. (default: False)

    Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at
    Washington University in St. Louis: http://gordonlab.wustl.edu.


## Authors
**Nick Semenkovich** | https://nick.semenkovich.com/

## License
Copyright 2014, Nick Semenkovich \<semenko@alum.mit.edu\>

Released under the MIT License. See LICENSE for details.


(Contains code adapted from Pyteomics, under the Apache license.)
