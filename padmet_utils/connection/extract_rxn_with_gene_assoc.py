# -*- coding: utf-8 -*-
"""
Description:
    From a given sbml file, create a sbml with only the reactions associated to a gene.

    Need for a reaction, in section 'note', 'GENE_ASSOCIATION': ....

::

    usage:
        extract_rxn_with_gene_assoc.py --sbml=FILE --output=FILE [-v]
    
    options:
        -h --help     Show help.
        --sbml=FILE    path to the sbml file
        --output=FILE    path to the sbml output (with only rxn with genes assoc)
        -v   print info
"""
from padmet.utils.connection import extract_rxn_with_gene_assoc
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml = args["--sbml"]
    output = args["--output"]
    verbose = args["-v"]

    extract_rxn_with_gene_assoc.extract_rxn_with_gene_assoc(sbml, output, verbose)


if __name__ == "__main__":
    main()