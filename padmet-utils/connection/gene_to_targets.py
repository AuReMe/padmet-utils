# -*- coding: utf-8 -*-
"""
Description:
    From a list of genes, get from the linked reactions the list of products.

    R1 is linked to G1, R1 produces M1 and M2.  output: M1,M2. Takes into account reversibility

::

    usage:
        gene_to_targets.py --padmetSpec=FILE --genes=FILE --output=FILE [-v]
    
    option:
        -h --help     Show help
        --padmetSpec=FILE    path to the padmet file
        --genes=FILE   path to the file containing gene ids, one id by line
        --output=FILE    path to the output file containing all tagerts which can by produced by all reactions associated to the given genes
        -v   print info
"""
from padmet.classes import PadmetSpec
from padmet.utils.connection import gene_to_targets
import docopt

def main():
    #recovering args
    args = docopt.docopt(__doc__)
    padmet = PadmetSpec(args["--padmetSpec"])
    genes_file = args["--genes"]
    output = args["--output"]
    verbose = args["-v"]
    gene_to_targets.gene_to_targets(padmet, genes_file, output, verbose)
    
        
                
if __name__ == "__main__":
    main()