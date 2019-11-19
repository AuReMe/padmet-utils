# -*- coding: utf-8 -*-
"""
Description:
    Create a stoichiometry matrix from a padmet file.

    The columns represent the reactions and rows represent metabolites.

    S[i,j] contains the quantity of metabolite 'i' produced (negative for consumed)
    by reaction 'j'.

::

    usage:
        padmet_to_matrix.py --padmet=FILE --output=FILE
    
    option:
        -h --help    Show help.
        --padmet=FILE    path to the padmet file to convert.
        --output=FILE    path to the output file, col: rxn, row: metabo, sep = "\t".
"""
import docopt
from padmet.classes import PadmetSpec
from padmet.utils.connection import padmet_to_matrix

def main():
    args = docopt.docopt(__doc__) 
    padmet_file = args["--padmet"]
    output = args["--output"]
    padmet = PadmetSpec(padmet_file)
    padmet_to_matrix.padmet_to_matrix(padmet, output)

if __name__ == "__main__":
    main()        
    
