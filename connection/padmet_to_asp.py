# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Convert PADMet to ASP following this predicats:
direction(reaction_id, reaction_direction). #LEFT-TO-RIGHT or REVERSIBLE
ec_number(reaction_id, ec(x,x,x)).
catalysed_by(reaction_id, enzyme_id).
uniprotID(enzyme_id, uniprot_id). #if has has_xref and db = "UNIPROT"
common_name({reaction_id or enzyme_id or pathway_id or compound_id} , common_name)

usage:
    padmet_to_ASP.py --padmet=FILE --output=FILE [-v]

options:
    -h --help     Show help.
    --padmet=FILE    pathname of the padmet file to convert
    --output=FILE    pathname of the lp file
    -v   print info
"""
from padmet.aspGenerator import padmet_to_asp
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmet_file = args["--padmet"]
    output = args["--output"]
    verbose = args["-v"]

    padmet_to_asp(padmet_file, output, verbose)

if __name__ == "__main__":
    main()