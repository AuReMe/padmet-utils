#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Convert PADMet to ASP following these predicats:
    common_name({reaction_id or enzyme_id or pathway_id or compound_id} , common_name)
    direction(reaction_id, reaction_direction). reaction_direction in[LEFT-TO-RIGHT,REVERSIBLE]
    ec_number(reaction_id, ec(x,x,x)).
    catalysed_by(reaction_id, enzyme_id).
    uniprotID(enzyme_id, uniprot_id). #if has has_xref and db = "UNIPROT"
    in_pathway(reaction_id, pathway_id).
    reactant(reaction_id, compound_id, stoechio_value).
    product(reaction_id, compound_id, stoechio_value).
    is_a(compound_id, class_id).
    is_a(pathway_id, pathway_id).

::
    
    usage:
        padmet_to_asp.py --padmet=FILE --output=FILE [-v]
    
    option:
        -h --help     Show help.
        --padmet=FILE    path to padmet file to convert.
        --output=FILE    path to output file in lp format.
        -v    print info.
    
"""
from padmet.utils.connection import padmet_to_asp
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmet_file = args["--padmet"]
    output = args["--output"]
    verbose = args["-v"]

    padmet_to_asp.padmet_to_asp(padmet_file, output, verbose)

if __name__ == "__main__":
    main()
