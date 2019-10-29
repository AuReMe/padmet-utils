# -*- coding: utf-8 -*-
"""
Description:
    extract 1 reaction (if rxn_id) or a list of reactions (if rxn_file) 
    from a sbml file to the form used in aureme for curation.
    For example use this script to extract specific missing reaction of a model to
    a just created metabolic network.

::

    usage:
        sbml_to_curation_form.py --sbml=FILE --output=FILE --rxn_id=ID [--comment=STR] [--extract-gene] [-v]
        sbml_to_curation_form.py --sbml=FILE --output=FILE --rxn_file=FILE [--comment=STR] [--extract-gene] [-v]
    
    options:
        -h --help     Show help.
        --sbml=FILE    path of the sbml.
        --output=FILE    form containing the reaction extracted, form used for manual curation in aureme.
        --rxn_id=FILE    id of one reaction to extract
        --rxn_file=FILE    file of reactions ids to extract, 1 id by line.
        --extract-gene    If true, extract also genes associated to reactions.
        --comment=STR    comment associated to the reactions in the form. Used to track sources of curation in aureme [default: "N.A"].
        -v   print info

"""

import docopt
from padmet.utils.connection import sbml_to_curation_form

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    if args["--rxn_id"]:
        rxn_list = [args["--rxn_id"]]
    if args["--rxn_file"]:
        with open(args["--rxn_id"], 'r') as f:
            rxn_list = f.read.splitlines()
    output = args["--output"]
    comment = args["--comment"]
    verbose = args["-v"]
    extract_gene = args["--extract-gene"]
    sbml_to_curation_form.sbml_to_curation(sbml_file, rxn_list, output, extract_gene=extract_gene, comment=comment, verbose=verbose)

if __name__ == "__main__":
    main()