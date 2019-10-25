#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    The module sbmlGenerator contains functions to generate sbml files from padmet and txt
    usign the libsbml package

::

    usage:
        sbmlGenerator.py --padmet=FILE --output=FILE --sbml_lvl=STR [--model_id=STR] [--obj_fct=STR] [--mnx_chem_prop=FILE] [--mnx_chem_xref=FILE] [-v]
        sbmlGenerator.py --padmet=FILE --output=FILE [--init_source=STR] [-v]
        sbmlGenerator.py --compound=FILE --output=FILE [--padmetRef=FILE] [-v]
        sbmlGenerator.py --reaction=FILE --output=FILE --padmetRef=FILE [-v]
    
    option:
        -h --help    Show help.
        --padmet=FILE    path of the padmet file to convert into sbml
        --output=FILE    path of the sbml file to generate.
        --mnx_chem_prop=FILE    path of the MNX chemical compounds properties.
        --mnx_chem_xref=FILE    path of the mnx dict of chemical compounds id mapping.
        --reaction=FILE    path of file of reactions ids, one by line to convert to sbml.
        --compound=FILE    path of file of compounds ids, one by line to convert to sbml.
        --init_source=STR    Select the reactions of padmet to convert on sbml based on the source of the reactions, check relations rxn has_reconstructionData.
        --sbml_lvl=STR    sbml level of output. [default 3]
        --obj_fct=STR    id of the reaction objective.
        -v   print info.
"""
from padmet.utils.connection import sbmlGenerator
import docopt

def main():
    args = docopt.docopt(__doc__)
    output = args["--output"]
    obj_fct = args["--obj_fct"]
    mnx_chem_xref = args["--mnx_chem_xref"]
    mnx_chem_prop = args["--mnx_chem_prop"]
    sbml_lvl = args["--sbml_lvl"]
    model_id = args["--model_id"]
    verbose = args["-v"]

    if args["--padmet"]:
        padmet_file = args["--padmet"]
        if args["--init_source"]:
            init_source = args["--init_source"]
            sbmlGenerator.from_init_source(padmet_file, init_source, output, verbose)
        else:
            sbmlGenerator.padmet_to_sbml(padmet_file, output, model_id, obj_fct, sbml_lvl, mnx_chem_prop, mnx_chem_xref, verbose)
    elif args["--reaction"]:
        padmetRef = args["--padmetRef"]
        reactions = args["--reaction"]
        sbmlGenerator.reaction_to_sbml(reactions, output, padmetRef, verbose)
    elif args["--compound"]:
        species_compart = args["--compound"]
        sbmlGenerator.compound_to_sbml(species_compart, output, verbose)


if __name__ == "__main__":
    main()