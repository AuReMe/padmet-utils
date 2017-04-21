# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
compare reactions in sbml and padmet files

usage:
    compare_sbml_padmet.py --padmet=FILE --sbml=FILE

option:
    -h --help    Show help.
    --padmet=FILE    pathname of the padmet file
    --sbml=FILE    pathanme of the sbml file
"""
from padmet.padmetSpec import PadmetSpec
import padmet.sbmlPlugin as sp
from libsbml import *
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    padmetSpec = PadmetSpec(args["--padmet"])
    
    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    model = document.getModel()
    sbml_listOfReactions = set([sp.convert_from_coded_id(r.getId()) for r in model.getListOfReactions()])
    
    padmet_reaction = set([node.id for node in padmetSpec.dicOfNode.itervalues() if node.type == "reaction"])
    
    diff = sbml_listOfReactions.difference(padmet_reaction)
    diff_inv = padmet_reaction.difference(sbml_listOfReactions)
    
    print(len(diff))
    print(len(diff_inv))

if __name__ == "__main__":
    main()