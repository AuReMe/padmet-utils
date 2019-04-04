# -*- coding: utf-8 -*-
"""
Description:
    compare reactions in sbml and padmet file

::
    
    usage:
        compare_sbml_padmet.py --padmet=FILE --sbml=FILE
    
    option:
        -h --help    Show help.
        --padmet=FILE    path of the padmet file
        --sbml=FILE    path of the sbml file
"""
from padmet.classes import PadmetSpec
import padmet.utils.sbmlPlugin as sp
import libsbml
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    sbml_model = document.getModel()

    padmet = PadmetSpec(args["--padmet"])
    compare_sbml_padmet(sbml_model, padmet)
    
def compare_sbml_padmet(sbml_model, padmet):
    """
    compare reactions ids in sbml vs padmet, return nb of reactions in both
    and reactions id not in sbml or not in padmet

    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to udpate
    sbml_file: libsbml.document.model
        sbml model

    """

    sbml_listOfReactions = set([sp.convert_from_coded_id(r.getId())[0] for r in sbml_model.getListOfReactions()])
    padmet_reaction = set([node.id for node in padmet.dicOfNode.values() if node.type == "reaction"])
    diff = sbml_listOfReactions.difference(padmet_reaction)
    diff_inv = padmet_reaction.difference(sbml_listOfReactions)
    
    print("%s reaction in sbml" %len(sbml_listOfReactions))
    print("%s reaction in padmet" %len(padmet_reaction))
    print("%s reaction in sbml not in padmet" %len(diff))
    for i in diff:
        print("\t%s" %i)
    print("%s reaction in padmet not in sbml" %len(diff_inv))
    for i in diff_inv:
        print("\t%s" %i)

if __name__ == "__main__":
    main()