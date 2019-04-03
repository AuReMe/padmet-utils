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
import libsbml
from padmet.utils.sbmlPlugin import parseNotes
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    output = args["--output"]
    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_file)
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    
    reactions_to_remove = []
    for reaction in listOfReactions:
        if "GENE_ASSOCIATION" not in list(parseNotes(reaction).keys()):
            reactions_to_remove.append(reaction.getId())
    for rId in reactions_to_remove:
        listOfReactions.remove(rId)
    
    libsbml.writeSBMLToFile(document, output)

if __name__ == "__main__":
    main()