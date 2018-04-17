# -*- coding: utf-8 -*-
"""
This file is part of padmet-utils.

padmet-utils is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

padmet-utils is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with padmet-utils. If not, see <http://www.gnu.org/licenses/>.

@author: Meziane AITE, meziane.aite@inria.fr
Description:
From a given sbml file, create a sbml with only the reactions associated to a gene.
Need for a reaction, in section 'note', 'GENE_ASSOCIATION': ....

usage:
    extract_rxn_with_gene_assoc.py --sbml=FILE --output=FILE [-v]

options:
    -h --help     Show help.
    --sbml=FILE    pathname to the sbml file
    --output=FILE    pathname to the sbml output (with only rxn with genes assoc)
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
        if "GENE_ASSOCIATION" not in parseNotes(reaction).keys():
            reactions_to_remove.append(reaction.getId())
    for rId in reactions_to_remove:
        listOfReactions.remove(rId)
    
    libsbml.writeSBMLToFile(document, output)

if __name__ == "__main__":
    main()