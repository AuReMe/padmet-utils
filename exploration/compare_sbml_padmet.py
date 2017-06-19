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