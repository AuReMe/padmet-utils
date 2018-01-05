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
from cobra import *
from cobra.io.sbml import create_cobra_model_from_sbml_file
import docopt
"""

"""
def main():
    args = docopt.docopt(__doc__)
    sbml_1 = create_cobra_model_from_sbml_file(args["--sbml1"])
    sbml_2 = create_cobra_model_from_sbml_file(args["--sbml2"])
    
    print("sbml1:")
    print("metabolites: %s" %(len(sbml_1.metabolites)))
    print("reactions: %s" %(len(sbml_1.reactions)))
    print("sbml2:")
    print("metabolites: %s" %(len(sbml_2.metabolites)))
    print("reactions: %s" %(len(sbml_2.reactions)))
    print("reactions not in sbml1: %s" %len([i for i in sbml_2.reactions if i not in sbml_1.reactions]))
    print("reactions not in sbml2: %s" %len([i for i in sbml_1.reactions if i not in sbml_2.reactions]))
    
    all_diff = set()
    for rxn1 in sbml_1.reactions:
        rxn_id = rxn1.id
        try:
            rxn2 = sbml_2.reactions.get_by_id(rxn_id)
            if not compare_rxn(rxn1, rxn2):
                print("%s is diff" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass
                
    for rxn2 in sbml_2.reactions:
        rxn_id = rxn2.id
        try:
            rxn1 = sbml_1.reactions.get_by_id(rxn_id)
            if not compare_rxn(rxn1, rxn2) and rxn_id not in all_diff:
                print("%s is diff" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass

def compare_rxn(rxn1,rxn2):
    same_cpd, same_rev = False, False
    rxn1_dict,rxn2_dict = {}, {}
    for k,v in rxn1.metabolites.items():
        rxn1_dict[k.id] = v
    for k,v in rxn2.metabolites.items():
        rxn2_dict[k.id] = v
    if rxn1_dict == rxn2_dict: 
        same_cpd = True
    if rxn1.reversibility == rxn2.reversibility:

       same_rev = True

    return same_cpd and same_rev
        
        
    
    
if __name__ == "__main__":
    main()