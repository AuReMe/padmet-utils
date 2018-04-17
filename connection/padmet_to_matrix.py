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
Create a stoichiometry matrix from a padmet file.
The columns represent the reactions and rows represent metabolites.
S[i,j] contains the quantity of metabolite 'i' produced (negative for consumed)
by reaction 'j'.

usage:
    padmet_to_matrix.py --padmet=FILE --output=FILE

option:
    -h --help    Show help.
    --padmet=FILE    pathname to the padmet file to convert.
    --output=FILE    pathname to the output file, col: rxn, row: metabo, sep = "\t".
"""
import docopt
from padmet.classes import PadmetSpec

def main():
    args = docopt.docopt(__doc__) 
    padmet_file = args["--padmet"]
    output = args["--output"]
    padmet = PadmetSpec(padmet_file)

    all_metabolites = sorted(set([rlt.id_out for rlt in padmet.getAllRelation() if rlt.type in ["consumes","produces"]]))
    all_reactions = sorted([node.id for node in padmet.dicOfNode.values() if node.type == "reaction"]) 
    #col = reactions, row = metabolites
    matrix_dict = {}
    for met in all_metabolites:
        matrix_dict[met] = {}
        all_cp_rlts = [rlt for rlt in padmet.dicOfRelationOut[met] if rlt.type in ["consumes","produces"]]
        for rlt in all_cp_rlts:
            try:
                stoich = int(rlt.misc["STOICHIOMETRY"][0])
            except ValueError:
                stoich = 1
            if rlt.type == "consumes":
                stoich = stoich * -1
            matrix_dict[met][rlt.id_in] = stoich
            
    with open(output, 'w') as f:
        header = ["metabolites-reactions"] + all_reactions
        f.write("\t".join(header)+"\n")
        for met_id in all_metabolites:
            line = [met_id]
            for rxn_id in all_reactions:
                try:
                    line.append(str(matrix_dict[met_id][rxn_id]))
                except KeyError:
                    line.append("0")
            f.write("\t".join(line)+"\n")
                
            
if __name__ == "__main__":
    main()        
    
