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
From a file containing a list of reaction, return the pathways where this reactions 
are involved.
ex: if rxn-a in pwy-x => return, pwy-x; all rxn ids in pwy-x; all rxn ids in pwy-x FROM the list; ratio

usage:
    get_pwy_from_rxn.py --reaction_file=FILE --padmetRef=FILE  --output=FILE

options:
    -h --help     Show help.
    --reaction_file=FILE    pathname of the file containing the reactions id, 1/line 
    --padmetRef=FILE    pathname of the padmet representing the database.
    --output=FILE    pathname of the file with line = pathway id, all reactions id, reactions ids from reaction file, ratio. sep = "\t"
"""
from padmet.padmetSpec import PadmetSpec
import docopt

def main():
    args = docopt.docopt(__doc__)
    reaction_file = args["--reaction_file"]
    padmet_file = args["--padmetRef"]
    output = args["--output"]
    padmet = PadmetSpec(padmet_file)
    
    dict_pwy = {}
    pwy_rlts = [rlt for rlt in padmet.getAllRelation() if rlt.type == "is_in_pathway"]
    for rlt in pwy_rlts:
        try:
            dict_pwy[rlt.id_out].add(rlt.id_in)
        except KeyError:
            dict_pwy[rlt.id_out] = set([rlt.id_in])
    
    with open(reaction_file, 'r') as f:
        rxn_to_check = set(f.read().splitlines())
    with open(output, 'w') as f:
        header = ["pathway_id","total_rxn","rxn_from_list","ration"]
        f.write("\t".join(header)+"\n")
        for pwy_id,all_rxn in dict_pwy.iteritems():
            rxn_in_pwy = all_rxn.intersection(rxn_to_check)
            if len(rxn_in_pwy) != 0:
                ratio = round(float(len(rxn_in_pwy))/float(len(all_rxn)),2)
                line = [pwy_id, ";".join(all_rxn), ";".join(rxn_in_pwy),str(ratio)]
                f.write("\t".join(line)+"\n")
    
if __name__ == "__main__":
    main()    
