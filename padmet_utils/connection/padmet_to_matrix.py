# -*- coding: utf-8 -*-
"""
Description:
    Create a stoichiometry matrix from a padmet file.

    The columns represent the reactions and rows represent metabolites.

    S[i,j] contains the quantity of metabolite 'i' produced (negative for consumed)
    by reaction 'j'.

::

    usage:
        padmet_to_matrix.py --padmet=FILE --output=FILE
    
    option:
        -h --help    Show help.
        --padmet=FILE    path to the padmet file to convert.
        --output=FILE    path to the output file, col: rxn, row: metabo, sep = "\t".
"""
import docopt
from padmet.classes import PadmetSpec

def main():
    args = docopt.docopt(__doc__) 
    padmet_file = args["--padmet"]
    output = args["--output"]
    padmet = PadmetSpec(padmet_file)

    all_metabolites = sorted(set([rlt.id_out for rlt in padmet.getAllRelation() if rlt.type in ["consumes","produces"]]))
    all_reactions = sorted([node.id for node in list(padmet.dicOfNode.values()) if node.type == "reaction"])
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
    