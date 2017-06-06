# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:

usage:
    get_pwy_from_rxn.py --padmet=FILE --cdt_file=FILE  --output=FILE

options:
    -h --help     Show help.
    --reaction_file=FILE    pathname of the file containing the reactions id, 1/line 
    --padmetRef=FILE    pathname of the padmet representing the database.
    --output=FILE    pathname of the file with line = pathway id, all reactions id, reactions ids from reaction file, ratio. sep = "\t"
"""
from padmet.padmetSpec import PadmetSpec
import docopt
import csv


def main():
    """    
    args = docopt.docopt(__doc__)
    cdt_file = args["--condition_file"]
    padmet_file = args["--padmet"]
    output = args["--output"]
    """
    padmet_file = "/home/maite/padmet_tools/test/arnaud/draft.padmet"
    cdt_file = "/home/maite/padmet_tools/test/arnaud/diff_TenorD3_vs_TenorD0_T2_SANS_T2MD6TB.csv"
    output = "/home/maite/padmet_tools/test/arnaud/output.csv"

    padmet = PadmetSpec(padmet_file)
    
    #read cdt_file:
    with open(cdt_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        #k: gene_id, v = float, FC value
        dict_gene_exp = dict([(dict_line[""],float((dict_line["log2FoldChange"].replace(",",".")))) for dict_line in reader])

    #k: pwy_id, v = set of gene_id
    all_pwys = [node.id for node in padmet.dicOfNode.values() if node.type == "pathway"]

    with open(output, 'w') as f:
        header = ["pathway_id", "common_name","dict_reactions_genes","dict_reactions_FC","UP/DOWN/BOTH"]
        line = ";".join(header)+"\n"
        f.write(line)
        for pwy_id in all_pwys:
            try:
                pwy_cname = padmet.dicOfNode[pwy_id].misc["COMMON_NAME"][0]
            except KeyError:
                pwy_cname = ""
            #pwy_id = all_pwys[0]
            all_gene_associated = set()
            rxn_associated = [rlt.id_in for rlt in padmet.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway"]
            dict_rxn_gene = {}
            for rxn_id in rxn_associated:
                #rxn_id = rxn_associated[0]
                gene_associated = set([rlt.id_out for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"])
                all_gene_associated.update(gene_associated)
                dict_rxn_gene[rxn_id] = list(gene_associated)
            if len(all_gene_associated) > 0:
                dict_specific_gene_exp = dict([(k,v) for k,v in dict_gene_exp.items() if k in all_gene_associated])
                if dict_specific_gene_exp:
                    nbUp = len([fc_value for fc_value in dict_specific_gene_exp.values() if fc_value > 0])
                    if nbUp == len(dict_specific_gene_exp.keys()):                
                        direction = "UP"
                    elif nbUp == 0:
                        direction = "DOWN"
                    else:
                        direction = "BOTH"
                    line = [pwy_id, pwy_cname, str(dict_rxn_gene),str(dict_specific_gene_exp),direction]
                    line = ";".join(line)+"\n"
                    f.write(line)
            
if __name__ == "__main__":
    main() 