# -*- coding: utf-8 -*-
"""
Description:
    Extract all data from sgs output


::

    usage:
        enhanced_sgs_output.py --sgs_output=FILE --padmetRef=FILE --output=FILE [-v]
    
    options:
        -h --help     Show help.
        --sgs_output=FILE    path of a sgs run result
        --padmetRef=FILE    path to padmet file corresponding to the database of reference (the repair network)
        --output=FILE    path to tsv output file
        -v   print info
"""
from padmet.classes import PadmetRef
import csv
import re
import docopt

def main():
    args = docopt.docopt(__doc__)
    
    sgs_output = args["--sgs_output"]
    output = args["--output"]
    verbose = args["-v"]
    padmetRef = PadmetRef(args["--padmetRef"])
    output = args["--output"]
    enhanced_sgs_output(sgs_output, padmetRef, output, verbose)

def enhanced_sgs_output(sgs_output, padmetRef, output, verbose=False):
    """
    Extract all data from sgs output

    Parameters
    ----------
    sgs_output: str
        path of a sgs run result
    padmetRef: padmet.classes.PadmetRef
        padmet containing the database of reference, need to calculat pathway completion rate
    output: str
        pathname of the output folder
    verbose: bool
        if True print information
    """
    #k = sgs_id, v = set of genes
    assoc_sgs_id_genes = {}
    #k = sgs_id, v = set of reactions
    assoc_sgs_id_reactions = {}
    #k = pwy_id, v = set of reactions
    assoc_pathways_reactions = {}
    #k = pwy_id, v = {sgs_id:["x/X; rxn_ids],...}
    assoc_pathways_sgs = {}
    all_sgs = []
    #1: extract data from shogen output.
    count = 1
    if verbose: print("Reading sgs output")
    with open(sgs_output, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            genes, reactions = row["gene_set"], row["reaction_set"]
            genes = re.sub("ac:|nc:|na:|\"","",genes).split(" ")
            reactions = re.sub("nc:|\"|R_","",reactions).split(" ")
            sgs_id = "sgs-"+str(count)+"_"+str(len(genes))
            assoc_sgs_id_genes[sgs_id] = set(genes)
            assoc_sgs_id_reactions[sgs_id] = set(reactions)
            count += 1
            all_sgs.append(sgs_id)
    
    #2: extract all pathways from database Ref
    if verbose: print("extract all pathways from database Ref")
    all_pathways_id = (node.id for node in list(padmetRef.dicOfNode.values()) if node.type == "pathway")
    for pwy_id in all_pathways_id:
        rxns_assoc = [rlt.id_in for rlt in padmetRef.dicOfRelationOut[pwy_id] 
        if rlt.type == "is_in_pathway" and padmetRef.dicOfNode[rlt.id_in].type == "reaction"]
        assoc_pathways_reactions[pwy_id] = set(rxns_assoc)
    
     
    if verbose: print("creating output %s" %output)
    with open(output, 'w') as f:
        header = "\t".join(["SGS ID","GENES ID","RXN ID"])+"\n"
        f.write(header)
        for sgs_id in all_sgs: 
            genes_id = assoc_sgs_id_genes[sgs_id]
            line = [sgs_id,";".join(genes_id)]
            rxns_id_bigg = ";".join(assoc_sgs_id_reactions[sgs_id])
            line.append(rxns_id_bigg)
            rxns_id_metacyc = ";".join(assoc_sgs_id_reactions[sgs_id])
            line.append(rxns_id_metacyc)
            line = "\t".join(line)+"\n"
            f.write(line)
        f.write("\n\n\n")
    
        header = ["PATHWAYS ID\SGS ID"]+all_sgs
        f.write("\t".join(header)+"\n")
        for pwy_id,sgs_dict in assoc_pathways_sgs.items():
            line = [""]*len(header)
            try:
                line[0] = pwy_id+"/"+padmetRef.dicOfNode[pwy_id].misc["COMMON_NAME"][0]
            except KeyError:
                line[0] = pwy_id
            for sgs_id, data in sgs_dict.items():
                sgs_index = header.index(sgs_id)
                line[sgs_index] = data
            line = "\t".join(line)+"\n"
            f.write(line)
                        
if __name__ == "__main__":
    main()                  
