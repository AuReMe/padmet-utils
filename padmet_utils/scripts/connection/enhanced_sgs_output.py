# -*- coding: utf-8 -*-
"""
Description:
    1./ extract all data from sgs output

    2./ If needed, convert reaction id to metacyc id

::

    usage:
        enhanced_sgs_output.py --sgs_output=FILE --padmetRef=FILE --output=FILE [--db=ID] [--mnx=FILE] [-v]
    
    options:
        -h --help     Show help.
        --sgs_output=FILE    path of a sgs run' result
        --db=ID    database origin of the reactions. will be converted to metacyc
        --padmetRef=FILE    path to padmet file corresponding to the database of reference (the repair network)
        --mnx=FILE    path to metanetx file for reactions (reac_xref.tsv)
        --output=FILE    path to tsv output file
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
    mnx_file = args["--mnx"]
    db_origin = args["--db"]
    

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
    
    #3:extract db ref from mnxref
    if mnx_file:
        if verbose: print("extract association %s to metacyc from mnxref" %db_origin)
        with open(mnx_file,'r') as f:
            data = [line.split("\t")[:2] for line in f.read().splitlines() if not line.startswith("#")]
        #k = mnxref, v: dict k' in [bigg,metacyc], v = [list]
        dict_mnxref = dict([(v[1],dict()) for v in data])
        for value in data:
            xref, mnx = value
            if (xref.startswith(db_origin) or xref.startswith("metacyc")):
                db, _id = xref.split(":")
                try:
                    dict_mnxref[mnx][db].add(_id)
                except KeyError:
                    dict_mnxref[mnx].update({db: set([_id])})
        for k,v in list(dict_mnxref.items()):
            if len(list(v.keys())) < 2:
                dict_mnxref.pop(k)
                
        
        #map reactions id
        assoc_sgs_id_reactions_mapped = {}
        for sgs_id, set_rxn in assoc_sgs_id_reactions.items():
            for assoc_dict in list(dict_mnxref.values()):
                for _id in assoc_dict[db_origin]:
                    if _id in set_rxn:
                        try:
                            assoc_sgs_id_reactions_mapped[sgs_id].update(assoc_dict["metacyc"])
                        except KeyError:
                            assoc_sgs_id_reactions_mapped[sgs_id] = assoc_dict["metacyc"]
     
        """
        for each set of reactions associated to a sgs, get the intersection of reactions associated to each pathways.
        if an intersection length is > 0 ==> the sgs covere a part of this pwy.
        Extract which reactions and the ratio on all the reactions in the pathways.
        """
        for sgs_id, rxns_in_sgs in assoc_sgs_id_reactions_mapped.items():
            for pwy_id, rxns_in_pwy in assoc_pathways_reactions.items():
                rxns_inter = rxns_in_sgs.intersection(rxns_in_pwy)
                len_rxns_inter = len(rxns_inter)
                if len_rxns_inter > 0 and float(len_rxns_inter)/float(len(rxns_in_pwy)) > float(1)/float(3):
                    data = str(len_rxns_inter)+"/"+str(len(rxns_in_pwy))+": "
                    data += ";".join([rxn_id for rxn_id in rxns_inter])
                    try:
                        assoc_pathways_sgs[pwy_id][sgs_id] = data
                    except KeyError:
                        assoc_pathways_sgs[pwy_id] = {sgs_id:data}
    
   

    if verbose: print("creating output %s" %output)
    with open(output, 'w') as f:
        header = "\t".join(["SGS ID","GENES ID","RXN ID_"+db_origin ,"RXN ID_metacyc"])+"\n"
        f.write(header)
        for sgs_id in all_sgs: 
            genes_id = assoc_sgs_id_genes[sgs_id]
            line = [sgs_id,";".join(genes_id)]
            rxns_id_bigg = ";".join(assoc_sgs_id_reactions[sgs_id])
            line.append(rxns_id_bigg)
            if mnx_file:
                try:
                    rxns_id_metacyc = ";".join(assoc_sgs_id_reactions_mapped[sgs_id])
                except KeyError:
                    rxns_id_metacyc = "NA"
            else:
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
