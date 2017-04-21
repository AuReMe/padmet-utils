# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr

Description:
1./From padmetRef and padmetSpec get all pathways in 2 dict: all_pathways_dict & network_pathways_dict
K = pathway id, v = set of reactions id (and only reactions, no pathways)
Function: extract_pathways(padmet), return dict
2./ From padmetSpec get the name of all the origin files: all_origin
Function: get_all_origin(padmetSpec), return list
2./ In some case, we can merge some origin file in categories
    given a dictionnary file with line = origin_file/category id 
        arg: categories = ""
3./ For each pwy in the padmetSpec:
3.a/ Count nb all reactions in pwy, nb reactions in padmetSpec, rate
3.b/ For each rxn in pwy:
get the orignin(s) in 2 dict
    category_rxn_dict, K = category of the origin, V = set of reactions id
    rxn_category_dict, K = rxn_id, V = set of origin
count the nb of reactions with unique origin (rxn with len(v) == 1)
count the nb of categories involved (len(category_rxn_dict.keys()))
#Particular case: in tiso we dont want to take in accout some category (like pwtool and primary network)
#Particular case: use EC number specific in each category to filtre real specific reactions
return file with:
header = ["PWY_ID","NbTotalRXN","NbRXNFound","Rate"] + list(all_sources (full-specific)) + ["RXN from unique source","NbRXN from unique source","Involved sources","Nb sources"]

usage:
    extract_pwy_sources.py --padmetSpec=FILE --padmetRef=FILE --categories=FILE --output=FILE

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathname of the padmet representing the network to analyse
    --padmetRef=FILE    pathname of the padmet representing the database
    --categories=FILE    pathname of the file in which the origins are categorized, line = 'origin' 'category' (sep = '\t')
    --output=FILE    pathname of the output
    -v   print info.
"""
from lib.padmetRef import PadmetRef
from lib.padmetSpec import PadmetSpec
from time import time
import docopt
import os

def main():
    args = docopt.docopt(__doc__)
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    padmetRef = PadmetRef(args["--padmetRef"])
    origins_caterogies_file = args["--categories"]
    output = args["--output"]
    chronoDepart = time()
    all_pathways_dict = extract_pwy(padmetRef)
    network_pathways_dict = extract_pwy(padmetSpec)
    origins_caterogies_dict = create_origins_caterogies_dict(origins_caterogies_file)
    
    with open(output, 'w') as f:
        header = ["PWY_ID","NbTotalRXN","NbRXNFound","Rate"]
        for c in set(origins_caterogies_dict.values()):
            header += [c+"-Full",c+"-Specific"]
        header += ["RXN from unique source","NbRXN from unique source"]
        header += ["Involved sources","Nb involved sources"]
        header += ["Unique involved sources","Nb unique involved sources"]
        f.write("\t".join(header)+"\n")
        for pwy_id, rxns_set in network_pathways_dict.items():
            #pwy_id = "PWY-5497"
            #rxns_set = network_pathways_dict["PWY-5497"]
            line = [""] * len(header)
            all_rxns_set = all_pathways_dict[pwy_id]
            nb_all_rxns = len(all_rxns_set)
            nb_in_network = len(rxns_set.intersection(all_rxns_set))
            #rate nb rxn in network / nb total rxn in pwy
            rate = round(float(nb_in_network)/float(nb_all_rxns),2)
            rate = str(rate).replace(".",",")
            line[:4] = [pwy_id,str(nb_all_rxns),str(nb_in_network),rate]
    
            #dict K= organism, v = list of rxn
            category_rxn_dict = {}
            rxn_category_dict = {}
            for rxn_id in rxns_set:
                #get all sources of a rxn
                rxn_origins = [padmetSpec.dicOfNode[rlt.id_out].misc["ORIGIN_FILE"][0] for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "has_suppData"]
                rxn_category_dict[rxn_id] = [origins_caterogies_dict[origin_id] for origin_id in rxn_origins]
                for origin_id in rxn_origins:
                    try:
                        category_rxn_dict[origins_caterogies_dict[origin_id]].add(rxn_id)
                    except KeyError:
                        category_rxn_dict[origins_caterogies_dict[origin_id]] = set([rxn_id])
            #add in line in the correct index of each source, all rxn sep = ";"
            rxn_from_unique_category = set([k for k,v in rxn_category_dict.items() if len(v) == 1 and v[0] not in ["PWT","PRIMARY"]])
            involved_category = category_rxn_dict.keys()                        
            
            for category, rxn_set in category_rxn_dict.items():
                category_full_index = header.index(category+"-Full")
                category_spe_index = header.index(category+"-Specific")
                rxn_c_specific = rxn_from_unique_category.intersection(rxn_set)
                line[category_full_index:category_spe_index+1] = [";".join(rxn_set),";".join(rxn_c_specific)] 
            unique_involved_category = set([rxn_category_dict[rxn_id][0] for rxn_id in rxn_from_unique_category])
            nb_unique_involved_category = len(unique_involved_category)
            line[header.index("RXN from unique source")] = ";".join(rxn_from_unique_category)
            line[header.index("NbRXN from unique source")] = str(len(rxn_from_unique_category))
            line[header.index("Involved sources")] = ";".join(involved_category)
            line[header.index("Nb involved sources")] = str(len(involved_category))
            line[header.index("Unique involved sources")] = ";".join(unique_involved_category)
            line[header.index("Nb unique involved sources")] = str(nb_unique_involved_category)
            line = "\t".join(line)+"\n"
            f.write(line)

    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    print "done in: ", chrono, "s !"

def extract_pwy(padmet):
    #get all pathways in ref in dict: k=pwy_id, v = set(rxn_id in pwy)
    #get all pathways in in spec dict: k=pwy_id, v = set(rxn_id in pwy)
    pathways_rxn_dict = dict([(node.id, set([rlt.id_in for rlt in padmet.dicOfRelationOut[node.id] 
    if rlt.type == "is_in_pathway" and padmet.dicOfNode[rlt.id_in].type == "reaction"])) for node in padmet.dicOfNode.values() if node.type == "pathway"])
    #pop k,v if len(v) == 0 (if not v)
    [pathways_rxn_dict.pop(k) for k,v in pathways_rxn_dict.items() if not v]
    return pathways_rxn_dict

def create_origins_caterogies_dict(origins_caterogies_file):
    with open(origins_caterogies_file, 'r') as f:
        origins_caterogies_dict = dict([line.split("\t") for line in f.read().splitlines()])
    return origins_caterogies_dict

       

"""
with open(origins_caterogies_file, 'w') as f:
    all_origins = get_all_origins(padmetSpec)
    f.write("\n".join(all_sources))

def get_all_origins(padmetSpec):
    all_origins = list(set([node.misc["ORIGIN_FILE"][0] for node in padmetSpec.dicOfNode.values() if node.type == "suppData"]))
    return all_origins
"""
if __name__ == "__main__":
    main()