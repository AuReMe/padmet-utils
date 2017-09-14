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
convert a padmet representing a database (padmetRef) and/or a padmet representing a model (padmetSpec)
to tsv files for askomics.

1./ Folder creation
given the output directory. Create this directory if required and create a folder
padmetRef filename and/or padmetSpec filename

2./
2.1/ For padmetRef:
2.1.a/ Nodes
get all reactions nodes => extract data from misc with extract_nodes(rxn_nodes, "reaction", "../rxn.tsv")
get all compounds nodes => extract data from misc with extract_nodes(cpd_nodes, "compounds", "../cpd.tsv")
get all pathways nodes => extract data from misc with extract_nodes(pwy_nodes, "pathway", "../pwy.tsv")
get all xrefs nodes => extract data from misc with extract_nodes(xref_nodes, "xref", "../xref.tsv")
2.1.b/ Relations
for each rxn in rxn nodes:
    get all rlt consumes/produces => create list of data with extract_rxn_cpd(rxn_cpd_rlt)
        fieldnames = "rxn_cpd","concerns@reaction","consumes@compound","produces@compound","stoichiometry","compartment"
    get all rlt is_in_pathway => create list of data with extract_rxn_pwy(rxn_pwy_rlt)
        fieldnames = "rxn_pwy","concerns@reaction","in_pwy@pathway"
    get all rlt has_xref => create list of data with extract_entity_xref(rxn_xref_rlt)
for each cpd in cpd nodes:
    get all rlt has_xref => update previous list of data with extract_entity_xref(cpd_xref_rlt)
        fieldnames = "entity_xref","concerns@reaction","concerns@compound","has_xref@xref"

usage:
    padmet_to_tsv.py --padmetSpec=FILE [--padmetRef=FILE] --output_dir=DIR [--url=FILE] [-v]
    padmet_to_tsv.py --padmetRef=FILE [--padmetSpec=FILE] --output_dir=DIR [--url=FILE] [-v]

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathname of the padmet representing the network to convert
    --padmetRef=FILE    pathname of the padmet representing the database
    --url=FILE    file with associated url to each element. line = url,entity id. Sep = "\t"
    --output_dir=DIR
    -v
"""
from padmet.padmetSpec import PadmetSpec
from padmet.padmetRef import PadmetRef
from itertools import izip_longest
import os
import csv
import docopt

def main():
    global url_dict
    args = docopt.docopt(__doc__)
    padmetSpec_file = args["--padmetSpec"]
    padmetRef_file = args["--padmetRef"]
    output_dir = args["--output_dir"]
    url_file = args["--url"]
    verbose = args["-v"]
    
    #check if output_dir exist, else create it
    if not output_dir.endswith("/"):
        output_dir += "/"
    if not os.path.isdir(output_dir):
        if verbose: print("Creating folder %s" %output_dir)
        os.makedirs(output_dir)
    #loading padmetSpec
    if padmetSpec_file:
        if verbose: print("Loading %s" %padmetSpec_file)
        padmetSpec = PadmetSpec(padmetSpec_file)
        padmetSpec_name = os.path.splitext(os.path.basename(padmetSpec_file))[0]
        padmetSpec_folder = output_dir+padmetSpec_name+"/"
        if not os.path.isdir(padmetSpec_folder):
            if verbose: print("Creating folder %s" %padmetSpec_folder)
            os.makedirs(padmetSpec_folder)

    #if padmetRef given, create folder for padmetRef
    if padmetRef_file:
        if verbose: print("Loading %s" %padmetRef_file)
        padmetRef = PadmetRef(padmetRef_file)
        padmetRef_name = os.path.splitext(os.path.basename(padmetRef_file))[0]
        padmetRef_folder = output_dir+padmetRef_name+"/"
        if not os.path.isdir(padmetRef_folder):
            if verbose: print("Creating folder %s" %padmetRef_folder)
            os.makedirs(padmetRef_folder)

    #if url file. Create dict: k = entity id, v = url. Default = ???
    url_dict = {"default" : "http://www.semanticweb.org/irisa/ontologies/2016/1/igepp-ontology#"}
    if url_file:
        with open(url_file, 'r') as f:
            url_dict.update(dict([line.split("\t") for line in f.read().splitlines()]))
    url_default = url_dict["default"]
    #NODES
    #Converting nodes data to tsv format
    if padmetRef_file:
        if verbose: print("Extracting nodes from %s" %padmetRef_name)
        with open(padmetRef_folder+"metabolic_network.tsv", 'w') as f:
            fieldnames = ["metabolic_network","name"]
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(fieldnames)
            try:
                url_padmetRef_name = url_dict[padmetRef_name]
            except KeyError:
                url_padmetRef_name = url_dict["default"]
            writer.writerow([url_padmetRef_name+padmetRef_name, padmetRef_name])

        if verbose: print("\tExtracting reactions")
        all_rxn_nodes = [node for node in padmetRef.dicOfNode.values() if node.type == "reaction"]
        if all_rxn_nodes: extract_nodes(all_rxn_nodes, "reaction", padmetRef_folder+"rxn.tsv", {"in@metabolic_network":padmetRef_name})
        if verbose: print("\t%s reactions" %len(all_rxn_nodes))

        if verbose: print("\tExtracting compounds")
        all_cpd_nodes = set([padmetRef.dicOfNode[rlt.id_out] for rlt in padmetRef.getAllRelation() if rlt.type in ["consumes","produces"]])
        if all_cpd_nodes: extract_nodes(all_cpd_nodes, "compound", padmetRef_folder+"cpd.tsv")
        if verbose: print("\t%s compounds" %len(all_cpd_nodes))
    
        if verbose: print("\tExtracting pathways")
        all_pwy_nodes = [node for node in padmetRef.dicOfNode.values() if node.type == "pathway"]
        if all_pwy_nodes: extract_nodes(all_pwy_nodes, "pathway", padmetRef_folder+"pwy.tsv")
        if verbose: print("\t%s pathways" %len(all_pwy_nodes))

        if verbose: print("\tExtracting xrefs")
        all_xref_nodes = [node for node in padmetRef.dicOfNode.values() if node.type == "xref"]
        if all_xref_nodes: extract_nodes(all_xref_nodes, "xref", padmetRef_folder+"xref.tsv")
        if verbose: print("\t%s xrefs" %len(all_xref_nodes))

    else:
        if verbose: print("No given padmetRef")
        all_rxn_nodes, all_cpd_nodes, all_pwy_nodes, all_xref_nodes=[[]]*4

    if padmetSpec_file:
        with open(padmetSpec_folder+"metabolic_network.tsv", 'w') as f:
            fieldnames = ["metabolic_network","name"]
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(fieldnames)
            try:
                url_padmetSpec_name = url_dict[padmetSpec_name]
            except KeyError:
                url_padmetSpec_name = url_dict["default"]
            writer.writerow([url_padmetSpec_name+padmetSpec_name,padmetSpec_name])

        if verbose: print("Extracting nodes from %s" %padmetSpec_name)

        if verbose: print("\tExtracting all reactions and sources")
        all_sources = set()
        spec_all_rxn_nodes = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"] 
        [all_sources.update(set(node.misc.get("SOURCE",[]))) for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
        [all_sources.update(set(rlt.misc.get("ASSIGNMENT",[]))) for rlt in padmetSpec.getAllRelation() if rlt.type == "is_linked_to"]
        if spec_all_rxn_nodes:
            fieldnames = ["reaction","in@metabolic_network"]
            with open(padmetSpec_folder+"rxn.tsv", 'w') as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(fieldnames)
                try:
                    url_rxn = url_dict["reaction"]
                except KeyError:
                    url_rxn = url_dict["default"]                
                for rxn_node in spec_all_rxn_nodes:
                    writer.writerow([url_rxn+rxn_node.id, url_padmetSpec_name+padmetSpec_name])
        if all_sources:
            fieldnames = ["reconstruction_information","name"]
            with open(padmetSpec_folder+"sources.tsv", 'w') as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(fieldnames)
                try:
                    url_src = url_dict["reconstruction_information"]
                except KeyError:
                    url_src = url_dict["default"]                
                for src in all_sources:
                    writer.writerow([url_src+src, src])
            
        if verbose: print("\t%s reactions" %len(spec_all_rxn_nodes))

        if verbose: print("\tExtracting reactions not in padmetRef")
        spec_unique_rxn_nodes = [padmetSpec.dicOfNode[node_id] for node_id in set([node.id for node in spec_all_rxn_nodes]).difference(set([node.id for node in all_rxn_nodes]))] 
        if spec_unique_rxn_nodes: extract_nodes(spec_unique_rxn_nodes, "reaction", padmetSpec_folder+"unique_rxn.tsv", {"in@metabolic_network":padmetSpec_name})
        if verbose: print("\t%s reactions" %len(spec_unique_rxn_nodes))

        if verbose: print("\tExtractingcompounds not in padmetRef")
        spec_unique_cpd_nodes = [padmetSpec.dicOfNode[node_id] for node_id in set([rlt.id_out for rlt in padmetSpec.getAllRelation() if rlt.type in ["consumes","produces"]]).difference(set([node.id for node in all_cpd_nodes]))]
        if spec_unique_cpd_nodes: extract_nodes(spec_unique_cpd_nodes, "compound", padmetSpec_folder+"cpd.tsv")
        if verbose: print("\t%s compounds" %len(spec_unique_cpd_nodes))

        if verbose: print("\tExtracting pathways not in padmetRef")
        spec_unique_pwy_nodes = [padmetSpec.dicOfNode[node_id] for node_id in set([node.id for node in padmetSpec.dicOfNode.values() if node.type == "pathway"]).difference(set([node.id for node in all_pwy_nodes]))]
        if spec_unique_pwy_nodes: extract_nodes(spec_unique_pwy_nodes, "pathway", padmetSpec_folder+"pwy.tsv")
        if verbose: print("\t%s pathways" %len(spec_unique_pwy_nodes))
        
        if verbose: print("\tExtracting all genes")
        spec_genes_nodes = [node for node in padmetSpec.dicOfNode.values() if node.type == "gene"]
        if spec_genes_nodes: extract_nodes(spec_genes_nodes, "gene", padmetSpec_folder+"gene.tsv", opt_col = {"in@metabolic_network":padmetSpec_name})
        if verbose: print("\t%s genes" %len(spec_genes_nodes))




    #RELATIONS
    #Converting relations data to tsv format
    if padmetRef_file:
        if verbose: print("Extracting relations from %s" %padmetRef_name)
        rxn_cpd_data = []
        rxn_pwy_data = []
        entity_xref_data = []
        if verbose: print("\tExtracting relations reaction-[consumes/produces]-compound")
        if verbose: print("\tExtracting relations reaction-is_inclued_in-pathway")
        if verbose: print("\tExtracting relations reactions-has_xref-xref")

        for rxn_node in all_rxn_nodes:
            rxn_id = rxn_node.id
            #if verbose: print("Reaction %s" %rxn_id)

            #all consumes/produces relations
            cp_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]]
            rxn_cpd_data += extract_rxn_cpd(cp_rlt)
            #all is_in_pathway relations
            pwy_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway"]
            if pwy_rlt: rxn_pwy_data += extract_rxn_pwy(pwy_rlt)
            #all has_xref relations
            rxn_xref_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[rxn_id] if rlt.type == "has_xref"]
            if rxn_xref_rlt: entity_xref_data += extract_entity_xref(rxn_xref_rlt, padmetRef)

        if verbose: print("\tExtracting relations compound-has_xref-xref")
        for cpd_node in all_cpd_nodes:
            cpd_id = cpd_node.id
            try:
                cpd_xref_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[cpd_id] if rlt.type == "has_xref"]
                if cpd_xref_rlt: entity_xref_data += extract_entity_xref(cpd_xref_rlt, padmetRef)
            except KeyError:
                pass

        if verbose: print("\tExtracting relations pwy-has_xref-xref")
        for pwy_node in all_pwy_nodes:
            pwy_id = pwy_node.id
            try:
                pwy_xref_rlt = [rlt for rlt in padmetRef.dicOfRelationIn[pwy_id] if rlt.type == "has_xref"]
                if pwy_xref_rlt: entity_xref_data += extract_entity_xref(pwy_xref_rlt, padmetRef)
            except KeyError:
                pass
        
        if rxn_cpd_data:
            if verbose: print("\t\tCreating rxn_cpd.tsv")
            rxn_cpd_file(rxn_cpd_data, padmetRef_folder+"rxn_cpd.tsv")
        if rxn_pwy_data:
            if verbose: print("\t\tCreating rxn_pwy.tsv")
            rxn_pwy_file(rxn_pwy_data, padmetRef_folder+"rxn_pwy.tsv")
        if entity_xref_data:
            if verbose: print("\t\tCreating entity_xref.tsv")
            entity_xref_file(entity_xref_data, padmetRef_folder+"entity_xref.tsv")

    if padmetSpec_file:
        if verbose: print("Extracting relations from %s" %padmetSpec_name)
        rxn_cpd_data = []
        rxn_pwy_data = []
        rxn_gene_data = []
        entity_xref_data = []

        fieldnames = ["rxn_reconstruction_info","concers@reaction","has_metadata@reconstruction_information","concerns@metabolic_network"]
        with open(padmetSpec_folder+"rxn_sources.tsv", 'w') as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(fieldnames)
            for rxn_node in spec_all_rxn_nodes:
                for src in rxn_node.misc.get("SOURCE",[]):
                    line = [url_rxn+rxn_node.id, url_src+src, url_padmetSpec_name+padmetSpec_name]
                    line.insert(0,"_".join(line))
                    writer.writerow(line)            

        if padmetRef:
            if verbose: print("\tPadmetRef given, extracting relations for unique reactions only")
            if verbose: print("\tExtracting relations reaction-[consumes/produces]-compound")
            if verbose: print("\tExtracting relations reaction-is_in_pathway-pathway")
            if verbose: print("\tExtracting relations reactions-has_xref-xref")

            for rxn_node in spec_unique_rxn_nodes:
                rxn_id = rxn_node.id
                #if verbose: print("Reaction %s" %rxn_id)

                #all consumes/produces relations
                rxn_cpd_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]]
                rxn_cpd_data += extract_rxn_cpd(rxn_cpd_rlt)
                #all is_in_pathway relations
                rxn_pwy_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_in_pathway"]
                rxn_pwy_data += extract_rxn_pwy(rxn_pwy_rlt)
                #all has_xref relations
                rxn_xref_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "has_xref"]
                entity_xref_data += extract_entity_xref(rxn_xref_rlt, padmetSpec)
                #all sources in relations
            if verbose: print("\tExtracting relations compound-has_xref-xref")
            for cpd_node in spec_unique_cpd_nodes:
                cpd_id = cpd_node.id
                #if verbose: print("Compound %s" %cpd_id)
                
                try:
                    cpd_xref_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[cpd_id] if rlt.type == "has_xref"]
                except KeyError:
                    pass
                entity_xref_data += extract_entity_xref(cpd_xref_rlt, padmetSpec)

            if verbose: print("\tExtracting relations reactions-is_linked_to-gene")
            for rxn_node in spec_all_rxn_nodes:
                rxn_id = rxn_node.id
                #if verbose: print("reaction %s" %rxn_id)

                #all is_linked_to relations
                rxn_gene_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]
                rxn_gene_data += extract_rxn_gene(rxn_gene_rlt)

            if verbose: print("\tExtracting pathways's completion rate and creating pwy_rate.tsv")
            pwy_rate(padmetRef, padmetSpec, padmetSpec_name, padmetSpec_folder+"pwy_rate.tsv")

        else:
            if verbose: print("No padmetRef given, extracting relations for all reactions")
                
            for rxn_node in spec_all_rxn_nodes:
                rxn_id = rxn_node.id
                #all is_linked_to relations
                rxn_gene_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]
                rxn_gene_data += extract_rxn_gene(rxn_gene_rlt)
         

        if rxn_cpd_data:
            if verbose: print("\t\tCreating rxn_cpd.tsv")
            rxn_cpd_file(rxn_cpd_data, padmetSpec_folder+"rxn_cpd.tsv")
        if rxn_pwy_data:
            if verbose: print("\t\tCreating rxn_pwy.tsv")
            rxn_pwy_file(rxn_pwy_data, padmetSpec_folder+"rxn_pwy.tsv")
        if entity_xref_data:
            if verbose: print("\t\tCreating entity_xref.tsv")
            entity_xref_file(entity_xref_data, padmetSpec_folder+"entity_xref.tsv")
        if rxn_gene_data:
            if verbose: print("\t\tCreating rxn_gene.tsv")
            rxn_gene_file(rxn_gene_data, padmetSpec_folder+"rxn_gene.tsv")


def extract_nodes(nodes, entity_id, output, opt_col = {}):
    """
    for n nodes in nodes. for each node.misc = {A:['x'],B:['y','z']}
    create a file with line = [node.id,A[0],B[0]],[node.id,"",B[1]]
    the order is defined in fieldnames.
    """
    all_keys = set()
    [all_keys.update(node.misc.keys()) for node in nodes]
    try:
        url_entity = url_dict[entity_id]
    except KeyError:
        url_entity = url_dict["default"]
    if opt_col:
        fieldnames = [url_entity+entity_id] + opt_col.keys() + list(all_keys)
    else:
        fieldnames = [entity_id] + list(all_keys)
    try:
        fieldnames.remove("SOURCE")
    except ValueError:
        pass
    
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        for node in nodes:
            data = [node.misc.get(col,[""]) for col in fieldnames]
            data = izip_longest(*data,fillvalue="")
            for d in data:
                d = list(d)
                d[0] = node.id
                if opt_col:
                    for k,v in opt_col.items():
                        k_index = fieldnames.index(k)
                        d[k_index] = v
                writer.writerow(d)
            
def extract_pwy(padmet):
    """
    from padmet return a dict, k = pwy_id, v = set of rxn_id in pwy
    """
    #get all pathways in ref in dict: k=pwy_id, v = set(rxn_id in pwy)
    #get all pathways in in spec dict: k=pwy_id, v = set(rxn_id in pwy)
    pathways_rxn_dict = dict([(node.id, set([rlt.id_in for rlt in padmet.dicOfRelationOut[node.id] 
    if rlt.type == "is_in_pathway" and padmet.dicOfNode[rlt.id_in].type == "reaction"])) for node in padmet.dicOfNode.values() if node.type == "pathway"])
    #pop k,v if len(v) == 0 (if not v)
    [pathways_rxn_dict.pop(k) for k,v in pathways_rxn_dict.items() if not v]
    return pathways_rxn_dict

def extract_rxn_cpd(rxn_cpd_rlt):
    """
    for rlt in rxn_cpd_rlt, append in data: [rxn_id,cpd_id(consumed),'',stoich,compartment]
    and/or [rxn_id,'',cpd_id(produced),stoich,compartment]. The value in index 0
    is a merge of all data to create a unique relation id
    """
    data = []
    try:
        url_rxn = url_dict["reaction"]
    except KeyError:
        url_rxn = url_dict["default"]
    try:
        url_cpd = url_dict["compound"]
    except KeyError:
        url_cpd = url_dict["default"]
    
    for rlt in rxn_cpd_rlt:
        rxn_id = url_rxn + rlt.id_in
        cpd_id = url_cpd + rlt.id_out
        stoich = rlt.misc["STOICHIOMETRY"][0]
        if "n" in stoich: stoich = "1"
        if rlt.type == "consumes":
            line = [rxn_id, cpd_id, "", stoich, rlt.misc["COMPARTMENT"][0]]
            line.insert(0, "_".join(line))
        elif rlt.type == "produces":
            line = [rxn_id, "", cpd_id, stoich, rlt.misc["COMPARTMENT"][0]]
            line.insert(0, "_".join(line))
        data.append(line)
    return data

def rxn_cpd_file(data, output):
    """
    from data obtained with extract_rxn_cpd(), create file rxn_cpd
    """
    fieldnames = ["rxn_cpd","concerns@reaction","consumes@compound","produces@compound","stoichiometry","compartment"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]

def extract_rxn_pwy(rxn_pwy_rlt):
    """
    for rlt in rxn_pwy_rlt, append in data: [rxn_id,pwy_id]. The value in index 0
    is a merge of all data to create a unique relation id
    """
    data = []
    try:
        url_rxn = url_dict["reaction"]
    except KeyError:
        url_rxn = url_dict["default"]
    try:
        url_pwy = url_dict["pathway"]
    except KeyError:
        url_pwy = url_dict["default"]
    for rlt in rxn_pwy_rlt:
        rxn_id = url_rxn + rlt.id_in
        pwy_id = url_pwy + rlt.id_out
        line = [rxn_id, pwy_id]
        line.insert(0,"_".join(line))
        data.append(line)
    return data

def rxn_pwy_file(data, output):
    fieldnames = ["rxn_pwy","concerns@reaction","in_included_in@pathway"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]


def pwy_rate(padmetRef, padmetSpec, metabolic_network, output):
    """
    pwy rate in padmetSpec is calculated based on padmetRef
    """
    all_pathways_dict = extract_pwy(padmetRef)    
    network_pathways_dict = extract_pwy(padmetSpec)
    try:
        url_pwy = url_dict["pathway"]
    except KeyError:
        url_pwy = url_dict["default"]
    try:
        url_padmetSpec = url_dict[padmetSpec_name]
    except KeyError:
        url_padmetSpec = url_dict["default"]

    if network_pathways_dict:
        fieldnames = ["pathway","is_part_of@metabolic_network","RATE"]
        with open(output, 'w') as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(fieldnames)
            for pwy_id, rxns_set in network_pathways_dict.items():
                #pwy_id = "PWY-5497"
                #rxns_set = network_pathways_dict["PWY-5497"]
                all_rxns_set = all_pathways_dict[pwy_id]
                nb_all_rxns = len(all_rxns_set)
                nb_in_network = len(rxns_set.intersection(all_rxns_set))
                #rate nb rxn in network / nb total rxn in pwy
                rate = round(float(nb_in_network)/float(nb_all_rxns),2)
                rate = str(rate).replace(",",".")
                writer.writerow([url_pwy + pwy_id, url_padmetSpec + metabolic_network, rate])

def extract_entity_xref(entity_xref_rlt, padmet):
    data = []
    for rlt in entity_xref_rlt:
        entity_id = rlt.id_in
        xref_id = rlt.id_out
        if padmet.dicOfNode[entity_id].type == "reaction":
            try:
                url_rxn = url_dict["reaction"]
            except KeyError:
                url_rxn = url_dict["default"] ~#TODO IICCCCCCCCIIIII !!!!!!!!
            line = [xref_id,entity_id, "", ""]
            line.insert(0, "_".join(line))
        elif padmet.dicOfNode[entity_id].type == "pathway":
            line = [xref_id, "", entity_id, ""]
            line.insert(0, "_".join(line))
        else:
            line = [xref_id, "", "", entity_id]
            line.insert(0, "_".join(line))
        data.append(line)
    return data

def entity_xref_file(data, output):
    fieldnames = ["entity_xref", "has_metadata@xref", "concerns@reaction","concerns@pathway","concerns@compound"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]

def extract_rxn_gene(rxn_gene_rlt):
    data = []
    for rlt in rxn_gene_rlt:
        rxn_id = rlt.id_in
        gene_id = rlt.id_out
        for src in rlt.misc["ASSIGNMENT"]:
            line = [rxn_id, gene_id, src]
            line.insert(0, "_".join(line))
            data.append(line)
    return data

def rxn_gene_file(data, output):
    fieldnames = ["rxn_gene","concerns@reaction","is_linked_to@gene","has_metadata@reconstruction_information"]
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        [writer.writerow(d) for d in data]

if __name__ == "__main__":
    main()
