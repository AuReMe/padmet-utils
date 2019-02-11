#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 15:53:06 2018

@author: maite
"""

from padmet.classes import PadmetSpec
import os
import csv

fold = "/home/maite/Documents/data/algues_last_gbk_ec/networks"
network_ids = next(os.walk(fold))[2]

dict_specie_data= {}
for file_id in network_ids:
    file_path = "{0}/{1}".format(fold, file_id)
    padmet = PadmetSpec(file_path)
    dict_src_rxn = {}
    for rxn_id in [node.id for node in padmet.dicOfNode.values() if node.type == "reaction"]:
        reconstruct_sources = [padmet.dicOfNode[rlt.id_out].misc['SOURCE'][0] for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"]
        for src in reconstruct_sources:
            if src.startswith("OUTPUT_ORTHOFINDER"):
                src = src.replace("OUTPUT_ORTHOFINDER_%s_FROM_"%(file_id.replace(".padmet","").upper()),"")
                if src in ['E_COLI_SSTR_K12_SUBSTR_MG1655', 'ATHALIANA', 'CLADOSIPHON_OKAMURANUS']:
                    src = "TEMPLATES"
            elif src == "GENOME":
                src = "ANNOTATION"
            try:
                dict_src_rxn[src].add(rxn_id)
            except KeyError:
                dict_src_rxn[src] = set([rxn_id])
    print(file_id)
    print([(k,len(v)) for k,v in dict_src_rxn.items()])
    print("---------------------")
    dict_specie_data[file_id.replace(".padmet","").upper()] = dict_src_rxn


output = "/home/maite/Forge/docker/comparison_workspace/workdir/sources.csv"
with open(output, 'w',) as csvfile:
    fieldnames = ['Algaes\Sources', "ANNOTATION", "TEMPLATES"]+[i.replace(".padmet","").upper() for i in network_ids]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for network_id, dict_network in dict_specie_data.items():
        line = {'Algaes\Sources':network_id}
        for src, rxns in dict_network.items():
            line[src] = len(rxns)
        writer.writerow(line)


output = "/home/maite/Forge/docker/comparison_workspace/workdir/sources_unique.csv"
with open(output, 'w',) as csvfile:
    fieldnames = ['Algaes\Sources', "ANNOTATION+TEMPLATES"]+[i.replace(".padmet","").upper() for i in network_ids]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for network_id, dict_network in dict_specie_data.items():
        print(network_id, dict_network)
        line = {'Algaes\Sources':network_id}
        rxns_base = dict_network["TEMPLATES"].union(dict_network["ANNOTATION"])
        line["ANNOTATION+TEMPLATES"] = len(rxns_base)
        dict_network.pop("ANNOTATION")
        dict_network.pop("TEMPLATES")
        for src, rxns in dict_network.items():
            unique_rxns = rxns.difference(rxns_base)
            for _src in [i for i in dict_network.keys() if i != src]:
                unique_rxns = unique_rxns.difference(dict_network[_src])
            line[src] = len(unique_rxns)
        writer.writerow(line)
