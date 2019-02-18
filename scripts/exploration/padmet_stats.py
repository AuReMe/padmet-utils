#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Create a padmet stats file containing the number of pathways, reactions,
    genes and compounds inside the padmet.

    The input is a padmet file or a folder containing multiple padmets.
    
    Create a tsv file named padmet_stats.tsv where the script have been
    launched.
"""

import argparse
import csv
import os
import pandas as pa
import sys

from padmet.classes import PadmetSpec

def main():

    parser = argparse.ArgumentParser(usage="python padmet_stats.py -p padmet_file_folder")
    parser.add_argument("-p", "--padmet", dest = "padmet_file_folder", help = "Padmet file or folder containing padmet.")

    parser_args = parser.parse_args(sys.argv[1:])

    padmet_file_folder = parser_args.padmet_file_folder

    if os.path.isdir(padmet_file_folder):
        padmet_type = "dir"
    elif os.path.isfile(padmet_file_folder):
        padmet_type = "file"
    else:
        raise TypeError("%s is not a dir or a file" %(padmet_file_folder))

    output_file = open('padmet_stats.tsv', 'w')
    output_writer = csv.writer(output_file, delimiter='\t')
    output_writer.writerow(['padmet_file', 'pathways', 'reactions', 'genes', 'compounds'])

    df_orthologs = pa.DataFrame()
    if padmet_type == "dir":
        padmet_names = [padmet_file.replace('.padmet', '').upper() for padmet_file in os.listdir(padmet_file_folder)]
        for padmet_file in os.listdir(padmet_file_folder):
            padmet_path = padmet_file_folder + '/' + padmet_file
            stats = padmet_stat(padmet_path)
            output_writer.writerow(stats)
            padmet_name = padmet_file.replace('.padmet', '').upper()
            df_temp = orthology_result(padmet_path, padmet_names, padmet_name)
            df_orthologs = df_orthologs.append(df_temp)

    if padmet_type == "file":
        stats = padmet_stat(padmet_file_folder)
        output_writer.writerow(stats)
        padmet_name = padmet_file_folder.replace('.padmet', '').upper()
        padmet_names = [padmet_name]
        df_temp = orthology_result(padmet_file_folder, padmet_names, padmet_name)
        df_orthologs = df_orthologs.append(df_temp)

    df_orthologs.to_csv('padmet_orthologs_stats.tsv', sep='\t')

    output_file.close()

def padmet_stat(padmet_file):
    padmetSpec = PadmetSpec(padmet_file)

    total_pwy_id = set()
    total_cpd_id = set()

    all_rxns = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node for node in padmetSpec.dicOfNode.values() if node.type == "gene"]

    for rxn_node in all_rxns:
        total_cpd_id.update([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type in ["consumes","produces"]])
        pathways_ids = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_node.id] if rlt.type == "is_in_pathway"])
        total_pwy_id.update(pathways_ids)

    all_pwys = [node for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_pwy_id]
    all_cpds = [node for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_cpd_id]

    return [padmet_file, len(all_pwys), len(all_rxns), len(all_genes), len(all_cpds)] 

def orthology_result(padmet_file, padmet_names, padmet_name):
    ortholog_species_counts = dict.fromkeys(padmet_names)
    for species, count in ortholog_species_counts.items():
        if count is None:
            ortholog_species_counts[species] = 0

    padmetSpec = PadmetSpec(padmet_file)
    
    for node in padmetSpec.dicOfNode.values():
        if node.type == 'suppData':
            ortholog_species = node.id.split('FROM_')[1]
            ortholog_species_counts[ortholog_species] += 1

    df = pa.DataFrame([ortholog_species_counts],columns=ortholog_species_counts.keys(), index=[padmet_name])

    return df

if __name__ == "__main__":
    main()
