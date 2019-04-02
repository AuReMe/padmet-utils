#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Use reactions.csv file from compare_padmet.py to create a dendrogram using a Jaccard distance.

From the matrix absence/presence of reactions in different species computes a Jaccard distance between these species.
Then create a dendrogram using these distances.

usage:
    dendrogram_reactions_distance.py --reactions=FILE --output=FILE [-v]

option:
    -h --help    Show help.
    -r --reactions=FILE    pathname of the file containing reactions in each species of the comparison.
    -o --output=FOLDER    path to the output folder.
    -v    verbose mode.

"""

import docopt
import itertools
import pandas as pa
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import subprocess
sns.set_style("white")
sns.set('poster', rc={'figure.figsize':(100,80)}, font_scale=4)

from collections import defaultdict
from scipy.cluster.hierarchy import dendrogram, fcluster
from scipy.spatial import distance
from sklearn.metrics.pairwise import pairwise_distances
from fastcluster import linkage
from scipy.spatial.distance import squareform, pdist

def main():
    global verbose

    args = docopt.docopt(__doc__)
    reaction_pathname = args["--reactions"]
    output_pathname = args["--output"]

    if args['-v']:
        verbose = args["-v"]
    else:
        verbose = None

    reaction_figure_creation(reaction_pathname, output_pathname)

def reaction_figure_creation(reaction_file, output_folder):
    # Check if output_folder exists, if not create it.
    output_folder_data = output_folder + '/data'
    output_folder_data_intersect = output_folder + '/data/intersect'
    output_folder_data_unique = output_folder + '/data/unique'
    output_folder_upset = output_folder + '/upset_graph'
    temp_data_folder = output_folder + '/upset_graph/temp_data/'

    folders = [output_folder, output_folder_data, output_folder_data_intersect,
                output_folder_data_unique, output_folder_upset, temp_data_folder]

    for folder in folders:
        if not os.path.isdir("{0}".format(folder)):
            os.mkdir("{0}".format(folder))

    path_to_intervene = 'intervene'

    # Read the reactions file with pandas.
    all_reactions_dataframe = pa.read_csv(reaction_file, sep='\t')
    # Keep column containing absence-presence of reactions.
    # (columns with (sep=;) are column with gene name linked to reactions)
    # (columns with _formula contain the reaction formula)
    columns = [column for column in all_reactions_dataframe.columns if '(sep=;)' not in column]
    columns = [column for column in columns if '_formula' not in column]
    reactions_dataframe = all_reactions_dataframe[columns].copy()

    reactions_dataframe.set_index('reaction', inplace=True)

    # Translate 'present'/(nan) data into a True/False absence-presence matrix.
    for column in reactions_dataframe.columns.tolist():
        reactions_dataframe[column] = [True if data == "present" else False for data in reactions_dataframe[column]]

    # Transpose the matrix to have species as index and reactions as columns.
    absence_presence_matrix = reactions_dataframe.transpose()

    # Compute a distance matrix using the Jaccard distance between species. Then condense it.
    distance_matrix_jaccard = distance.squareform(pairwise_distances(absence_presence_matrix, metric="jaccard"))

    # Hierarchical clustering on the condensed distance matrix.
    linkage_matrix = linkage(distance_matrix_jaccard, method="average")

    # Draw a dendrogram of the clustering.
    reaction_dendrogram = dendrogram(linkage_matrix, labels=absence_presence_matrix.index, leaf_font_size=40)
    plt.savefig(output_folder+'/reaction_dendrogram.png')

    # Extract species in each clusteR.
    k = len(set(reaction_dendrogram['color_list']))
    results = fcluster(linkage_matrix, k, criterion='maxclust')

    species = absence_presence_matrix.index.tolist()

    cluster_species = dict(zip(species, results))
    cluster_classes = defaultdict(list)

    for key, value in cluster_species.items():
        cluster_classes[value].append(key)

    # Extract reactions in each cluster.
    cluster_reactions = {}
    for cluster in cluster_classes:
        reactions_temp = []
        for species in cluster_classes[cluster]:
            species_reactions_dataframe = reactions_dataframe[reactions_dataframe[species] == True]
            reactions_temp.extend(species_reactions_dataframe.index.tolist())
        cluster_reactions[cluster] = set(reactions_temp)

    all_reactions = [reactions for reactions in cluster_reactions.values()]

    cluster_intersections = set.intersection(*all_reactions)
    # Create file containing the intersection of the reactions for all cluster.
    df = pa.DataFrame({'all_species': list(cluster_intersections)})
    df.set_index('all_species', inplace=True)
    all_reactions_dataframe.set_index('reaction', inplace=True)
    gene_assoc_columns = [column for column in all_reactions_dataframe.columns if '(sep=;)' in column]
    gene_assoc_reactions = all_reactions_dataframe[gene_assoc_columns]
    df = df.join(gene_assoc_reactions)
    df.to_csv(output_folder_data+'/'+'all_species.tsv', sep='\t', index=True)

    cluster_subintersection = {}
    cluster_subintersection_name = {}
    # Extract intersection between clusters.
    for cluster_number in reversed(range(len(cluster_reactions))):
        if cluster_number != 0 and cluster_number != 1:
            for set_list in itertools.combinations(cluster_reactions, cluster_number):
                tmp_reactions = [cluster_reactions[cluster] for cluster in set_list]
                cltemp = set.intersection(*tmp_reactions)
                intersection_temp = cltemp - cluster_intersections
                cluster_subintersection[set_list] = intersection_temp
                cluster_species_name = '&'.join(['_'.join(cluster_classes[cluster]) for cluster in set_list])
                cluster_subintersection_name[cluster_species_name] = list(intersection_temp)
                # Create a file containing intersection between each cluster.
                df = pa.DataFrame({'reaction': list(intersection_temp)})
                df.set_index('reaction', inplace=True)
                gene_assoc_columns = [column for column in all_reactions_dataframe.columns if '(sep=;)' in column]
                column_species = [ species for cluster in set_list for species in cluster_classes[cluster]]
                temp_gene_assoc_columns = [gene_assoc for gene_assoc in gene_assoc_columns if gene_assoc.split('_genes_assoc')[0] in column_species]
                gene_assoc_reactions = all_reactions_dataframe[temp_gene_assoc_columns]
                df = df.join(gene_assoc_reactions)
                df.to_csv(output_folder_data_intersect+'/'+cluster_species_name+'_intersect.tsv', sep='\t', index=True)

    # Create reactions which intersect for each cluster.
    cluster_subsubintersection = {}
    for cluster in cluster_classes:
        species_intersections = []
        for set_list in cluster_subintersection:
            if cluster in set_list:
                species_intersections.append(cluster_subintersection[set_list])
        cluster_subsubintersection[cluster] = set([j for i in species_intersections for j in i])

    # Extract reactions unique for each cluster.
    cluster_unique = {}
    for cluster in cluster_classes:
        cluster_unique[cluster] = cluster_reactions[cluster]-cluster_intersections-cluster_subsubintersection[cluster]
        # Create a file containing reactions unique for each cluster.
        df = pa.DataFrame({'reaction': list(cluster_unique[cluster])})
        df.set_index('reaction', inplace=True)
        gene_assoc_columns = [column for column in all_reactions_dataframe.columns if '(sep=;)' in column]
        column_species = [ species for species in cluster_classes[cluster]]
        temp_gene_assoc_columns = [gene_assoc for gene_assoc in gene_assoc_columns if gene_assoc.split('_genes_assoc')[0] in column_species]
        gene_assoc_reactions = all_reactions_dataframe[temp_gene_assoc_columns]
        df = df.join(gene_assoc_reactions)
        df.to_csv(output_folder_data_unique+'/'+'_'.join(cluster_classes[cluster])+'_unique.tsv', sep='\t', index=True)

    # Create data for creating upset graph using intervene.
    for cluster in cluster_classes:
        df = pa.DataFrame({'_'.join(cluster_classes[cluster]): list(cluster_reactions[cluster])})
        df.to_csv(temp_data_folder+'/'+'_'.join(cluster_classes[cluster])+'.tsv', sep='\t', index=None, header=None)

    cmd = '{0} upset -i  {1}/*.tsv --type list -o {2} --figtype svg'.format(path_to_intervene, temp_data_folder, output_folder_upset)
    if verbose:
        subprocess.call(cmd, shell=True)
    else:
        FNULL = open(os.devnull, 'w')
        subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

main()
