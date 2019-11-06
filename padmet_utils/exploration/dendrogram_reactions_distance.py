#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
    Use reactions.csv file from compare_padmet.py to create a dendrogram using a Jaccard distance.
    
    From the matrix absence/presence of reactions in different species computes a Jaccard distance between these species.
    Apply a hierarchical clustering on these data with a complete linkage. Then create a dendrogram.
    Apply also intervene to create an upset graph on the data.


::

    usage:
        dendrogram_reactions_distance.py --reactions=FILE --output=FILE [--padmetRef=STR] [--pvclust] [--upset=INT] [-v]
    
    option:
        -h --help    Show help.
        -r --reactions=FILE    pathname of the file containing reactions in each species of the comparison.
        -o --output=FOLDER    path to the output folder.
        --pvclust    launch pvclust dendrogram using R
        --padmetRef=STR    path to the padmet Ref file
        -u --upset=INT    number of cluster in the upset graph.
        -v    verbose mode.

"""

import docopt
from padmet.utils.exploration import dendrogram_reactions_distance


def main():
    args = docopt.docopt(__doc__)
    reaction_pathname = args['--reactions']
    upset_cluster = int(args['--upset']) if args['--upset'] else None
    output_pathname = args['--output']
    padmet_ref_file = args['--padmetRef']
    pvclust = args['--pvclust']
    #verbose = args['-v']

    dendrogram_reactions_distance.reaction_figure_creation(reaction_pathname, output_pathname, upset_cluster, padmet_ref_file, pvclust)



if __name__ == '__main__':
    main()
