#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Create a padmet stats file containing the number of pathways, reactions,
    genes and compounds inside the padmet.

    The input is a padmet file or a folder containing multiple padmets.
    
    Create a tsv file named padmet_stats.tsv where the script have been
    launched.

usage:
    padmet_stats.py --padmet=FILE

option:
    -h --help    Show help.
    -p --padmet=FILE    padmet file or folder containing padmet(s).

"""

import docopt
from padmet.utils.exploration import padmet_stats

def main():
    args = docopt.docopt(__doc__)

    padmet_file_folder = args['--padmet']

    padmet_stats.compute_stats(padmet_file_folder)



if __name__ == "__main__":
    main()
