# -*- coding: utf-8 -*-
"""
Description:
    Require internet access !

    Allows to extract the bigg database from the API to create a padmet.

    1./ Get all reactions universal id from http://bigg.ucsd.edu/api/v2/universal/reactions, escape reactions of biomass.

    2./ Using async_list, extract all the informations for each reactions (compounds, stochio, name ...)

    3./ Need to use sleep time to avoid to lose the server access.

    4./ Because the direction fo the reaction is not set by default in bigg. 
    We get all the models where the reaction is and the final direction will the one found
    in more than 75%

    5./ Also extract xrefs

::

    usage:
        biggAPI_to_padmet.py --output=FILE [--pwy_file=FILE] [-v]
    
    options:
        -h --help     Show help.
        --output=FILE    path to output, the padmet file.
        --pwy_file=FILE   add kegg pathways from pathways file, line:'pwy_id, pwy_name, x, rxn_id'.
        -v   print info.
"""

import docopt
from padmet.utils.connection import biggAPI_to_padmet

def main():
    global list_of_relation
    #parsing args
    args = docopt.docopt(__doc__)
    output = args["--output"]
    verbose = args["-v"]
    pwy_file = args["--pwy_file"]
    biggAPI_to_padmet.biggAPI_to_padmet(output, pwy_file, verbose)
            
    
if __name__ == "__main__":
    main()

