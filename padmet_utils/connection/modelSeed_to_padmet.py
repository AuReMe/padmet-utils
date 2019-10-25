# -*- coding: utf-8 -*-
"""
Description:
    #TODO

::

    usage:
        modelSeed_to_padmet.py --output=FILE --rxn_file=FILE --pwy_file=FILE [-v]
    
    options:
        -h --help     Show help.
        --output=FILE    path of the padmet file to create
        --rxn_file=FILE   path to json file of modelSeed reactions
        --pwy_file=FILE   path to pathway reactions association from modelSeed
        -v   print info.
"""
from padmet.utils.connection import modelSeed_to_padmet
import docopt

def main():
    #parsing args
    args = docopt.docopt(__doc__)
    output = args["--output"]
    verbose = args["-v"]
    rxn_file = args["--rxn_file"]
    pwy_file = args["--pwy_file"]
    modelSeed_to_padmet.modelSeed_to_padmet(rxn_file, pwy_file, output, verbose)
     
if __name__ == "__main__":
    main()

