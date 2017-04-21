# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Convert a list of reactions to SBML file.

usage:
    reactions_to_sbml.py --padmetRef=FILE --reactions=FILE --output=FILE [-v]

option:
    -h --help     Show help.
    --padmetRef=FILE   pathname to the padmet used as reference (ex: metacyc.padmet)   
    --reactions=FILE    pathname to the file of reaction id. one id by line.
    --output=FILE    pathname to the sbml output.
    -v    print info
"""
from padmet.sbmlGenerator import reactions_to_SBML
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmetRef_file = args["--padmetRef"]
    reactions_file = args["--reactions"]     
    output = args["--output"]
    verbose = args["-v"]

    reactions_to_SBML(reactions_file, output, padmetRef_file, verbose=verbose)
            
if __name__ == "__main__":
    main()