# -*- coding: utf-8 -*-
"""

Description:
    Allows to merge 1-n padmet.
    1./ Update the 'init_padmet' with the 'to_add' padmet(s).
    to_add can be a file or a folder with only padmet files to add.
    
    padmetRef can be use to ensure data uniformization.

::
    
    usage:
        padmet_to_padmet.py --to_add=FILE/DIR --output=FILE [--padmetRef=FILE]  [-v]
    
    options:
        -h --help     Show help.
        --to_add=FILE/DIR    path to the padmet file to add (sep: ;) or path to folder of padmet files.
        --output=FILE   path to the new padmet file
        --padmetRef=FILE    path to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
        -v   print info
"""
from padmet.classes import PadmetRef
from padmet.utils.connection import padmet_to_padmet
import docopt

def main():
    args = docopt.docopt(__doc__)
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    output = args["--output"]
    verbose = args["-v"]
    to_add = args["--to_add"]
    padmet_to_padmet.padmet_to_padmet(to_add, output, padmetRef, verbose)

if __name__ == "__main__":
    main()
