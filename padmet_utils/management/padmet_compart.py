# -*- coding: utf-8 -*-
"""
Description:
    For a given padmet file, check and update compartment.
    
    1./ Get all compartment with 1st usage
    
    2./ Remove a compartment with 2nd usage.
    Remove all reactions acting in the given compartment
    
    3./ change compartment id with 3rd usage

::
    
    usage:
        padmet_compart.py --padmet=FILE
        padmet_compart.py --padmet=FILE --remove=STR [--output=FILE] [-v]
        padmet_compart.py --padmet=FILE --old=STR --new=STR [--output=FILE] [-v]
    
    options:
        -h --help     Show help.
        --padmet=FILE    pathname of the padmet file
        --remove=STR    compartment id to remove
        --old=STR    compartment id to change to new id
        --new=STR    new compartment id
        --output=FILE    new padmet pathname, if none, overwritting the original padmet
        -v   print info
"""
from padmet.classes import PadmetSpec
from padmet.utils.management import padmet_compart
import docopt


def main():
    #parsing args
    args = docopt.docopt(__doc__)
    padmet_file = args["--padmet"]
    old_compart = args["--old"]
    new_compart = args["--new"]
    new_padmet = args["--output"]
    to_remove = args["--remove"]
    verbose = args["-v"]
    if new_padmet is None:
        new_padmet = args["--padmet"]
    padmet = PadmetSpec(padmet_file)

    if to_remove:
        padmet = padmet_compart.remove_compart(padmet, to_remove, verbose = False)
        padmet.generateFile(new_padmet)            

    elif old_compart and new_compart:
        padmet = padmet_compart.remplace_compart(padmet, old_compart, new_compart, verbose)
        padmet.generateFile(new_padmet)            
    else:
        print("List of compartments:")
        print(list(padmet.get_all_compart()))

if __name__ == "__main__":
    main()
