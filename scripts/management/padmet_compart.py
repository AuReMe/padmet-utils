# -*- coding: utf-8 -*-
"""
Description:
    #TODO

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
    padmetSpec = PadmetSpec(padmet_file)

    if to_remove:
        if "," in to_remove:
            to_remove = to_remove.split(",")
        else:
            to_remove = [to_remove]
        for compart in to_remove:
            if verbose:
                print("removing reaction(s) in compartment %s" %compart)
            padmetSpec.delCompart(compart, verbose)
        padmetSpec.generateFile(new_padmet)            
    elif old_compart and new_compart:
        if verbose:
            print("Converting compartment %s to %s" %(old_compart, new_compart))
        padmetSpec.change_compart(old_compart, new_compart, verbose)
        padmetSpec.generateFile(new_padmet)            
    else:
        print("List of compartments:")
        print(list(padmetSpec.get_all_compart()))

if __name__ == "__main__":
    main()
