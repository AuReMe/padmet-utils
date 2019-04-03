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

__all__ = ["remove_compart", "remplace_compart"]

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
        padmet = remove_compart(padmet, to_remove, verbose = False)
        padmet.generateFile(new_padmet)            

    elif old_compart and new_compart:
        padmet = remplace_compart(padmet, old_compart, new_compart, verbose)
        padmet.generateFile(new_padmet)            
    else:
        print("List of compartments:")
        print(list(padmet.get_all_compart()))

def remove_compart(padmet, to_remove, verbose = False):
    """
    Remove all reaction associated to a compound in the compartment to remove.

    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to udpate
    to_remove: str
        compartment id to remove, if many separate compartment id by ','
    verbose: bool
        if True print information

    Returns
    -------
    padmet.classes.PadmetSpec:
        New padmet after removing compartment(s)
    """
    if "," in to_remove:
        to_remove = to_remove.split(",")
    else:
        to_remove = [to_remove]
    for compart in to_remove:
        if verbose:
            print("removing reaction(s) in compartment %s" %compart)
        padmet.delCompart(compart, verbose)
    return padmet

def remplace_compart(padmet, old_compart, new_compart, verbose = False):    
    """
    Remplace compartment 'old_compart' by 'new_compart'.

    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to udpate
    old_comaprt: str
        compartment id to remplace
    new_compart: str
        new compartment id
    verbose: bool
        if True print information

    Returns
    -------
    padmet.classes.PadmetSpec:
        New padmet after remplacing compartment

    """
    if verbose:
        print("Converting compartment %s to %s" %(old_compart, new_compart))
    padmet.change_compart(old_compart, new_compart, verbose)
    return padmet

if __name__ == "__main__":
    main()
