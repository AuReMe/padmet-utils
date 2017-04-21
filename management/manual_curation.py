# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
This script allows to combine rxn_creator.py and update_padmetSpec.py.
This script was created specially for AuReMe and the default metabolic network
reconstruction workflow.
If the file reaction_to_add_delete exist: calls update_padmetSpec
If the file new_reaction_data exist: calls rxn_creator.py
Update padmetSpec and create a new padmet (new_padmet) or overwritte the input

usage:
    manual_curation.py --padmetSpec=FILE --reaction_to_add_delete=FILE  [--new_reaction_data=FILE] [--padmetRef=FILE] [--new_padmet=FILE] [-v]
    manual_curation.py --padmetSpec=FILE [--reaction_to_add_delete=FILE]  --new_reaction_data=FILE [--padmetRef=FILE] [--new_padmet=FILE] [-v]

option:
    -h --help    Show help.
    --padmetSpec=FILE    pathname to the padmet to update
    --padmetRef=FILE    pathname of the padmet representing the reference database
    --reaction_to_add_delete=FILE    pathname to the file used for update_padmetSpec.py
    --reaction_creator=FILE    pathname to the file used for rxn_creator.py
    --new_padmet=FILE    pathname to the ouput. if None. Overwritting padmetSpec
    -v    print info
"""
import docopt
import os
from time import time


def main():
    #set the path to the scripts based on the current path.
    dir_path_reaction_creator = os.path.dirname(os.path.realpath(__file__))+"/rxn_creator.py"
    dir_path_update_padmetSpec = os.path.dirname(os.path.realpath(__file__))+"/misc/update_padmetSpec.py"
    args = docopt.docopt(__doc__)
    chronoDepart = time()
    reaction_to_add_delete = args["--reaction_to_add_delete"]
    new_reaction_data = args["--new_reaction_data"]
    padmetSpec = args["--padmetSpec"]
    padmetRef = args["--padmetRef"]
    new_padmet = args["--new_padmet"]
    verbose = args["-v"]
    if os.path.isfile(new_reaction_data):
        rxn_creator_cmd = "python "+dir_path_reaction_creator+" --padmetSpec="+padmetSpec+" --padmetRef="+padmetRef+" --new_reaction_data="+new_reaction_data+\
        " --new_padmet="+new_padmet
        if verbose:
            rxn_creator_cmd += " -v"
        os.system(rxn_creator_cmd)
        if os.path.isfile(reaction_to_add_delete):
            update_padmet_cmd = "python "+dir_path_update_padmetSpec+" --padmetSpec="+new_padmet+" --padmetRef="+padmetRef+" --updateFile="+reaction_to_add_delete+\
            " --new_padmet="+new_padmet
            if verbose:
                update_padmet_cmd += " -v"
            os.system(update_padmet_cmd)
    elif os.path.isfile(reaction_to_add_delete):
        update_padmet_cmd = "python "+dir_path_update_padmetSpec+" --padmetSpec="+padmetSpec+" --padmetRef="+padmetRef+" --updateFile="+reaction_to_add_delete+\
        " --new_padmet="+new_padmet
        if verbose:
            update_padmet_cmd += " -v"
        os.system(update_padmet_cmd)

    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    print "done in: ", chrono, "s !"


    
if __name__ == "__main__":
    main()