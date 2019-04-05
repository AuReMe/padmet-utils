# -*- coding: utf-8 -*-
"""

Description:
    Allows to merge 1-n padmet.
    1./ Update the 'init_padmet' with the 'to_add' padmet(s).
    to_add can be a file or a folder with only padmet files to add.
    
    padmetRef can be use to ensure data uniformization.

::
    
    usage:
        padmet_to_padmet.py --to_add=FILE/DIR --output=FILE [--init_padmet=FILE] [--padmetRef=FILE]  [-v]
    
    options:
        -h --help     Show help.
        --to_add=FILE/DIR    path to the padmet file to add or path to folder of padmet files.
        --output=FILE   path to the new padmet file
        --init_padmet=FILE    path to the padmet file to update.
        --padmetRef=FILE    path to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
        -v   print info
"""
from padmet.classes import PadmetSpec, PadmetRef
from datetime import datetime
import os
import docopt

def main():
    args = docopt.docopt(__doc__)
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    if args["--init_padmet"]:
        padmet_init = PadmetSpec(args["--init_padmet"])
    else:
        padmet_init = None
    output = args["--output"]
    verbose = args["-v"]
    to_add = args["--to_add"]
    padmet_to_padmet(to_add, output, padmet_init, padmetRef, verbose)

def padmet_to_padmet(to_add, output, padmet_init=None, padmetRef=None, verbose=False):
    """
    
    """
    if os.path.isdir(to_add):
        padmet_type = "path"
    elif os.path.isfile(to_add):
        padmet_type = "file"
    else:
        raise TypeError("%s is not a dir or a file" %(to_add))

    if padmet_type == "path":
        path = to_add
        all_files = [i for i in next(os.walk(path))[2] if not i.startswith(".~lock")]
        padmetFiles = [os.path.join(path, i) for i in all_files if i.endswith(".padmet")]
        if len(padmetFiles) == 0:
            raise IOError("No padmet found in %s" %path)
    else:
        padmetFiles = list([to_add])


    if not padmet_init:
        now = datetime.now()
        today_date = now.strftime("%Y-%m-%d")

        padmet_init = PadmetSpec()
        if padmetRef:
            padmet_init.setInfo(padmetRef)
            padmet_init.info["PADMET"]["creation"] = today_date
            padmet_init.setPolicy(padmetRef)
        else:
            POLICY_IN_ARRAY = [['class','is_a_class','class'], ['class','has_name','name'], ['class','has_xref','xref'], ['class','has_suppData','suppData'],
                            ['compound','is_a_class','class'], ['compound','has_name','name'], ['compound','has_xref','xref'], ['compound','has_suppData','suppData'],
                            ['gene','is_a_class','class'], ['gene','has_name','name'], ['gene','has_xref','xref'], ['gene','has_suppData','suppData'], ['gene','codes_for','protein'],
                            ['pathway','is_a_class','class'], ['pathway','has_name','name'], ['pathway','has_xref','xref'], ['pathway','is_in_pathway','pathway'], 
                            ['protein','is_a_class','class'], ['protein','has_name','name'], ['protein','has_xref','xref'], ['protein','has_suppData','suppData'], ['protein','catalyses','reaction'],
                            ['protein','is_in_species','class'], 
                            ['reaction','is_a_class','class'], ['reaction','has_name','name'], ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','has_reconstructionData','reconstructionData'], ['reaction','is_in_pathway','pathway'],  
                            ['reaction','consumes','class','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','class','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                            ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                            ['reaction','consumes','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                            ['reaction','is_linked_to','gene','SOURCE:ASSIGNMENT','X:Y']]
            dbNotes = {"PADMET":{"creation":today_date,"version":"2.6"},"DB_info":{"DB":'NA',"version":'NA'}}
            padmet_init.setInfo(dbNotes)
            padmet_init.setPolicy(POLICY_IN_ARRAY)


    for padmet_update_file in padmetFiles:
        if verbose:
            print("Updating %s from %s" %(os.path.basename(padmet_init),os.path.basename(padmet_update_file)))
        padmet_update = PadmetSpec(padmet_update_file)
        padmet_init.updateFromPadmet(padmet_update)

    if verbose:
        print("Generated file: %s" %output)
    padmet_init.generateFile(output)
        

if __name__ == "__main__":
    main()
