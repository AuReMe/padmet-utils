# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Allows to add or delete reactions/metabolites/pathways based on the given id.
From the updateFile (obtained with 2nd usage) the script extract the id (col1),
why adding/deleting this element (comment col2) and the action (col3)
THE ACTION IS IN: ['add','delete','']
When add: the element is added from padmetRef, and the COMMENT and SOURCE is stored in the node of
the element in padmetSpec.
When delete: the element is deleted from padmetSpec.
Whene '': Nothing to do, but allow to manage the output of enhancedMenecoOutput.py

usage:
    update_padmetSpec.py --padmetSpec=FILE --updateFile=FILE [--padmetRef=FILE] [--source=STR] [--new_padmet=FILE] [-v]
    update_padmetSpec.py --getTemplate --output=FILE

option:
    -h --help     Show help.
    --padmetSpec=FILE    pathanme of the padmet file to update
    --padmetRef=FILE  pathname of the padmet used as database of reactions
    --updateFile=FILE    pathname of the file containing elements ids to add or delete
    --new_padmet=FILE    pathanme of the new padmetSpec after update [default: --padmetSpec]
    --source=STR   source of the reactions to add [default: manual]
    -v   print info
"""
from padmet.padmetSpec import PadmetSpec
from padmet.padmetRef import PadmetRef
from time import time
import csv
import docopt

def main():
    args = docopt.docopt(__doc__)
    
    get_template = args["--getTemplate"]
    if get_template:
        output = args["--output"]
        getTemplate(output)
        return
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    if args["--padmetRef"] is not None:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    updateFile = args["--updateFile"]
    source = args["--source"]
    verbose = args["-v"]
    new_padmet = args["--new_padmet"]
    if new_padmet is None:
        new_padmet = args["--padmetSpec"]
    with open(updateFile,'r') as f:
        file_in_array = f.read().splitlines()
        data = dict([(line.split("\t")[0],line.split("\t")[-1]) for line in file_in_array
        if line.split("\t")[0] != "idRef"])
    
    nb_elements = len(data.keys())
    count = 0

    chronoDepart = time()    

    with open(updateFile, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            element_id, comment, action = row["idRef"], row["Comment"], row["Action"]
            count += 1
            if action == "add":
                if padmetRef is None:
                    if verbose: print("No given padmetRef, unable to copy %s" %element_id)
                else:
                    if verbose: print("copy: %s %s/%s" %(element_id, count, nb_elements))
                    padmetSpec.copyNode(padmetRef, element_id)
                    try:
                        padmetSpec.dicOfNode[element_id].misc["SOURCE"].append(source)
                    except KeyError:
                        padmetSpec.dicOfNode[element_id].misc["SOURCE"] = [source]
                    if len(comment) != 0:
                        try:
                            padmetSpec.dicOfNode[element_id].misc["COMMENT"].append(comment)
                        except KeyError:
                            padmetSpec.dicOfNode[element_id].misc["COMMENT"] = [comment]
            elif action == "delete":
                if verbose:
                    print("delete: %s %s/%s" %(element_id, count, nb_elements))
                padmetSpec.delNode(element_id)
            elif action == "":
                print("Nothing to do for: %s %s/%s" %(element_id, count, nb_elements))
                pass
            else:
                print("Action: %s unknown for %s" %(action, element_id))
                print("action must be = 'add' or 'delete' or ''")
                exit()
        padmetSpec.generateFile(new_padmet)
    
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose:
        print "done in: ", chrono, "s !"

def getTemplate(output):
    with open(output, 'w') as csvfile:
        fieldnames = ["idRef","Comment", "Action"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow({"idRef": 'rxn_id_1', 'Comment': 'Reaction deleted for x reason', "Action":"delete"})
        writer.writerow({"idRef": 'rxn_id_2', 'Comment': 'Reaction added for x reason', "Action":"add"})


if __name__ == "__main__":
    main()