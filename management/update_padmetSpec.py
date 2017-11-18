# -*- coding: utf-8 -*-
"""
This file is part of padmet-utils.

padmet-utils is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

padmet-utils is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with padmet-utils. If not, see <http://www.gnu.org/licenses/>.

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
    --new_padmet=FILE    pathanme of the new padmetSpec after update
    --source=STR   source of the reactions to add [default: MANUAL]
    -v   print info
"""
from padmet.padmetSpec import PadmetSpec
from padmet.padmetRef import PadmetRef
from padmet.sbmlPlugin import parseGeneAssoc
from padmet.node import Node
from padmet.relation import Relation
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

    chronoDepart = time()
    with open(updateFile, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            element_id, comment, action, genes_assoc = row["idRef"], row["Comment"], row["Action"], row.get("Genes",None)
            if action == "add":
                if padmetRef is None:
                    if verbose: print("No given padmetRef, unable to copy %s" %element_id)
                else:
                    if verbose: print("Adding: %s" %(element_id))
                    padmetSpec.copyNode(padmetRef, element_id)
                    
                    #reconstructionData:
                    reconstructionData_id = element_id+"_reconstructionData_"+source
                    if reconstructionData_id in padmetSpec.dicOfNode.keys() and verbose:
                        print("Warning: The reaction %s seems to be already added from the same source %s" %(element_id, source))
                    reconstructionData =  {"SOURCE":[source], "COMMENT":[comment]}
                    reconstructionData_rlt = Relation(element_id,"has_reconstructionData",reconstructionData_id)
                    padmetSpec.createNode("reconstructionData", reconstructionData_id, reconstructionData, [reconstructionData_rlt])
                    if verbose: print("Creating reconstructionData %s" %reconstructionData_id)                
                    
                    if genes_assoc:
                        #suppData:
                        suppData_id = element_id+"_SuppData_"+source
                        if suppData_id in padmetSpec.dicOfNode.keys() and verbose:
                            print("Warning: The reaction %s seems to be already added from the same source %s" %(element_id, source))
                        suppData = {"GENE_ASSOCIATION":[genes_assoc]}
                        #create the node suppData and the relation has_suppData
                        suppData_rlt = Relation(element_id,"has_suppData",suppData_id)
                        padmetSpec.createNode("suppData", suppData_id, suppData,[suppData_rlt])
                        if verbose: print("Creating suppData %s" %suppData_id)                

                        all_genes = parseGeneAssoc(genes_assoc)
                        nbGenes = len(all_genes)
                        if verbose: print("%s is linked to %s genes" %(element_id, nbGenes))
                        for gene_id in all_genes:
                           try:
                               #check if gene already in the padmet
                               padmetSpec.dicOfNode[gene_id]
                           except KeyError:
                               padmetSpec.createNode("gene",gene_id)
                           #check if rxn already linked to gene x
                           try:
                               linked_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[element_id] if rlt.type == "is_linked_to"
                               and rlt.id_out == gene_id][0]
                               #rxn already linked to gene x, update misc
                               try:
                                   linked_rlt.misc["SOURCE:ASSIGNMENT"].append(source)
                               except KeyError:
                                   linked_rlt.misc["SOURCE:ASSIGNMENT"] = [source]
                           #rxn not linked to gene x
                           except IndexError:
                               linked_rlt = Relation(element_id, "is_linked_to", gene_id,{"SOURCE:ASSIGNMENT":[source]})
                           padmetSpec._addRelation(linked_rlt)

            elif action == "delete":
                if verbose:
                    print("deleting: %s" %(element_id))
                padmetSpec.delNode(element_id)
            elif action == "":
                print("Nothing to do for: %s" %(element_id))
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
        writer.writerow({"idRef": 'rxn_id_1', 'Comment': 'Reaction deleted for x reason', "Genes":"", "Action":"delete"})
        writer.writerow({"idRef": 'rxn_id_2', 'Comment': 'Reaction added for x reason', "Genes":"(gene1 and gene2)", "Action":"add"})
        writer.writerow({"idRef": 'rxn_id_3', 'Comment': 'Reaction added for x reason', "Genes":"", "Action":"add"})


if __name__ == "__main__":
    main()