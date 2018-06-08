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
This script allows to combine rxn_creator.py and update_padmetSpec.py.
This script was created specially for AuReMe and the default metabolic network
reconstruction workflow.
If the file reaction_to_add_delete exist: calls update_padmetSpec
If the file new_reaction_data exist: calls rxn_creator.py
Update padmetSpec and create a new padmet (new_padmet) or overwritte the input

usage:
    manual_curation.py --padmetSpec=FILE --data=FILE [--padmetRef=FILE] [--output=FILE] [--tool=STR] [--category=STR] [-v]
    manual_curation.py --template_new_rxn --output=FILE
    manual_curation.py --template_add_delete --output=FILE

option:
    -h --help    Show help.
    --padmetSpec=FILE    pathname to the padmet to update
    --padmetRef=FILE    pathname of the padmet representing the reference database
    --data=FILE    pathname to the file used for rxn_creator.py
    --output=FILE    pathname to the ouput. if None. Overwritting padmetSpec
    -v    print info
"""
import docopt
import os
from time import time
from padmet.classes import PadmetSpec
from padmet.classes import PadmetRef
from padmet.classes import Relation
from padmet.utils.sbmlPlugin import parseGeneAssoc
import re
import csv

def main():
    args = docopt.docopt(__doc__)
    global tool, category, source
    data_file = args["--data"]
    output = args["--output"]
    if not output:
        output = args["--padmetSpec"]
    verbose = args["-v"]

    if data_file:
        filename = os.path.splitext(os.path.basename(data_file))[0]
        source = filename

    tool = args["--tool"]
    if tool:
        tool = tool.upper()
    if args["--category"]:
        category = args["--category"].upper()
    else:
        category = "MANUAL"

    if args["--template_new_rxn"]:
        template_new_rxn(output)
        return()
    elif args["--template_add_delete"]:
        template_add_delete(output)
        return()

    padmetSpec =  PadmetSpec(args["--padmetSpec"])
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None

        
    chronoDepart = time()
    with open(data_file, 'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        header = reader.next()
        if len(header) == 2:
            to_do = "rxn_creator"
        elif len(header) > 2:
            to_do = "add_delete_rxn"
        else:
            raise TypeError("Unable to read the file")
    if to_do == "rxn_creator":
        rxn_creator(data_file, padmetSpec, padmetRef, output, verbose)
    elif to_do == "add_delete_rxn":
        add_delete_rxn(data_file, padmetSpec, padmetRef, output, verbose)

    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    print "done in: ", chrono, "s !"

def rxn_creator(data_file, padmetSpec, padmetRef = None, output = None, verbose = False):
    dict_data = {}
    with open(data_file, 'r') as f:
        all_read = f.read()
    sep = csv.Sniffer().sniff(all_read).delimiter
    data = (line for line in all_read.splitlines() if len(line) != 0 and not line.startswith("#"))
    for line in data:
        #if len of value is 0 then TypeError raised
        try:
            attrib, value = line.split(sep)
        except TypeError:
            continue
        if attrib == "reaction_id":
            current_id = value
            dict_data[current_id] = {}
        else:
            try:
                dict_data[current_id][attrib] .append(value)
            except KeyError:
                dict_data[current_id][attrib] = [value]

    if verbose: print("%s reactions to add" %len(dict_data.keys()))
    for reaction_id, reaction_data in dict_data.iteritems():
        if verbose: print("check if the id %s is already used" %reaction_id)
        if reaction_id in padmetSpec.dicOfNode.keys():
            print("the id : %s is already associated to an other reaction in padmetSpec, choose an other" %reaction_id)
            continue
        if padmetRef is not None and reaction_id in padmetRef.dicOfNode.keys():
            print("the id : %s is already associated to an other reaction in padmetRef, choose an other" %reaction_id)
            continue

        if verbose: print("adding reaction %s" %reaction_id)
        reaction_rev = reaction_data["reversible"][0].lower()
        if reaction_rev.upper() == "TRUE":
            reaction_rev = "REVERSIBLE"
        elif reaction_rev.upper() == "FALSE":
            reaction_rev = "LEFT-TO-RIGHT"
        else:
            print("Please choose a value in ['true','false'] for the reversibility of the reaction: %s" %reaction_id)
            continue
        comment = reaction_data["comment"][0]
        node_misc = {"DIRECTION":[reaction_rev]}
        padmetSpec.createNode("reaction", reaction_id, node_misc)

        #reconstructionData:
        if tool:
            reconstructionData_id = reaction_id+"_reconstructionData_"+tool
            reconstructionData =  {"SOURCE": [source], "CATEGORY":[category], "TOOL":[tool], "COMMENT":[comment]}
            if reconstructionData_id in padmetSpec.dicOfNode.keys() and verbose:
                print("Warning: The reaction %s seems to be already added from the same source %s" %(reaction_id,tool))
        else:
            reconstructionData_id = reaction_id+"_reconstructionData_MANUAL"
            reconstructionData =  {"SOURCE": [source], "CATEGORY":["MANUAL"], "COMMENT":[comment]}
            if reconstructionData_id in padmetSpec.dicOfNode.keys() and verbose:
                print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %reaction_id)

        reconstructionData_rlt = Relation(reaction_id,"has_reconstructionData",reconstructionData_id)
        padmetSpec.createNode("reconstructionData", reconstructionData_id, reconstructionData, [reconstructionData_rlt])
    
        genes_assoc = reaction_data["linked_gene"][0]
        if genes_assoc:
            #suppData:
            if tool:
                suppData_id = reaction_id+"_SuppData_"+tool
                if suppData_id in padmetSpec.dicOfNode.keys() and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source %s" %(reaction_id,tool))
            else:
                suppData_id = reaction_id+"_SuppData_MANUAL"
                if suppData_id in padmetSpec.dicOfNode.keys() and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %reaction_id)
            suppData = {"GENE_ASSOCIATION":[genes_assoc]}
            #create the node suppData and the relation has_suppData
            suppData_rlt = Relation(reaction_id,"has_suppData",suppData_id)
            padmetSpec.createNode("suppData", suppData_id, suppData,[suppData_rlt])

            all_genes = parseGeneAssoc(genes_assoc)
            nbGenes = len(all_genes)
            if verbose: print("%s is linked to %s genes" %(reaction_id, nbGenes))
            for gene_id in all_genes:
               try:
                   #check if gene already in the padmet
                   padmetSpec.dicOfNode[gene_id]
               except KeyError:
                   padmetSpec.createNode("gene",gene_id)
               #check if rxn already linked to gene x
               try:
                   linked_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[reaction_id] if rlt.type == "is_linked_to"
                   and rlt.id_out == gene_id][0]
                   #rxn already linked to gene x, update misc
                   try:
                       linked_rlt.misc["SOURCE:ASSIGNMENT"].append(source)
                   except KeyError:
                       linked_rlt.misc["SOURCE:ASSIGNMENT"] = [source]
               #rxn not linked to gene x
               except IndexError:
                   linked_rlt = Relation(reaction_id, "is_linked_to", gene_id,{"SOURCE:ASSIGNMENT":[source]})
               padmetSpec._addRelation(linked_rlt)

        if verbose: print("check if all metabolites are already in the network")
        try:
            for reactant_data in reaction_data["reactant"]:
                stoechio, metabo_id, compart = reactant_data.split(":")
                stoechio = stoechio.replace(",",".") #in case comma for sep
                try:
                    padmetSpec.dicOfNode[metabo_id]
                except KeyError:
                    if verbose: print("%s not in the network" %metabo_id)
                    try:
                        if padmetRef is not None:
                            if verbose: print("Try to copy from dbref")
                            padmetSpec._copyNodeExtend(padmetRef, metabo_id)
                        else:
                            raise KeyError
                    except KeyError:
                        if padmetRef is not None and verbose:
                            print("%s not in the padmetRef" %metabo_id)
                        if verbose:
                            print("creating a new compound")
                        padmetSpec.createNode("compound", metabo_id)
                        if verbose: print("new compound created: id = %s" %metabo_id)
                rlt = Relation(reaction_id,"consumes",metabo_id)
                rlt.misc.update({"STOICHIOMETRY":[stoechio],"COMPARTMENT":[compart]})
                padmetSpec._addRelation(rlt)
        except KeyError:
            if verbose: print("No reactants defined")

        try:
            for product_data in reaction_data["product"]:
                stoechio, metabo_id, compart = product_data.split(":")
                stoechio = stoechio.replace(",",".") #in case comma for sep
                try:
                    padmetSpec.dicOfNode[metabo_id]
                except KeyError:
                    if verbose: print("%s not in the network" %metabo_id)
                    try:
                        if padmetRef is not None:
                            if verbose: print("Try to copy from dbref")
                            padmetSpec._copyNodeExtend(padmetRef, metabo_id)
                        else:
                            raise KeyError
                    except KeyError:
                        if padmetRef is not None and verbose:
                            print("%s not in the padmetRef" %metabo_id)
                        if verbose:
                            print("creating a new compound")
                        padmetSpec.createNode("compound", metabo_id)
                        print("new compound created: id = %s" % metabo_id)
                rlt = Relation(reaction_id,"produces",metabo_id)
                rlt.misc.update({"STOICHIOMETRY":[stoechio],"COMPARTMENT":[compart]})
                padmetSpec._addRelation(rlt)
        except KeyError:
            if verbose: print("No products defined")

    if verbose: print("Creating output: %s" % output)
    padmetSpec.generateFile(output)

def add_delete_rxn(data_file, padmetSpec, padmetRef = None, output = None, verbose = False):
    with open(data_file, 'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)        
        file_name = os.path.basename(data_file)
        file_name = os.path.splitext(file_name)[0]

        reader = csv.DictReader(csvfile, delimiter= dialect.delimiter)
        for row in reader:
            element_id, comment, action, genes_assoc = row["idRef"], row["Comment"], row["Action"], row.get("Genes",None)
            if action.upper() == "ADD":
                if padmetRef is None:
                    if verbose: print("No given padmetRef, unable to copy %s" %element_id)
                else:
                    if verbose: print("Adding: %s" %(element_id))
                    padmetSpec.copyNode(padmetRef, element_id)
                    
                    #reconstructionData:
                    if tool:
                        reconstructionData_id = element_id+"_reconstructionData_"+tool
                        reconstructionData =  {"SOURCE": [source], "CATEGORY":[category], "TOOL":[tool], "COMMENT":[comment]}
                        if reconstructionData_id in padmetSpec.dicOfNode.keys() and verbose:
                            print("Warning: The reaction %s seems to be already added from the same source %s" %(element_id,tool))
                    else:
                        reconstructionData_id = element_id+"_reconstructionData_MANUAL"
                        reconstructionData =  {"SOURCE": [source], "CATEGORY":["MANUAL"], "COMMENT":[comment]}
                        if reconstructionData_id in padmetSpec.dicOfNode.keys() and verbose:
                            print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %element_id)
                        
                    reconstructionData_rlt = Relation(element_id,"has_reconstructionData",reconstructionData_id)
                    padmetSpec.createNode("reconstructionData", reconstructionData_id, reconstructionData, [reconstructionData_rlt])
                    
                    if genes_assoc:
                        #suppData:
                        if tool:
                            suppData_id = element_id+"_SuppData_"+tool
                            if suppData_id in padmetSpec.dicOfNode.keys() and verbose:
                                print("Warning: The reaction %s seems to be already added from the same source %s" %(element_id,tool))
                        else:
                            suppData_id = element_id+"_SuppData_MANUAL"
                            if suppData_id in padmetSpec.dicOfNode.keys() and verbose:
                                print("Warning: The reaction %s seems to be already added from the same source 'MANUAL'" %element_id)
                        suppData = {"GENE_ASSOCIATION":[genes_assoc]}
                        #create the node suppData and the relation has_suppData
                        suppData_rlt = Relation(element_id,"has_suppData",suppData_id)
                        padmetSpec.createNode("suppData", suppData_id, suppData,[suppData_rlt])

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

            elif action.upper() == "DELETE":
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
        padmetSpec.generateFile(output)


def template_new_rxn(output):
    with open(output, 'w') as f:
        line = "\t".join(["reaction_id","my_rxn"])+"\n"
        f.write(line)
        line = "\t".join(["comment","reaction added for X reason"])+"\n"
        f.write(line)
        line = "\t".join(["reversible","false"])+"\n"
        f.write(line)
        line = "\t".join(["linked_gene","(gene_a or gene_b) and gene_c"])+"\n"
        f.write(line)
        line = "\t".join(["#reactant/product","#stoichio:compound_id:compart"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","1.0:compound_a:c"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","2.0:compound_b:c"])+"\n"
        f.write(line)
        line = "\t".join(["product","1.0:compound_c:c"])+"\n"
        f.write(line)
        f.write("\n")
        line = "\t".join(["reaction_id","my_rxn_2"])+"\n"
        f.write(line)
        line = "\t".join(["comment","reaction added for X reason"])+"\n"
        f.write(line)
        line = "\t".join(["reversible","true"])+"\n"
        f.write(line)
        line = "\t".join(["linked_gene",""])+"\n"
        f.write(line)
        line = "\t".join(["#reactant/product","#stoichio:compound_id:compart"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","1.0:compound_a:c"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","2.0:compound_d:c"])+"\n"
        f.write(line)
        line = "\t".join(["product","1.0:compound_c:c"])+"\n"
        f.write(line)

def template_add_delete(output):
    with open(output, 'w') as csvfile:
        fieldnames = ["idRef","Comment", "Action", "Genes"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow({"idRef": 'rxn_id_1', 'Comment': 'Reaction deleted for x reason', "Genes":"", "Action":"delete"})
        writer.writerow({"idRef": 'rxn_id_2', 'Comment': 'Reaction added for x reason', "Genes":"(gene1 and gene2)", "Action":"add"})
        writer.writerow({"idRef": 'rxn_id_3', 'Comment': 'Reaction added for x reason', "Genes":"", "Action":"add"})
 
if __name__ == "__main__":
    main()