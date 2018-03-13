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
#TODO
usage:
    modelSeed_to_padmet.py --output=FILE --rxn_file=FILE --pwy_file=FILE [-v]

options:
    -h --help     Show help.
    --output=FILE    pathname of the padmet file to create
    --rxn_file=FILE   json file of modelSeed reactions
    --pwy_file=FILE   pathway reactions association from modelSeed
    -v   print info.
"""
from padmet.node import Node
from padmet.relation import Relation
from padmet.padmetRef import PadmetRef
from datetime import datetime
from time import time

import csv
import json
import docopt
import re

def main():
    global list_of_relation
    chronoDepart = time()
    #parsing args
    args = docopt.docopt(__doc__)
    output = args["--output"]
    verbose = args["-v"]
    rxn_file = args["--rxn_file"]
    pwy_file = args["--pwy_file"]
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")
    #print(verbose,today_date,version, output, classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file)
    policyInArray = [['compound','has_name','name'], ['compound','has_xref','xref'], ['compound','has_suppData','suppData'],
                    ['gene','has_name','name'], ['gene','has_xref','xref'], ['gene','has_suppData','suppData'], ['gene','codes_for','protein'],
                    ['pathway','has_name','name'], ['pathway','has_xref','xref'], ['pathway','is_in_pathway','pathway'], 
                    ['protein','has_name','name'], ['protein','has_xref','xref'], ['protein','has_suppData','suppData'], ['protein','catalyses','reaction'],
                    ['reaction','has_name','name'], ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','has_reconstructionData','reconstructionData'], ['reaction','is_in_pathway','pathway'],  
                    ['reaction','consumes','class','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','class','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                    ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                    ['reaction','consumes','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                    ['reaction','is_linked_to','gene','SOURCE:ASSIGNMENT','X:Y']]
    dbNotes = {"PADMET":{"Creation":today_date, "version":"2.6"}, "DB_info":{"DB":"MODELSEED", "version":"1.0"}}
    padmetRef = PadmetRef()
    if verbose: print("setting policy")
    padmetRef.setPolicy(policyInArray)
    if verbose: print("setting dbInfo")
    padmetRef.setInfo(dbNotes)
    list_of_relation = []
    
    rxn_data = json.load(open(rxn_file))
    #remove biomass rxn:
    rxn_data.pop("rxn12985")
    if verbose: print("updating padmet")
    count = 0
    for rxn_id, rxn_dict in rxn_data.items():
        count += 1
        if verbose: print("reaction: %s, %s/%s" %(rxn_id, count, len(rxn_data)))
        try:
            if not rxn_dict["compound_ids"]:
                raise KeyError
        except KeyError:
            print(rxn_id)
            continue
        if rxn_id not in padmetRef.dicOfNode.keys():
            if rxn_dict["reversibility"] == ">":
                rxn_direction = "LEFT-TO-RIGHT"
            else:
                rxn_direction = "REVERSIBLE"
            rxn_name = rxn_dict["name"]
            padmetRef.createNode("reaction",rxn_id,{"COMMON_NAME":[rxn_name],"DIRECTION":[rxn_direction]})
                            
            rxn_metabolites = rxn_dict["stoichiometry"].split(";")
    
            for metabo_data in rxn_metabolites:
                metabo_data = metabo_data.replace("???","\"")
                try:
                    metabo_temp, metabo_name = metabo_data.split("\"")[:2]
                    metabo_stoich, metabo_id, metabo_compart = metabo_temp.split(":")[:3]
                except ValueError:
                    metabo_stoich, metabo_id, metabo_compart, metabo_name = metabo_data.split(":")[:4]
                    
                metabo_stoich = float(metabo_stoich)
                #from modelSeed github
                if metabo_compart == "0":
                    metabo_compart = "c"
                elif metabo_compart == "1":
                    metabo_compart = "e"
                elif metabo_compart == "2":
                    metabo_compart = "p"
                try:
                    padmetRef.dicOfNode[metabo_id]
                except KeyError:
                    padmetRef.createNode("compound",metabo_id,{"COMMON_NAME":[metabo_name]})
                if metabo_stoich < 0:
                    consumes_rlt = Relation(rxn_id,"consumes",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                    list_of_relation.append(consumes_rlt)
                else:
                    produces_rlt = Relation(rxn_id,"produces",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                    list_of_relation.append(produces_rlt)
        else: 
            if verbose: print("%s already in padmet" %rxn_id)
            continue
    with open(pwy_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t' )
        pwy_raw_data = [row for row in reader]
    for pwy_raw in pwy_raw_data:
        pwy_id = pwy_raw["Source ID"]
        pwy_names = [pwy_raw["Name"],pwy_raw["Aliases"]]
        rxn_ids = pwy_raw["Reactions"].split("|")
        try:
            padmetRef.dicOfNode[pwy_id]
        except KeyError:
            padmetRef.createNode("pathway",pwy_id,{"COMMON_NAME":pwy_names})
        for rxn_id in rxn_ids:
            pwy_rlt = Relation(rxn_id,"is_in_pathway",pwy_id)
            list_of_relation.append(pwy_rlt)
        
        
    if verbose: print("Adding all relations")
    count = 0
    for rlt in list_of_relation:
        count += 1
        if verbose: print("relation %s/%s" %(count, len(list_of_relation)))
        try:
            padmetRef.dicOfRelationIn[rlt.id_in].append(rlt)
        except KeyError:
            padmetRef.dicOfRelationIn[rlt.id_in] = [rlt]
        try:
            padmetRef.dicOfRelationOut[rlt.id_out].append(rlt)
        except KeyError:
            padmetRef.dicOfRelationOut[rlt.id_out] = [rlt]
    
    """
    if pwy_file:
        add_kegg_pwy(pwy_file, padmetRef, verbose)
    """
    if verbose: print("Generating file: %s" %output)
    padmetRef.generateFile(output)
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print "done in: ", chrono, "s !"

def add_kegg_pwy(pwy_file, padmetRef, verbose = False):
    global list_of_relation
    with open(pwy_file, 'r') as f:
        for data in [line.split("\t") for line in f.read().splitlines()][1:]:
            pwy_id, name, ec, rxn_id = data
            try:
                pwy_node = padmetRef.dicOfNode[pwy_id]
            except KeyError:
                pwy_node = padmetRef.createNode("pathway", pwy_id)
            if name:
                try:
                    pwy_node.misc["COMMON_NAME"].append(name)
                except KeyError:
                    pwy_node.misc["COMMON_NAME"] = [name]
            if rxn_id:
                if rxn_id in padmetRef.dicOfNode.keys():
                    pwy_rlt = Relation(rxn_id,"is_in_pathway",pwy_id)
                    padmetRef._addRelation(pwy_rlt)
                else:
                    if verbose: print("%s in pwy %s but not in padmet" %(rxn_id, pwy_id))
    padmetRef.generateFile("/home/maite/Documents/data/bigg/bigg_v2.padmet")
            
    
if __name__ == "__main__":
    main()

