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
Require internet access !
Allows to extract the bigg database from the API to create a padmet.
1./ Get all reactions universal id from http://bigg.ucsd.edu/api/v2/universal/reactions
2./ Using async_list, extract all the informations for each reactions (compounds, stochio, name ...)
3./ Need to use sleep time to avoid to lose the server access.
4./ Because the direction fo the reaction is not set by default in bigg. 
We get all the models where the reaction is and the final direction will the one found
in more than 75%
5./ Also extract xrefs

usage:
    biggAPI_to_padmet.py

options:
    -h --help     Show help.
    --output=FILE    pathname of the padmet file to create
    --pwy_file=FILE   add kegg pathways from this pathways file 'pwy_id, pwy_name, x, rxn_id'.
    -v   print info.
"""
import json
from padmet.classes import Relation

import docopt
import requests
import grequests
from padmet.classes import PadmetSpec
from datetime import datetime
from time import time
import ast

def main():
    raw_to_padmet()
    
    
    
def raw_data():
    url_bigg = 'http://bigg.ucsd.edu/api/v2/'
    raw_data = requests.get(url_bigg + "universal/reactions").json()['results']
    all_reactions_ids = [rxn_dict['bigg_id'] for rxn_dict in raw_data if not rxn_dict['bigg_id'].startswith("BIOMASS")]
    print(("%s reactions to extract" %(len(all_reactions_ids))))

    with open("/home/maite/Documents/data/bigg/bigg_raw.txt", 'w') as f:
        count = 0
        for rxn_id in [i for i in all_reactions_ids if not i.startswith("BIOMASS")]:
            count += 1
            print(rxn_id)
            f.write("reaction: %s, %s/%s\n" %(rxn_id, count, len(all_reactions_ids)))
            rxn_response = requests.get(url_bigg + "universal/reactions/" +rxn_id)
            rxn_dict = rxn_response.json()
            rxn_metabolites = rxn_dict["metabolites"]
            if len(rxn_metabolites) > 1:
                json.dump(rxn_dict,f)
                f.write("\n")
                all_models_id = [i["bigg_id"] for i in rxn_dict["models_containing_reaction"]]
                async_list = []
                for model_id in all_models_id:
                    action_item = grequests.get(url_bigg + "models/"+ model_id +"/reactions/"+ rxn_id)
                    async_list.append(action_item)  
                models_responses = [r.json() for r in grequests.map(async_list)]
                all_lower_bound = [i["results"][0]["lower_bound"] for i in models_responses]
                ratio_not_rev = float(all_lower_bound.count(0))/float(len(all_lower_bound))
                f.write("Reaction not reversible in %s/%s model(s)\n" %(all_lower_bound.count(0), len(all_lower_bound)))
                if ratio_not_rev >= 0.75:
                    rxn_direction = "LEFT-TO-RIGHT"
                else:
                    rxn_direction = "REVERSIBLE"
                f.write(rxn_direction+"\n")

def raw_to_padmet():
    verbose = True
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
    dbNotes = {"PADMET":{"Creation":today_date, "version":"2.6"}, "DB_info":{"DB":"BIGG", "version":"2.0"}}
    padmet = PadmetSpec()
    if verbose: print("setting policy")
    padmet.setPolicy(policyInArray)
    if verbose: print("setting dbInfo")
    padmet.setInfo(dbNotes)
    list_of_relation = []

    with open("/home/maite/Documents/data/bigg/bigg_raw.txt", "r") as f:
        dataInArray = [line for line in f.read().splitlines() if not line.upper().startswith("REACTION")]

#    while index < len(dataInArray):
    index = 0
    raw_data = []
    while index < len(dataInArray):
        line = dataInArray[index]
        try:
            rxn_dict = ast.literal_eval(line)
            direction = dataInArray[index+1]
            rxn_dict.update({"direction":direction})
            raw_data.append(rxn_dict)
            index += 2
        except ValueError:
            index +=1

    count = 0
    for rxn_dict in raw_data:
        count += 1
        print("%s/%s" %(count, len(raw_data)))
        rxn_id = rxn_dict["bigg_id"]
        if rxn_id not in list(padmet.dicOfNode.keys()):
            rxn_metabolites = rxn_dict["metabolites"]
            rxn_name = rxn_dict["name"]
            rxn_direction = rxn_dict["direction"]
            padmet.createNode("reaction",rxn_id,{"COMMON-NAME":[rxn_name],"DIRECTION":[rxn_direction]})
            rxn_xrefs = rxn_dict["database_links"]
    
            xref_id = rxn_id+"_xrefs"
            xref_node = padmet.createNode("xref", xref_id)
            has_xref_rlt = Relation(rxn_id, "has_xref", xref_id)
            list_of_relation.append(has_xref_rlt)
    
            for db, k in list(rxn_xrefs.items()):
                _id = k[0]["id"]
                if db in list(xref_node.misc.keys()) and _id not in xref_node.misc[db]:
                    xref_node.misc[db].append(_id)
                else:
                    xref_node.misc[db] = [_id]
    
            for metabo_dict in rxn_metabolites:
                metabo_id = metabo_dict["bigg_id"]
                metabo_name = metabo_dict["name"]
                metabo_compart = metabo_dict["compartment_bigg_id"]
                metabo_stoich = metabo_dict["stoichiometry"]
                try:
                    padmet.dicOfNode[metabo_id]
                except KeyError:
                    padmet.createNode("compound",metabo_id,{"COMMON-NAME":[metabo_name]})
                if metabo_stoich < 0:
                    consumes_rlt = Relation(rxn_id,"consumes",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                    list_of_relation.append(consumes_rlt)
                else:
                    produces_rlt = Relation(rxn_id,"produces",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                    list_of_relation.append(produces_rlt)
        else: 
            if verbose: print("%s already in padmet" %rxn_id)
            continue                                        
    if verbose: print("Adding all relations")
    count = 0
    for rlt in list_of_relation:
        count += 1
        if verbose: print("relation %s/%s" %(count, len(list_of_relation)))
        try:
            padmet.dicOfRelationIn[rlt.id_in].append(rlt)
        except KeyError:
            padmet.dicOfRelationIn[rlt.id_in] = [rlt]
        try:
            padmet.dicOfRelationOut[rlt.id_out].append(rlt)
        except KeyError:
            padmet.dicOfRelationOut[rlt.id_out] = [rlt]
        
    pwy_file = "/home/maite/Documents/data/bigg/bigg_kegg_pwy.txt"
    add_kegg_pwy(pwy_file, padmet, verbose)
        

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
                if rxn_id in list(padmetRef.dicOfNode.keys()):
                    pwy_rlt = Relation(rxn_id,"is_in_pathway",pwy_id)
                    padmetRef._addRelation(pwy_rlt)
                else:
                    if verbose: print("%s in pwy %s but not in padmet" %(rxn_id, pwy_id))
    padmetRef.generateFile("/home/maite/Documents/data/bigg/bigg_v2.padmet")


if __name__ == "__main__":
    main()



