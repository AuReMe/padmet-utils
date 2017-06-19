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
    biggAPI_to_padmet.py --output=FILE [-v]

options:
    -h --help     Show help.
    --output=FILE    pathname of the padmet file to create
    -v   print info.
"""
from padmet.node import Node
from padmet.relation import Relation
from padmet.padmetRef import PadmetRef
from datetime import datetime
from time import time

import docopt
import requests
import grequests

def main():
    chronoDepart = time()
    #parsing args
    args = docopt.docopt(__doc__)
    output = args["--output"]
    verbose = args["-v"]
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")
    #print(verbose,today_date,version, output, classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file)
    policyInArray = [['compound','has_xref','xref'], ['compound','has_suppData','suppData'], 
                    ['gene','has_xref','xref'], ['gene','has_suppData','suppData'],
                    ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','is_linked_to','gene'],
                    ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y']] 
    dbNotes = {"PADMET":{"Creation":today_date, "version":"2.4"}, "DB_info":{"DB":"BIGG", "version":"2.0"}}
    padmetRef = PadmetRef()
    if verbose: print("setting policy")
    padmetRef.setPolicy(policyInArray)
    if verbose: print("setting dbInfo")
    padmetRef.setInfo(dbNotes)
    new_id_prefix = "META"
    meta_max_id = 0
    list_of_relation = []
    if verbose: print("Getting all reactions ids")
    url_bigg = 'http://bigg.ucsd.edu/api/v2/'
    raw_data = requests.get(url_bigg + "universal/reactions").json()['results']
    all_reactions_ids = [rxn_dict['bigg_id'] for rxn_dict in raw_data]
    if verbose: print("%s reactions to extract" %(len(all_reactions_ids)))

    """
    if verbose: print("Extracting informations... Wait")
    step = 100
    rxn_lower_index = -(step)
    rxn_upper_index = 0
    rxn_responses = []
    all_range = len(all_reactions_ids)/step

    for i in range(all_range):
        async_list = []
        rxn_lower_index += step
        rxn_upper_index += step

        for rxn_id in all_reactions_ids[rxn_lower_index:rxn_upper_index]:
            action_item = grequests.get(url_bigg + "universal/reactions/" +rxn_id)
            async_list.append(action_item) 
        new_responses = [r.json() for r in grequests.map(async_list)]
        rxn_responses += new_responses
        print("%s/%s done" %(len(rxn_responses),len(all_reactions_ids)))

    if rxn_upper_index != len(all_reactions_ids):
        async_list = []
        last_index = len(all_reactions_ids) - rxn_upper_index
        rxn_lower_index += step
        rxn_upper_index += last_index
        for rxn_id in all_reactions_ids[rxn_lower_index:rxn_upper_index]:
            action_item = grequests.get(url_bigg + "universal/reactions/" +rxn_id)
            async_list.append(action_item) 
        new_responses = [r.json() for r in grequests.map(async_list)]
        rxn_responses += new_responses
    """
    if verbose: print("updating padmet")
    count = 0
    for rxn_id in all_reactions_ids:
        count += 1
        if verbose: print("reaction: %s, %s/%s" %(rxn_id, count, len(all_reactions_ids)))
        rxn_response = requests.get(url_bigg + "universal/reactions/" +rxn_id)
        rxn_dict = rxn_response.json()


        rxn_metabolites = rxn_dict["metabolites"]
        if len(rxn_metabolites) > 1:
            rxn_id = rxn_dict['bigg_id']
            rxn_name = rxn_dict["name"]
        
            all_models_id = [i["bigg_id"] for i in rxn_dict["models_containing_reaction"]]
            async_list = []
            for model_id in all_models_id:
                action_item = grequests.get(url_bigg + "models/"+ model_id +"/reactions/"+ rxn_id)
                async_list.append(action_item)  
            models_responses = [r.json() for r in grequests.map(async_list)]
            all_lower_bound = [i["results"][0]["lower_bound"] for i in models_responses]
            ratio_not_rev = float(all_lower_bound.count(0))/float(len(all_lower_bound))
            if verbose: print("Reaction not reversible in %s/%s model(s)" %(all_lower_bound.count(0), len(all_lower_bound)))
            if ratio_not_rev >= 0.75:
                rxn_direction = "LEFT-TO-RIGHT"
                if verbose: print("Reaction not reversible")
            else:
                rxn_direction = "REVERSIBLE"
                if verbose: print("Reaction reversible")
    
            rxn_node = Node("reaction",rxn_id,{"COMMON_NAME":[rxn_name],"DIRECTION":[rxn_direction]})
            padmetRef.dicOfNode[rxn_id] = rxn_node
        
            rxn_xrefs = rxn_dict["database_links"]
            for db, k in rxn_xrefs.items():
                _id = k[0]["id"]
                meta_max_id += 1
                xref_id = new_id_prefix+"_"+str(meta_max_id)
                xref_node = Node("xref", xref_id, {"DB": [db], "ID": [_id]})
                padmetRef.dicOfNode[xref_id] = xref_node
                has_xref_rlt = Relation(rxn_id, "has_xref", xref_id)
                list_of_relation.append(has_xref_rlt)
        
            
            for metabo_dict in rxn_metabolites:
                metabo_id = metabo_dict["bigg_id"]
                metabo_name = metabo_dict["name"]
                metabo_compart = metabo_dict["compartment_bigg_id"]
                metabo_stoich = metabo_dict["stoichiometry"]
                try:
                    metabo_node = padmetRef.dicOfNode[metabo_id]
                except KeyError:
                    metabo_node = Node("compound",metabo_id,{"COMMON_NAME":[metabo_name]})
                    padmetRef.dicOfNode[metabo_id] = metabo_node
                if metabo_stoich < 0:
                    consumes_rlt = Relation(rxn_id,"consumes",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                    list_of_relation.append(consumes_rlt)
                else:
                    produces_rlt = Relation(rxn_id,"produces",metabo_id,{"STOICHIOMETRY":[abs(metabo_stoich)],"COMPARTMENT":[metabo_compart]})
                    list_of_relation.append(produces_rlt)

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
    
    if verbose: print("Generating file: %s" %output)
    padmetRef.generateFile(output)
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print "done in: ", chrono, "s !"
    
if __name__ == "__main__":
    main()

