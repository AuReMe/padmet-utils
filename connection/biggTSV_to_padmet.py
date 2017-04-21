# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Create a padmet file from 2 tsv files (reactions.tsv and metabolites.tsv)
metabolites.tsv: col1 = metabolite ID, col2 = metabolite Name 
ex: 10fthf_c 10-Formyltetrahydrofolate
reactions.tsv: col1 = reaction ID, col2 = reaction Name, col3 = formula, col4 = reversibility
ex: 10FTHF5GLUtl 10FTHF5GLUtl 1.0:10fthf5glu_c => 1.0:10fthf5glu_l False
NB: the formula must be in this format: STOICHIOMETRY:METABOLITE_ID + ... => STOICHIOMETRY:METABOLITE_ID

for reaction in reactions:
    /!\ reaction id is modified with sbmlPlugin.converto_from_coded_id. delete the type information (R_)
    Do not copy reactions without reactant or without products
    create Node, class = reaction, id = id, misc.COMMON_NAME = name, misc.DIRECTION = reversible ('True' = 'REVERSIBLE', 'False' = 'LEFT-TO-RIGHT')
    Parse formula to recovere listOfReactants, listOfProducts
    for metabolite in listOfReactants:
        extract compart with convert_from_coded_id. reactant_id = id without compart (if compart not None)
        create relation, type = consumes, idIn = reaction_id, idOut = reactant_id, relation.misc["COMPARTMENT"] = [compart], relation.misc["STOICHIOMETRY"] = [stoic] 
    for metabolite in listOfProducts:
        extract compart with convert_from_coded_id. product_id = id without compart (if compart not None)
        create relation, type = produces, idIn = reaction_id, idOut = product_id, relation.misc["COMPARTMENT"] = [compart/Unknown], relation.misc["STOICHIOMETRY"] = [stoic] 

for metabolite in metabolites.tsv:
    /!\ metabolite id is modified. => delete the type information (M_) and the compart info (_c/e/p..)
    extract compart with convert_from_coded_id. product_id = id without compart (if compart not None)
    add in node.misc["COMMON_NAME"] = metaboName

usage:
    biggTSV_to_padmet.py --version=STR --output=FILE --rxn_file=FILE [--met_file=FILE] [-v]

options:
    -h --help     Show help.
    --version=STR    The tag associated to the origin of the information used
    --output=FILE    pathname of the padmet file to create
    --rxn_file=FILE    pathanme of the reactions file. col1 = reaction ID, col2 = reaction Name, col3 = formula, col4 = reversibility.
    --met_file=FILE    pathanme of the metabolites file. col1 = metabolite ID, col2 = metabolite Name 
    -v   print info

"""
from padmet.node import Node
from padmet.relation import Relation
from padmet.padmetRef import PadmetRef
from padmet.sbmlPlugin import convert_from_coded_id,convert_to_coded_id
from datetime import datetime
from time import time
import re
import docopt

def main():
    chronoDepart = time()
    #parsing args
    args = docopt.docopt(__doc__)
    version = args["--version"]
    output = args["--output"]
    rxn_file = args["--rxn_file"]
    met_file = args["--met_file"]
    verbose = args["-v"]
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")
    #print(verbose,today_date,version, output, classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file)
    #return    
    policyInArray = [['compound','has_xref','xref'], ['compound','has_suppData','suppData'], 
                    ['gene','has_xref','xref'], ['gene','has_suppData','suppData'],
                    ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','is_linked_to','gene'],
                    ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y']] 
    dbNotes = {"PADMET":{"Creation":today_date, "version":"2.3"}, "DB_ref":{"db":"BIGG", "version":version}}
    padmetRef = PadmetRef()
    if verbose: print("setting policy")
    padmetRef.setPolicy(policyInArray)
    if verbose: print("setting dbInfo")
    padmetRef.setInfo(dbNotes)

    if verbose: print ("loading rxn file: " + rxn_file)
    with open(rxn_file, 'r') as f:
        rxn_fileInArray = [line.split("\t") for line in f.read().splitlines()[1:]]
    if verbose: print(str(len(rxn_fileInArray))+" reactions")
    for rxn_id, rxn_name, rxn_formula, rxn_rev in rxn_fileInArray:
        reactant_part, product_part = rxn_formula.split(" => ")
        if reactant_part and product_part:
            reactants = [(r.split(":")) for r in reactant_part.split("+")]
            products = [(r.split(":")) for r in product_part.split("+")]
            #if reactants or prodcuts != 0        
            if rxn_rev == "True":
                reaction_dir = "REVERSIBLE"
            else:
                reaction_dir = "LEFT-TO-RIGHT"
            rxn_node = Node("reaction", rxn_id)
            rxn_node.misc.update({"COMMON_NAME":[rxn_name], "DIRECTION":[reaction_dir]})
            padmetRef._addNode(rxn_node)
            for reactant in reactants:
                reactant_stoich = reactant[0]
                reactant_id, _type, compart = convert_from_coded_id(reactant[1].strip())
                if compart is None:
                    compart = "unknown"
                if reactant_id not in padmetRef.dicOfNode.keys():
                    reactant_node = Node("compound",reactant_id)
                    padmetRef.dicOfNode[reactant_id] = reactant_node
                consumes_rlt = Relation(rxn_id,"consumes",reactant_id)
                consumes_rlt.misc.update({"STOICHIOMETRY":[reactant_stoich],"COMPARTMENT":[compart]})
                padmetRef._addRelation(consumes_rlt)
            try:
                for product in products:
                    product_stoich = product[0]
                    product_id, _type, compart = convert_from_coded_id(product[1].strip())
                if compart is None:
                    compart = "unknown"
                if product_id not in padmetRef.dicOfNode.keys():
                    product_node = Node("compound",product_id)
                    padmetRef.dicOfNode[product_id] = product_node
                produces_rlt = Relation(rxn_id,"produces",product_id)
                produces_rlt.misc.update({"STOICHIOMETRY":[product_stoich],"COMPARTMENT":[compart]})
                padmetRef._addRelation(produces_rlt)
            except IndexError:
                pass
    
    if met_file is not None:
        if verbose: print ("loading metabolites file: " + met_file)
        #a file with a header
        with open(met_file, 'r') as f:
            met_fileInArray = [line.split("\t") for line in f.read().splitlines()[1:]]
        if verbose: print("%s compounds" %len(met_fileInArray))
        if verbose: print("updating compounds node")
        for cpd_id, cpd_name in met_fileInArray:
            #ex: M_x_e => x,M,e , x => x,None,None
            uncoded_cpd_id, _type, cpd_compart = convert_from_coded_id(cpd_id)
            try:
                padmetRef.dicOfNode[uncoded_cpd_id].misc["COMMON_NAME"] = [cpd_name]
            except KeyError:
                pass

                
    if verbose: print("Generating file: %s" %output)
    padmetRef.generateFile(output)
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print "done in: ", chrono, "s !" 

if __name__ == "__main__":
    main()