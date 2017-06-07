# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
For a given set of compounds representing the growth medium (or seeds). Create 2 reactions
for each compounds to maintain consistency of the network for flux analysis.
For each compounds create 
    -an exchange reaction: this reaction consumes the compound in the
compartment 'C-BOUNDARY' and produces the compound in the compartment 'e' extracellular
    -a transport reaction: this reaction consumes the compound in the compartment
'e' extracellular' and produces the compound in the compartment 'c' cytosol
ex: for seed 'cpd-a'
1/ check if cpd-a in padmetSpec, if not, copy from padmetRef.
2/ create exchange reaction: ExchangeSeed_cpd-a_b: 1 cpd-a (C-BOUNDARAY) => 1 cpd-a (e)
3/ create transport reaction: TransportSeed_cpd-a_e: 1 cpd-a (e) => 1 cpd-a (c)
4/ create a new file if output not None, or overwritte padmetSpec

usage:
    add_seeds_rxn.py --padmetSpec=FILE --padmetRef=FILE --seeds=FILE [--output=FILE] [-v]

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathanme to the padmet file to update
    --padmetRef=FILE    pathname to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
    --seeds=FILE    the pathname to the file containing the compounds ids and the compart, line = cpd-id\tcompart.
    --output=FILE    If not None, pathname to the padmet file updated
    -v   print info

"""
from padmet.padmetSpec import PadmetSpec
from padmet.padmetRef import PadmetRef
from padmet.node import Node
from padmet.relation import Relation
import docopt

def main():
    args = docopt.docopt(__doc__)
    with open(args["--seeds"], 'r') as f:
        seeds = [line.split("\t")[0] for line in f.read().splitlines()]
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    padmetRef  = PadmetRef(args["--padmetRef"])
    output = args["--output"]
    verbose = args["-v"]
    if output is None:
        output = args["--padmetSpec"]
    
    for seed_id in seeds:
        #check if seed in padmetSpec or in padmetRef
        try:
            padmetSpec.dicOfNode[seed_id]
        except KeyError:
            if verbose: print("%s not in the network" %seed_id)
            try:
                if verbose: print("Try to copy from dbref")
                padmetSpec._copyNodeExtend(padmetRef, seed_id)
            except KeyError:
                if verbose:
                    print("%s not in the padmetRef" %seed_id)
                    print("creating a new compound")
                #not in padmetRef and padmetSpec, create compound and transport/exchange rxn
                seed_node = Node("compound", seed_id)
                padmetSpec.dicOfNode[seed_id] = seed_node
                if verbose:
                    print("new compound created: id = %s" %seed_id)
        exchange_rxn_id = "ExchangeSeed_"+seed_id+"_b"
        if exchange_rxn_id not in padmetSpec.dicOfNode.keys():
            if verbose: print("creating exchange reaction: id = ExchangeSeed_%s_b" %seed_id)
            exchange_rxn_node = Node("reaction", exchange_rxn_id, {"DIRECTION":["REVERSIBLE"],"SOURCE":["manual"],"COMMENT":["Added to manage seeds from boundary to extracellular compartment"]})
            padmetSpec.dicOfNode[exchange_rxn_id] = exchange_rxn_node
            consumption_rlt = Relation(exchange_rxn_id, "consumes", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":["C-BOUNDARY"]})
            padmetSpec._addRelation(consumption_rlt)        
            production_rlt = Relation(exchange_rxn_id, "produces", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":["e"]})
            padmetSpec._addRelation(production_rlt)        

        transport_rxn_id = "TransportSeed_"+seed_id+"_e"
        if transport_rxn_id not in padmetSpec.dicOfNode.keys():
            if verbose: print("creating trasnport reaction: id = TransportSeed_%s_e" %seed_id)
            transport_rxn_node = Node("reaction", transport_rxn_id, {"DIRECTION":["REVERSIBLE"],"SOURCE":["manual"],"COMMENT":["Added to manage seeds from extracellular to cytosol compartment"]})
            padmetSpec.dicOfNode[transport_rxn_id] = transport_rxn_node
            consumption_rlt = Relation(transport_rxn_id, "consumes", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":["e"]})
            padmetSpec._addRelation(consumption_rlt)        
            production_rlt = Relation(transport_rxn_id, "produces", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":["c"]})
            padmetSpec._addRelation(production_rlt)
    padmetSpec.generateFile(output)
        
if __name__ == "__main__":
    main()  

                
    
    