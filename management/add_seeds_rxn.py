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
For a given set of compounds representing the growth medium (or seeds). Create 2 reactions
for each compounds to maintain consistency of the network for flux analysis.
For each compounds create 
    -an exchange reaction: this reaction consumes the compound in the
compartment 'C-BOUNDARY' and produces the compound in the compartment 'e' extracellular
    -a transport reaction: this reaction consumes the compound in the compartment
'e' extracellular' and produces the compound in the compartment 'c' cytosol
ex: for seed 'cpd-a'
1/ check if cpd-a in padmetSpec, if not, copy from padmetRef.
2/ create exchange reaction: ExchangeSeed_cpd-a_b: 1 cpd-a (C-BOUNDARAY) <=> 1 cpd-a (e)
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
    boundary_compart = "C-BOUNDARY"
    e_compart = "e"
    c_compart = "c"
    source = "MANUAL"
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
            exchange_rxn_node = Node("reaction", exchange_rxn_id, {"DIRECTION":["REVERSIBLE"]})
            padmetSpec.dicOfNode[exchange_rxn_id] = exchange_rxn_node
            consumption_rlt = Relation(exchange_rxn_id, "consumes", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":[boundary_compart]})
            padmetSpec._addRelation(consumption_rlt)        
            production_rlt = Relation(exchange_rxn_id, "produces", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":[e_compart]})
            padmetSpec._addRelation(production_rlt)

            #reconstructionData:
            reconstructionData_id = exchange_rxn_id+"_reconstructionData_"+source
            if reconstructionData_id in padmetSpec.dicOfNode.keys() and verbose:
                print("Warning: The reaction %s seems to be already added from the same source %s" %(exchange_rxn_id, source))
            reconstructionData = {"SOURCE":[source],"COMMENT":["Added to manage seeds from boundary to extracellular compartment"]}
            reconstructionData_rlt = Relation(exchange_rxn_id,"has_reconstructionData",reconstructionData_id)
            padmetSpec.createNode("reconstructionData", reconstructionData_id, reconstructionData, [reconstructionData_rlt])
            if verbose: print("Creating reconstructionData %s" %reconstructionData_id)                

        transport_rxn_id = "TransportSeed_"+seed_id+"_e"
        if transport_rxn_id not in padmetSpec.dicOfNode.keys():
            if verbose: print("creating trasnport reaction: id = TransportSeed_%s_e" %seed_id)
            transport_rxn_node = Node("reaction", transport_rxn_id, {"DIRECTION":["LEFT-TO-RIGHT"]})
            padmetSpec.dicOfNode[transport_rxn_id] = transport_rxn_node
            consumption_rlt = Relation(transport_rxn_id, "consumes", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":[e_compart]})
            padmetSpec._addRelation(consumption_rlt)        
            production_rlt = Relation(transport_rxn_id, "produces", seed_id, {"STOICHIOMETRY":[1.0],"COMPARTMENT":[c_compart]})
            padmetSpec._addRelation(production_rlt)

            #reconstructionData:
            reconstructionData_id = transport_rxn_id+"_reconstructionData_"+source
            if reconstructionData_id in padmetSpec.dicOfNode.keys() and verbose:
                print("Warning: The reaction %s seems to be already added from the same source %s" %(transport_rxn_id, source))
            reconstructionData = {"SOURCE":[source],"COMMENT":["Added to manage seeds from extracellular to cytosol compartment"]}
            reconstructionData_rlt = Relation(transport_rxn_id,"has_reconstructionData",reconstructionData_id)
            padmetSpec.createNode("reconstructionData", reconstructionData_id, reconstructionData, [reconstructionData_rlt])
            if verbose: print("Creating reconstructionData %s" %reconstructionData_id)                

    padmetSpec.generateFile(output)
        
if __name__ == "__main__":
    main()  

                
    
    