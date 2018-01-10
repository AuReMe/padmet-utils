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
Run flux balance analyse with cobra package. If the flux is >0. Run also FVA
and return result in standard output

usage:
    flux_analysis.py --sbml=FILE
    flux_analysis.py --sbml=FILE --seeds=FILE --targets=FILE
    flux_analysis.py --sbml=FILE --all_species

option:
    -h --help    Show help.
    --sbml=FILE    pathname to the sbml file to test for fba and fva.
    --seeds=FILE    pathname to the sbml file containing the seeds (medium).
    --targets=FILE    pathname to the sbml file containing the targets.
    --all_species    allow to make FBA on all the metabolites of the given model.
"""
from padmet.sbmlPlugin import convert_from_coded_id
from cobra import *
from cobra.io.sbml import create_cobra_model_from_sbml_file
import docopt
import subprocess

def main():
    args = docopt.docopt(__doc__) 
    sbml_file = args["--sbml"]
    seeds_file = args["--seeds"]
    targets_file = args["--targets"]
    if targets_file:
        targets = create_cobra_model_from_sbml_file(targets_file).metabolites
    all_species = args["--all_species"]
    model=create_cobra_model_from_sbml_file(sbml_file)
    try:
        biomassrxn = [rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0][0]
        biomassname = biomassrxn.id
    except IndexError:
        print("Need to set OBJECTIVE COEFFICIENT to '1.0' for the reaction to test")
        exit()
    
    #nb metabolites
    real_metabolites = set([i.id.replace("_"+i.compartment,"") for i in model.metabolites])
    rxn_with_ga = [i for i in model.reactions if i.gene_reaction_rule]
    print("#############")
    print("Model summary")
    print("Number of compounds: %s" %len(real_metabolites))
    print("Number of reactions: %s" %len(model.reactions))
    print("Number of genes: %s" %len(model.genes))
    print("Ratio rxn with genes/rxns: %s%%" %(100*len(rxn_with_ga)/len(model.reactions)))
    
    #
    if seeds_file and targets_file:
        print("#############")
        print("Analyzing targets")
        print("#Topological analysis")
        cmd = "menecheck.py -d %s -s %s -t %s" %(sbml_file, seeds_file, targets_file)
        out = subprocess.check_output(["/bin/bash", "-i", "-c", cmd])
        #prod_targets = (i for i in out.splitlines() if i.endswith(" producible targets:")).next()[:-1]
        print("Number of targets: %s" %(len(targets)))
        print out
        print("#Flux Balance Analysis")
        fba_on_targets(targets, model)
    if all_species:
        targets = model.metabolites
        print("#Flux Balance Analysis on all model metabolites (long process...)")
        fba_on_targets(targets, model)
    print("#############")
    print("Computing optimization")
    solution= model.optimize()
    print("Testing reaction %s" %biomassname)
    print("Growth rate: %s" %solution.objective_value)
    print("Status: %s" %solution.status)
    model.summary()
    if (solution.objective_value > 1e-5):
        blocked = flux_analysis.find_blocked_reactions(model, model.reactions)
        essRxns = flux_analysis.find_essential_reactions(model)
        essGenes = flux_analysis.find_essential_genes(model)

        print('FVA analysis:')
        print('\tBlocked reactions: %s' %len(blocked))
        print('\tEssential reactions: %s' %len(essRxns))
        print('\tEssential genes: %s' %len(essGenes))
    
    #get biomass rxn reactants
    bms_reactants = dict([(k,v) for k,v in biomassrxn.metabolites.items() if v < 0])
    bms_products = dict([(k,v) for k,v in biomassrxn.metabolites.items() if v > 0])
    dict_output = {"positive":{},"negative":{}}
    #for each metabolite in reactant, create a biomass rxn with only this metabolite in reactants
    biomassrxn.objective_coefficient = 0.0
    for reactant, stoich in bms_reactants.items():
        test_rxn = Reaction("test_rxn")
        test_rxn.lower_bound = 0
        test_rxn.upper_bound = 1000
        metabolitedict = dict(bms_products)
        metabolitedict.update({reactant:stoich})

        model.add_reaction(test_rxn)
        test_rxn.add_metabolites(metabolitedict)
        test_rxn.objective_coefficient = 1.0

        solution = model.optimize()
        if (solution.objective_value > 1e-5):
            dict_output["positive"][reactant] = solution.objective_value
        else:
            dict_output["negative"][reactant] = solution.objective_value
        model.remove_reactions([test_rxn])
    print("%s/%s compounds with positive flux" %(len(dict_output["positive"].keys()), len(bms_reactants)))
    print("%s/%s compounds with negative flux" %(len(dict_output["negative"].keys()), len(bms_reactants)))

    for k,v in dict_output["positive"].items():
        print("%s // %s %s positive" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))
    for k,v in dict_output["negative"].items():
        print("%s // %s %s negative" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))


def fba_on_targets(allspecies, model):
    dict_output = {"positive":{},"negative":{}}
    for species in allspecies:
        #lets create a copy of the initial model
        model2 = model.copy()
        #Create a new reaction that consume the given species
        FBA_rxn = Reaction("FBA_TEST")
        FBA_rxn.lower_bound = 0
        FBA_rxn.upper_bound = 1000
        model2.add_reaction(FBA_rxn)
        FBA_rxn.objective_coefficient = 1.0
        metabolitedict = {}
        metabolitedict[species]=-1.0
        FBA_rxn.add_metabolites(metabolitedict)
                
        solution = model2.optimize()
        if (solution.objective_value > 1e-5):
            dict_output["positive"][species] = solution.objective_value
        else:
            dict_output["negative"][species] = solution.objective_value
    print("%s/%s compounds with positive flux" %(len(dict_output["positive"].keys()), len(allspecies)))
    print("%s/%s compounds with negative flux" %(len(dict_output["negative"].keys()), len(allspecies)))

    for k,v in dict_output["positive"].items():
        print("%s // %s %s positive" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))
    for k,v in dict_output["negative"].items():
        print("%s // %s %s negative" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))
    

if __name__ == "__main__":
    main()