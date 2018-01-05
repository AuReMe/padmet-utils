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
    flux_analysis.py --sbml=FILE [--full_report] --seeds=FILE --targets=FILE
    flux_analysis.py --sbml=FILE --all_species
    flux_analysis.py --sbml=FILE --species=FILE

option:
    -h --help    Show help.
    --sbml=FILE    pathname to the sbml file to test for fba and fva.
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

    targets = create_cobra_model_from_sbml_file(targets_file).metabolites
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
    print("ratio rxn_ga/rxns: %s%%" %(100*len(rxn_with_ga)/len(model.reactions)))
    cmd = "menecheck.py -d %s -s %s -t %s" %(sbml_file, seeds_file, targets_file)
    out = subprocess.check_output(["/bin/bash", "-i", "-c", cmd])
    prod_targets = (i for i in out.splitlines() if i.endswith(" producible targets:")).next()[:-1]
    print("%s targets" %(len(targets)))
    print(prod_targets)
    fba_given_species(targets, model)
    solution= model.optimize()
    print("Testing reaction %s" %biomassname)
    print( 'growth rate: %s\nStatus: %s.' %(solution.objective_value, solution.status))
    if (solution.objective_value > 1e-5):
        blocked = flux_analysis.find_blocked_reactions(model, model.reactions)
        essRxns = flux_analysis.find_essential_reactions(model)
        essGenes = flux_analysis.find_essential_genes(model)

        print('FVA analysis:')
        print('\tBlocked reactions: %s' %len(blocked))
        print('\tEssential reactions: %s' %len(essRxns))
        print('\tEssential genes: %s' %len(essGenes))
    
    else:
        # print(metabolites)
        # only append to the list compounds that are reactants
        bms_reactants = biomassrxn.reactants
        biomassrxn.objective_coefficient = 0
        # for metabo in metabolites:
        #     metabolites_list.append(metabo.id)
    
        #print(metabolites_list)
        dict_output = {"positive":{},"negative":{}}
        for metabolite in bms_reactants:
            #lets create a copy of the initial model
            model2 = model.copy()
            #get the biomass reaction
            biomassrxn2 = model2.reactions.get_by_id(biomassname) #R_Ec_biomass_iAF1260_core_59p81M
            biomassrxn2.objective_coefficient = 1
            # empty the biomass metabolites
            for i in range(len(bms_reactants)):
                biomassrxn2.reactants.pop(i)
            #create a clean dict with only compound of interest and stoichio to -1 (reactant)
            metabolitedict = {}
            metabolitedict[metabolite]=-1.0
            # replace the biomass metabolites with only metabolite1 inside
            biomassrxn2.add_metabolites(metabolitedict)
            # print(biomassrxn2._metabolites)
            #test flux balance analysis
            #print(biomassrxn2._metabolites)
            solution2 = model2.optimize()
            if (solution2.objective_value > 1e-5):
                dict_output["positive"][metabolite] = solution2.objective_value
            else:
                dict_output["negative"][metabolite] = solution2.objective_value

        for k,v in dict_output["positive"].items():
            print("%s %s positive" %(k,v))
        for k,v in dict_output["negative"].items():
            print("%s %s negative" %(k,v))
        print("%s/%s compounds with positive flux" %(len(dict_output["positive"].keys()), len(bms_reactants)))
        print("%s/%s compounds with negative flux" %(len(dict_output["negative"].keys()), len(bms_reactants)))


def fba_given_species(allspecies, model):
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
    """
    for k,v in dict_output["positive"].items():
        print("%s %s positive" %(k,v))
    for k,v in dict_output["negative"].items():
        print("%s %s negative" %(k,v))
    """
    print("%s/%s compounds with positive flux" %(len(dict_output["positive"].keys()), len(allspecies)))
    print("%s/%s compounds with negative flux" %(len(dict_output["negative"].keys()), len(allspecies)))
    

if __name__ == "__main__":
    main()