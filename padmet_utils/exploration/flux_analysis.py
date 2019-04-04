    # -*- coding: utf-8 -*-
"""
Description:
    Run flux balance analyse with cobra package. If the flux is >0. Run also FVA
    and return result in standard output

::
    
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
from padmet.utils.sbmlPlugin import convert_from_coded_id
from cobra import Reaction
from cobra import flux_analysis as cobra_flux_analysis
from cobra.io.sbml import create_cobra_model_from_sbml_file
import docopt
import subprocess

def main():
    args = docopt.docopt(__doc__) 
    sbml_file = args["--sbml"]
    seeds_file = args["--seeds"]
    targets_file = args["--targets"]
    all_species = args["--all_species"]
    flux_analysis(sbml_file, seeds_file, targets_file, all_species)

def flux_analysis(sbml_file, seeds_file = None, targets_file = None, all_species = False):
    """
    Run flux balance analyse with cobra package. If the flux is >0. Run also FVA
    and return result in standard output

    Parameters
    ----------
    sbml_file: str
        path to sbml file to analyse
    seeds_file: str
        path to sbml file with only compounds representing the seeds/growth medium
    targets_file: str
        path to sbml file with only compounds representing the targets to reach
    all_species: bool
        if True will try to create obj function for each compound and return which are reachable by flux.
        
    """
    if targets_file:
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
    print("Ratio rxn with genes/rxns: %s%%" %(100*len(rxn_with_ga)/len(model.reactions)))
    
    #
    if seeds_file and targets_file:
        try:
            print("#############")
            print("Analyzing targets")
            print("#Topological analysis")
            cmd = "menecheck.py -d %s -s %s -t %s" %(sbml_file, seeds_file, targets_file)
            out = subprocess.check_output(["/bin/bash", "-i", "-c", cmd])
            #prod_targets = (i for i in out.splitlines() if i.endswith(" producible targets:")).next()[:-1]
            print("Number of targets: %s" %(len(targets)))
            print(out.decode("UTF-8"))
        except subprocess.CalledProcessError:
            print("Menetools are not install, Can't run topological analysis")
        print("#Flux Balance Analysis")
        fba_on_targets(targets, model)
    if all_species:
        targets = model.metabolites
        print("#Flux Balance Analysis on all model metabolites (long process...)")
        fba_on_targets(targets, model)
        return
    print("#############")
    print("Computing optimization")
    solution= model.optimize()
    print("Testing reaction %s" %biomassname)
    print("Growth rate: %s" %solution.objective_value)
    print("Status: %s" %solution.status)
    model.summary()
    if (solution.objective_value > 1e-5):
        blocked = cobra_flux_analysis.find_blocked_reactions(model, model.reactions)
        essRxns = cobra_flux_analysis.find_essential_reactions(model)
        essGenes = cobra_flux_analysis.find_essential_genes(model)

        print('FVA analysis:')
        print('\tBlocked reactions: %s' %len(blocked))
        print('\tEssential reactions: %s' %len(essRxns))
        print('\tEssential genes: %s' %len(essGenes))
    
    #get biomass rxn reactants
    bms_reactants = dict([(k,v) for k,v in list(biomassrxn.metabolites.items()) if v < 0])
    bms_products = dict([(k,v) for k,v in list(biomassrxn.metabolites.items()) if v > 0])
    dict_output = {"positive":{},"negative":{}}
    #for each metabolite in reactant, create a biomass rxn with only this metabolite in reactants
    biomassrxn.objective_coefficient = 0.0
    for reactant, stoich in list(bms_reactants.items()):
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
    print("%s/%s compounds with positive flux" %(len(list(dict_output["positive"].keys())), len(bms_reactants)))
    print("%s/%s compounds without flux" %(len(list(dict_output["negative"].keys())), len(bms_reactants)))

    for k,v in list(dict_output["positive"].items()):
        print("%s // %s %s positive" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))
    for k,v in list(dict_output["negative"].items()):
        print("%s // %s %s NULL" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))


def fba_on_targets(allspecies, model):
    """
    for each specie in allspecies, create an objective function with the current species as only product
    and try to optimze the model and get flux.

    Parameters
    ----------
    allSpecies: list
        list of species ids to test
    model: cobra.model
        Cobra model from a sbml file
    """
    #dict_output = {"positive":{},"negative":{}}
    for species in allspecies:
        #lets create a copy of the initial model
        model2 = model.copy()
        #remove all obj coef
        for rxn in model2.reactions:
            if rxn.objective_coefficient == 1.0:
                rxn.objective_coefficient = 0.0
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
            print("%s // %s %s positive" %(species, convert_from_coded_id(species.id)[0]+"_"+convert_from_coded_id(species.id)[2], solution.objective_value))
        else:
            print("%s // %s %s NULL" %(species, convert_from_coded_id(species.id)[0]+"_"+convert_from_coded_id(species.id)[2], solution.objective_value))

if __name__ == "__main__":
    main()
