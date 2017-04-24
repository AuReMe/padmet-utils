#!/usr/bin/python
#-*- coding: utf-8 -*-
"""
@author: Clemence FRIOUX, clemence.frioux@inria.fr
Description:
For each compound in ListOfSpecies,the biomass reaction is (re)created with
this specie as a reactant (stoichio = -1), and the FBA calculated accordingly

usage:
    change_biomass_cobra.py --network=FILE

option:
    -h --help    Show help.
    --network=FILE    metabolic network in SBML format.
"""
import docopt
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra import *

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    inputFile = args["--network"]
    model=create_cobra_model_from_sbml_file(inputFile)
    #print(model.reactions)
    try:
        biomassrxn = [rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0][0]
        biomassname = biomassrxn.id
    except KeyError:
        print("biomass reaction must be at OBJECTIVE_COEFFICIENT=1")
        exit()

    allspecies = model.metabolites
    # species_list = []
    #print(len(allspecies))
    species_list = [metabo.id for metabo in allspecies]
    # for metabo in allspecies:
    #     species_list.append(metabo.id)
    #print(species_list)
    bmsmetabolites = biomassrxn._metabolites
    bmsmetabolites_list = [metabo.id for metabo in bmsmetabolites]
    # for metabo in bmsmetabolites:
    #     bmsmetabolites_list.append(metabo.id)

    #print(metabolites_list)
    dict_output = {"positive":{},"negative":{}}
    has_been_tested = []
    for metabolite1 in allspecies:
        if not metabolite1 in has_been_tested:
            #lets create a copy of the initial model
            model2 = model.copy()
            #get the biomass reaction
            biomassrxn2 = model2.reactions.get_by_id(biomassname) #R_Ec_biomass_iAF1260_core_59p81M
            #get the metabolites of the biomass
            # empty the biomass metabolites
            for m in bmsmetabolites_list:
                biomassrxn2.pop(m)
            #create a clean dict with only compound of interest and stoichio to -1 (reactant)
            metabolitedict = {}
            metabolitedict[metabolite1]=-1.0
            #create a biomass with only metabolite1 inside
            biomassrxn2.add_metabolites(metabolitedict)
            #test flux balance analysis
            #print(biomassrxn2._metabolites)
            #print(model2.solution.f)
            model2.optimize()
            if (model2.solution.f > 1e-5):
                dict_output["positive"][metabolite1] = model2.solution.f
            else:
                dict_output["negative"][metabolite1] = model2.solution.f
            has_been_tested.append(metabolite1)
        else:
            pass
    for k,v in dict_output["positive"].items():
        print("%s %s positive" %(k,v))
    for k,v in dict_output["negative"].items():
        print("%s %s negative" %(k,v))
    print("%s/%s compounds with positive flux" %(len(dict_output["positive"].keys()), len(has_been_tested)))
    print("%s/%s compounds with positive flux" %(len(dict_output["negative"].keys()), len(has_been_tested)))

