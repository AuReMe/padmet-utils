#!/usr/bin/python
#-*- coding: utf-8 -*-
"""
@author: Clemence FRIOUX, clemence.frioux@inria.fr
Description:
For a given compound, the biomass reaction is recreated with 
this compound only as a reactant (stoichio = -1), 
and the FBA calculated accordingly if FBA is positive, FVA is calculated

usage:
    change_biomass_cobra.py --network=FILE --compound=ID

option:
    -h --help    Show help.
    --network=FILE    metabolic network in SBML format.
    --compound=ID    ID of the compound of interest for which FBA must be tested
"""
import docopt
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra import *
from padmet.sbmlPlugin import convert_from_coded_id

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    inputFile = args["--network"]
    compoundofinterest = args["--compound"]
    model=create_cobra_model_from_sbml_file(inputFile)
    #print(model.reactions)
    try:
        biomassrxn = [rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0][0]
        biomassname = biomassrxn.id
    except KeyError:
        print("biomass reaction must be at OBJECTIVE_COEFFICIENT=1")
        exit()
    #print(model.reactions)

    #metabolites = biomassrxn._metabolites
    bms_metabolites = biomassrxn._metabolites
    bms_metabolites_list = [metabo.id for metabo in bms_metabolites]

    metabolites_list = [compoundofinterest]
    #for metabo in metabolites:
    #    metabolites_list.append(metabo.id)

    #print(metabolites_list)

    has_been_tested = []
    for metabolite1 in metabolites_list:
        if not metabolite1 in has_been_tested:
            #lets create a copy of the initial model
            model2 = model.copy()
            #get the biomass reaction
            biomassrxn2 = model2.reactions.get_by_id(biomassname) #R_Ec_biomass_iAF1260_core_59p81M
            #get the metabolites of the biomass
            metabolites2_list = []
            # empty the biomass metabolites
            for m in bms_metabolites_list:
                biomassrxn2.pop(m)
            #create a clean dict with only compound of interest and stoichio to -1 (reactant)
            metabolitedict = {}
            metabolitedict[metabolite1]=-1.0
            # replace the biomass metabolites with only metabolite1 inside
            biomassrxn2.add_metabolites(metabolitedict)
            #test flux balance analysis
            #print(biomassrxn2._metabolites)
            model2.optimize()
            if (model2.solution.f > 1e-5):
                print(metabolite1+' '+str(model2.solution.f)+' '+'positive')
                FVA_result = flux_analysis.variability.flux_variability_analysis(model2, fraction_of_optimum=1.0, tolerance_feasibility=1e-6, tolerance_integer=1e-6)
                essential, alternative, blocked = [[] for i in range(3)]
                for k,v in FVA_result.iteritems():
                    _min = float(v['minimum'])
                    _max = float(v['maximum'])
                    if (_min > 0.0001 and _max > 0.0001) or (_min < -0.0001 and _max < -0.0001):
                        essential.append(k)
                    elif _min > -0.0001 and _min < 0.0001 and _max > -0.0001 and _max < 0.0001:
                        blocked.append(k)
                    else:
                        alternative.append(k)
        
                print('FVA analysis:')
                print('\t'+'Essential reactions: '+str(len(essential)))
                print('\t'+'Alternative reactions: '+str(len(alternative)))
                print('\t'+'Blocked reactions: '+str(len(blocked)))
                
                print('Essential reactions: '+str(len(essential)))
                for r_id in essential:
                    r_id_decoded = convert_from_coded_id(r_id)[0]
                    print('\t'+r_id_decoded+' ('+r_id+')')
                print('Alternative reactions: '+str(len(alternative)))
                for r_id in alternative:
                    r_id_decoded = convert_from_coded_id(r_id)[0]
                    print('\t'+r_id_decoded+' ('+r_id+')')
                print('Blocked reactions: '+str(len(blocked)))
                for r_id in blocked:
                    r_id_decoded = convert_from_coded_id(r_id)[0]
                    print('\t'+r_id_decoded+' ('+r_id+')')
            else:
                print(metabolite1+' '+str(model2.solution.f)+' '+'negative')
            has_been_tested.append(metabolite1)
        else:
            pass
