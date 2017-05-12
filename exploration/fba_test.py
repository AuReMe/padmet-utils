# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Run flux balance analyse with cobra package. If the flux is >0. Run also FVA
and return result in standard output

usage:
    fba_test.py --sbml=FILE

option:
    -h --help    Show help.
    --sbml=FILE    pathname to the sbml file to test for fba and fva.
"""
from padmet.sbmlPlugin import convert_from_coded_id
from cobra import *
from cobra.io.sbml import create_cobra_model_from_sbml_file
import docopt

def main():
    args = docopt.docopt(__doc__) 
    sbml_file = args["--sbml"]
    model=create_cobra_model_from_sbml_file(sbml_file)
    model.optimize()
    print( 'growth rate: '+str(model.solution.f)+ '\n'+ 'status: '+ model.solution.status)
    if (model.solution.f > 1e-5):
        #FVA_result = flux_analysis.variability.flux_variability_analysis(model, fraction_of_optimum=1.0)
        FVA_result = flux_analysis.variability.flux_variability_analysis(model, fraction_of_optimum=1.0, tolerance_feasibility=1e-6, tolerance_integer=1e-6)
        #print(str(FVA_result)+"\n")    

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
        biomassrxn = [rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0][0]
        biomassname = biomassrxn.id
        metabolites = biomassrxn._metabolites
        # print(metabolites)
        # only append to the list compounds that are reactants
        bms_metabolites_list = [x.id for x in metabolites]
        bms_reactants = [x.id for x in metabolites if metabolites[x] < 0]
        # for metabo in metabolites:
        #     metabolites_list.append(metabo.id)
    
        #print(metabolites_list)
        dict_output = {"positive":{},"negative":{}}
        has_been_tested = []
        for metabolite1 in bms_reactants:
            if not metabolite1 in has_been_tested:
                #lets create a copy of the initial model
                model2 = model.copy()
                #get the biomass reaction
                biomassrxn2 = model2.reactions.get_by_id(biomassname) #R_Ec_biomass_iAF1260_core_59p81M
                # empty the biomass metabolites
                for m in bms_metabolites_list:
                    biomassrxn2.pop(m)
                #create a clean dict with only compound of interest and stoichio to -1 (reactant)
                metabolitedict = {}
                metabolitedict[metabolite1]=-1.0
                # replace the biomass metabolites with only metabolite1 inside
                biomassrxn2.add_metabolites(metabolitedict)
                # print(biomassrxn2._metabolites)
                #test flux balance analysis
                #print(biomassrxn2._metabolites)
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
        print("%s/%s compounds with negative flux" %(len(dict_output["negative"].keys()), len(has_been_tested)))




if __name__ == "__main__":
    main()