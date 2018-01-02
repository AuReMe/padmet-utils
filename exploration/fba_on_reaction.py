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
    fba_on_reaction.py --sbml=FILE

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
    solution = model.optimize()
    try:
        biomassrxn = [rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0][0]
        biomassname = biomassrxn.id
    except IndexError:
        print("Need to set OBJECTIVE COEFFICIENT to '1.0' for the reaction to test")
        exit()
    print("Testing reaction %s" %biomassname)
    print( 'growth rate: %s\nStatus: %s.' %(solution.objective_value, solution.status))
    if (solution.objective_value > 1e-5):
        #FVA_result = flux_analysis.variability.flux_variability_analysis(model, fraction_of_optimum=1.0)
        FVA_result = flux_analysis.variability.flux_variability_analysis(model, fraction_of_optimum=1.0, tolerance_feasibility=1e-6, tolerance_integer=1e-6)
        #print(str(FVA_result)+"\n")    

        essential, alternative, blocked = [[] for i in range(3)]
        for k,v in FVA_result.iterrows():
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




if __name__ == "__main__":
    main()