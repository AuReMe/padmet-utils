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
if __name__ == "__main__":
    main()