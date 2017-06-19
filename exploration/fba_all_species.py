#!/usr/bin/python
#-*- coding: utf-8 -*-
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

@author: Clemence FRIOUX, clemence.frioux@inria.fr
Description:
For each compound in ListOfSpecies,the biomass reaction is (re)created with
this specie as a reactant (stoichio = -1), and the FBA calculated accordingly

usage:
    fba_all_species.py --sbml=FILE

option:
    -h --help    Show help.
    --sbml=FILE    metabolic network in SBML format.
"""
import docopt
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra import *

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    inputFile = args["--sbml"]
    model=create_cobra_model_from_sbml_file(inputFile)
    #print(model.reactions)
    try:
        biomassrxn = [rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0][0]
        biomassrxn.objective_coefficient = 0.0
    except IndexError:
        pass

    allspecies = model.metabolites

    dict_output = {"positive":{},"negative":{}}
    for species in allspecies:
        #lets create a copy of the initial model
        model2 = model.copy()
        #Create a new reaction that consume the given species
        FBA_rxn = Reaction("FBA_TEST")
        FBA_rxn.lower_bound = 0
        FBA_rxn.upper_bound = 1000
        FBA_rxn.objective_coefficient = 1.0
        metabolitedict = {}
        metabolitedict[species]=-1.0
        FBA_rxn.add_metabolites(metabolitedict)

        model2.add_reaction(FBA_rxn)
        model2.optimize()
        if (model2.solution.f > 1e-5):
            dict_output["positive"][species] = model2.solution.f
        else:
            dict_output["negative"][species] = model2.solution.f

    for k,v in dict_output["positive"].items():
        print("%s %s positive" %(k,v))
    for k,v in dict_output["negative"].items():
        print("%s %s negative" %(k,v))
    print("%s/%s compounds with positive flux" %(len(dict_output["positive"].keys()), len(allspecies)))
    print("%s/%s compounds with negative flux" %(len(dict_output["negative"].keys()), len(allspecies)))

