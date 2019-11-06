    # -*- coding: utf-8 -*-
"""
Description:
    1./ Run flux balance analyse with cobra package on an already defined reaction.
    Need to set in the sbml the value 'objective_coefficient' to 1.
    If the reaction is reachable by flux: return the flux value and the flux value
    for each reactant of the reaction.
    If not: only return the flux value for each reactant of the reaction.
    If a reactant has a flux of '0' this means that it is not reachable by flux
    (and maybe topologically). To unblock the reaction it is required to fix the
    metabolic network by adding/removing reactions until all reactant are reachable.
    
    2./If seeds and targets given as sbml files with only compounds.
    Will also try to use the Menetools library to make a topologicall analysis.
    Topological reachabylity of the targets compounds from the seeds compounds.
    
    3./ If --all_species: will test flux reachability of all the compounds in the
    metabolic network (may take several minutes)

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
from padmet.utils.exploration import flux_analysis
import docopt

def main():
    args = docopt.docopt(__doc__) 
    sbml_file = args["--sbml"]
    seeds_file = args["--seeds"]
    targets_file = args["--targets"]
    all_species = args["--all_species"]
    flux_analysis.flux_analysis(sbml_file, seeds_file, targets_file, all_species)

if __name__ == "__main__":
    main()
