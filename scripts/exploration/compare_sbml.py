# -*- coding: utf-8 -*-
"""
Description:
    compare reactions in two sbml

::

    usage:
        compare_sbml.py --sbml1=FILE --sbml2=FILE
    
    option:
        -h --help    Show help.
        --sbml1=FILE    path of the first sbml file
        --sbml2=FILE    path of the second sbml file
"""
from cobra import *
from cobra.io.sbml import create_cobra_model_from_sbml_file
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml_1 = create_cobra_model_from_sbml_file(args["--sbml1"])
    sbml_2 = create_cobra_model_from_sbml_file(args["--sbml2"])
    
    print("sbml1:")
    print("metabolites: %s" %(len(sbml_1.metabolites)))
    print("reactions: %s" %(len(sbml_1.reactions)))
    print("sbml2:")
    print("metabolites: %s" %(len(sbml_2.metabolites)))
    print("reactions: %s" %(len(sbml_2.reactions)))
    not_in1 = [i for i in sbml_2.reactions if i not in sbml_1.reactions]
    print("reactions not in sbml1: %s" %len(not_in1))
    for i in not_in1:
        print("\t%s" %i)
    not_in2 = [i for i in sbml_1.reactions if i not in sbml_2.reactions]
    print("reactions not in sbml2: %s" %len(not_in2))
    for j in not_in2:
        print("\t%s" %j)
    
    all_diff = set()
    for rxn1 in sbml_1.reactions:
        rxn_id = rxn1.id
        try:
            rxn2 = sbml_2.reactions.get_by_id(rxn_id)
            if not compare_rxn(rxn1, rxn2):
                print("%s is diff" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass
                
    for rxn2 in sbml_2.reactions:
        rxn_id = rxn2.id
        try:
            rxn1 = sbml_1.reactions.get_by_id(rxn_id)
            same_cpd, same_rev =  compare_rxn(rxn1, rxn2)
            if rxn_id not in all_diff:
                if not same_cpd:
                    print("%s use diff cpds" %rxn_id)
                if not same_rev:
                    print("%s diff rev" %rxn_id)
                all_diff.add(rxn_id)
        except KeyError:
            pass

def compare_rxn(rxn1,rxn2):
    same_cpd, same_rev = False, False
    rxn1_dict,rxn2_dict = {}, {}
    for k,v in list(rxn1.metabolites.items()):
        rxn1_dict[k.id] = v
    for k,v in list(rxn2.metabolites.items()):
        rxn2_dict[k.id] = v
    if rxn1_dict == rxn2_dict: 
        same_cpd = True
    if rxn1.reversibility == rxn2.reversibility:

       same_rev = True

    return (same_cpd, same_rev)
        
        
    
    
if __name__ == "__main__":
    main()