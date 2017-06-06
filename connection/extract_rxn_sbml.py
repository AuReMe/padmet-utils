# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:

usage:
    extract_rxn_sbml.py --sbml=FILE --output=FILE [--rxn_id=ID]

options:
    -h --help     Show help.
    --reaction_file=FILE    pathname of the file containing the reactions id, 1/line 
    --padmetRef=FILE    pathname of the padmet representing the database.
    --output=FILE    pathname of the file with line = pathway id, all reactions id, reactions ids from reaction file, ratio. sep = "\t"
"""

from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra import *
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    rxn_id = args["--rxn_id"]
    output = args["--output"]

    model=create_cobra_model_from_sbml_file(sbml_file)

    if rxn_id is None:
        try:
            rxn = (rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0).next()
        except StopIteration:
            print("No reaction id given and no reaction with obj coefficient to 1.0, enable to get the reaction")
            exit()
    else:
        rxn = model.reactions.get_by_id(rxn_id)
    
    with open(output, 'w') as f:
        line = ["#reactant/product","#stoichio:compound_id:compart"]
        line = "\t".join(line)+"\n"
        f.write(line)
        reactants = rxn.reactants
        products = rxn.products
        for reactant in reactants:
            stoich = str(abs(rxn.get_coefficient(reactant)))
            reactant_id = reactant.id.replace("DASH","")[:-2]
            compart = reactant.compartment
            print reactant_id+"\t"+compart
            line = ":".join([stoich, reactant_id, compart])
            line = "reactant"+"\t"+line+"\n"
            f.write(line)
        for product in products:
            stoich = str(abs(rxn.get_coefficient(product)))
            product_id = product.id.replace("DASH","")[:-2]
            compart = product.compartment
            line = ":".join([stoich, product_id, compart])
            line = "product"+"\t"+line+"\n"
            f.write(line)

if __name__ == "__main__":
    main()