# -*- coding: utf-8 -*-
"""
Description:
    extract all reactions from a sbml file to the form used in aureme for curation.

::

    usage:
        sbml_to_curation_form.py --sbml=FILE --output=FILE --comment=STR [--rxn_id=ID]
    
    options:
        -h --help     Show help.
        --sbml=FILE    path of the sbml.
        --output=FILE    form containing the reaction extracted, form used for manual curation in aureme.
        --rxn_id=FILE    id of one reaction or n reactions sep by ';', if None try to extract the reaction with objective coefficient == 1.
        --comment=STR    comment associated to the reactions in the form. Used to track sources of curation in aureme.
"""

import docopt
from cobra.io.sbml import create_cobra_model_from_sbml_file
from padmet.utils.sbmlPlugin import convert_from_coded_id

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    rxn_id = args["--rxn_id"]
    output = args["--output"]
    comment = args["--comment"]

    model=create_cobra_model_from_sbml_file(sbml_file)

    if not args["--rxn_id"]:
        try:
            rxn = next(rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0)
        except StopIteration:
            print("No reaction id given and no reaction with obj coefficient to 1.0, enable to get the reaction")
            exit()
    else:
        rxns = args["--rxn_id"].split(";")
        rxns_list = []
        for rxn_id in rxns:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
            except KeyError:
                try:
                    rxn_id = rxn_id[2:]
                    rxn = model.reactions.get_by_id(rxn_id)
                except KeyError:
                    print("Can't find reaction %s" %rxn_id)
                    exit()
            rxns_list.append(rxn)
    
    with open(output, 'w') as f:
        for rxn in rxns_list:
            rxn_id = convert_from_coded_id(rxn.id)[0]
            line = ["reaction_id",rxn_id]
            line = "\t".join(line)+"\n"
            f.write(line)
            line = ["comment",comment]
            line = "\t".join(line)+"\n"
            f.write(line)
            if rxn.reversibility:
                line = ["reversible","true"]
            else:
                line = ["reversible","false"]
            line = "\t".join(line)+"\n"
            f.write(line)
            #check if have gene assoc
            try:
                k = [i for i in list(rxn.notes.keys()) if i.lower().startswith("gene")][0]
                gene_assoc = rxn.notes[k][0]
                line = ["linked_gene", gene_assoc]
            except IndexError:
                line = ["linked_gene", ""]
            line = "\t".join(line)+"\n"
            f.write(line)
            line = ["#reactant/product","#stoichio:compound_id:compart"]
            line = "\t".join(line)+"\n"
            f.write(line)
            reactants = rxn.reactants
            products = rxn.products
            for reactant in reactants:
                stoich = str(abs(rxn.get_coefficient(reactant)))
                reactant_id = convert_from_coded_id(reactant.id)[0]
                compart = reactant.compartment
                line = ":".join([stoich, reactant_id, compart])
                line = "reactant"+"\t"+line+"\n"
                f.write(line)
            for product in products:
                stoich = str(abs(rxn.get_coefficient(product)))
                product_id = convert_from_coded_id(product.id)[0]
                compart = product.compartment
                line = ":".join([stoich, product_id, compart])
                line = "product"+"\t"+line+"\n"
                f.write(line)
            f.write("\n")

if __name__ == "__main__":
    main()