# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Allows to create news reactions that are not in the database of reference and add
them in a padmet file. First, fill the template (to get the template use the 2nd usage).
Then, use the 1st usage to add the reaction(s) in new_reaction_data to padmetSpec.
The script will check if the choosen id are available and if the compounds used exist or not.
If 1-n compounds don't exist, they will be created. Checking the existence of compounds
in padmetSpec AND padmetRef. 

usage:
    rxn_creator.py --padmetSpec=FILE --new_rxn_data=FILE [--padmetRef=FILE] [--new_padmet=FILE] [-v]
    rxn_creator.py --getTemplate --output=FILE

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathname of the padmet file to update
    --padmetRef=FILE    pathname of the padmet file representing the database of reference
    --new_rxn_data=FILE    pathname to the file containing the data (template obtained from 2nd usage)
    --new_padmet=FILE    new name of padmet after update
    -v   print info

"""
from padmet.padmetSpec import PadmetSpec
from padmet.padmetRef import PadmetRef
from padmet.node import Node
from padmet.relation import Relation
import re
try:
    import docopt
except ImportError:
    print("package docopt needed, use this cmd:\n pip install "
          + "docopt")
    exit()

def main():
    args = docopt.docopt(__doc__)

    get_template = args["--getTemplate"]
    if get_template:
        output = args["--output"]
        getTemplate(output)
        return

    padmetSpec = PadmetSpec(args["--padmetSpec"])
    if args["--padmetRef"] is not None:
        padmetRef  = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    data_file = args["--new_rxn_data"]
    output = args["--new_padmet"]
    if output is None:
        output = args["--padmetSpec"]
    verbose = args["-v"]

    dict_data = {}
    #data_file = "/home/maite/Forge/docker/aureme_share/wenelen/manual_curation/test.csv"
    with open(data_file, 'r') as f:
        data = (line for line in f.read().splitlines() if len(line) != 0 and not line.startswith("#"))
    for line in data:
        #if len of value is 0 then ValueError raised
        try:
            attrib, value = line.split("\t")
        except TypeError:
            continue
        #delete all tags
        if attrib == "reaction_id":
            current_id = value
            dict_data[current_id] = {}
        else:
            try:
                dict_data[current_id][attrib] .append(value)
            except KeyError:
                dict_data[current_id][attrib] = [value]

    if verbose: print len(dict_data.keys()),"reactions to add"
    for reaction_id, reaction_data in dict_data.iteritems():
        if verbose: print("check if the id "+reaction_id+" is already used")
        if reaction_id in padmetSpec.dicOfNode.keys():
            print("the id : %s is already associated to an other reaction in padmetSpec, choose an other" %reaction_id)
            continue
        if padmetRef is not None and reaction_id in padmetRef.dicOfNode.keys():
            print("the id : %s is already associated to an other reaction in padmetRef, choose an other" %reaction_id)
            continue

        if verbose: print("adding reaction %s" %reaction_id)
        reaction_rev = reaction_data["reversible"][0].lower()
        if reaction_rev == "true":
            reaction_rev = "REVERSIBLE"
        elif reaction_rev == "false":
            reaction_rev = "LEFT-TO-RIGHT"
        else:
            print("Please choose a value in ['true','false'] for the reversibility of the reaction: %s" %reaction_id)
            continue
        comment = reaction_data["comment"][0]
        node_misc = {"DIRECTION":[reaction_rev], "SOURCE":["manual"], "COMMENT":[comment]}
        reaction_node = Node("reaction", reaction_id)
        reaction_node.misc = node_misc
        padmetSpec.dicOfNode[reaction_id] = reaction_node

        gene_assoc = reaction_data["linked_gene"][0]
        if len(gene_assoc) != 0:
            if verbose: print("adding gene assocations")
            gene_assoc = re.sub("\(|\)|\s","",gene_assoc).replace("and","or")
            gene_assoc = gene_assoc.split("or")
            if verbose: print("nb gene_assoc: %s" %len(gene_assoc))
            for gene in gene_assoc:
                try:
                    padmetSpec.dicOfNode[gene]
                except KeyError:
                    gene_node = Node("gene", gene)
                    padmetSpec.dicOfNode[gene] = gene_node
                assoc_rlt = Relation(reaction_id, "is_linked_to", gene)
                padmetSpec._addRelation(assoc_rlt)

        if verbose: print("check if all metabolites are already in the network")
        try:
            for reactant_data in reaction_data["reactant"]:
                stoechio, metabo_id, compart = reactant_data.split(":")
                stoechio = stoechio.replace(",",".") #in case comma for sep
                try:
                    padmetSpec.dicOfNode[metabo_id]
                except KeyError:
                    if verbose: print("%s not in the network" %metabo_id)
                    try:
                        if padmetRef is not None:
                            if verbose: print("Try to copy from dbref")
                            padmetSpec._copyNodeExtend(padmetRef, metabo_id)
                        else:
                            raise KeyError
                    except KeyError:
                        if padmetRef is not None and verbose:
                            print("%s not in the padmetRef" %metabo_id)
                        if verbose:
                            print("creating a new compound")
                        metabo_node = Node("compound", metabo_id)
                        padmetSpec.dicOfNode[metabo_id] = metabo_node
                        if verbose: print("new compound created: id = %s" %metabo_id)
                rlt = Relation(reaction_id,"consumes",metabo_id)
                rlt.misc.update({"STOICHIOMETRY":[stoechio],"COMPARTMENT":[compart]})
                padmetSpec._addRelation(rlt)
        except KeyError:
            if verbose: print("No reactants defined")

        try:
            for product_data in reaction_data["product"]:
                stoechio, metabo_id, compart = product_data.split(":")
                stoechio = stoechio.replace(",",".") #in case comma for sep
                try:
                    padmetSpec.dicOfNode[metabo_id]
                except KeyError:
                    if verbose: print("%s not in the network" %metabo_id)
                    try:
                        if padmetRef is not None:
                            if verbose: print("Try to copy from dbref")
                            padmetSpec._copyNodeExtend(padmetRef, metabo_id)
                        else:
                            raise KeyError
                    except KeyError:
                        if padmetRef is not None and verbose:
                            print("%s not in the padmetRef" %metabo_id)
                        if verbose:
                            print("creating a new compound")
                        metabo_node = Node("compound", metabo_id)
                        padmetSpec.dicOfNode[metabo_id] = metabo_node
                        print("new compound created: id = %s" % metabo_id)
                rlt = Relation(reaction_id,"produces",metabo_id)
                rlt.misc.update({"STOICHIOMETRY":[stoechio],"COMPARTMENT":[compart]})
                padmetSpec._addRelation(rlt)
        except KeyError:
            if verbose: print("No products defined")

    if verbose: print("Creating output: %s" % output)
    padmetSpec.generateFile(output)

def getTemplate(output):
    with open(output, 'w') as f:
        line = "\t".join(["reaction_id","my_rxn"])+"\n"
        f.write(line)
        line = "\t".join(["comment","reaction added for X reason"])+"\n"
        f.write(line)
        line = "\t".join(["reversible","false"])+"\n"
        f.write(line)
        line = "\t".join(["linked_gene","(gene_a or gene_b) and gene_c"])+"\n"
        f.write(line)
        line = "\t".join(["#reactant/product","#stoichio:compound_id:compart"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","1.0:compound_a:c"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","2.0:compound_b:c"])+"\n"
        f.write(line)
        line = "\t".join(["product","1.0:compound_c:c"])+"\n"
        f.write(line)
        f.write("\n")
        line = "\t".join(["reaction_id","my_rxn_2"])+"\n"
        f.write(line)
        line = "\t".join(["comment","reaction added for X reason"])+"\n"
        f.write(line)
        line = "\t".join(["reversible","true"])+"\n"
        f.write(line)
        line = "\t".join(["linked_gene",""])+"\n"
        f.write(line)
        line = "\t".join(["#reactant/product","#stoichio:compound_id:compart"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","1.0:compound_a:c"])+"\n"
        f.write(line)
        line = "\t".join(["reactant","2.0:compound_d:c"])+"\n"
        f.write(line)
        line = "\t".join(["product","1.0:compound_c:c"])+"\n"
        f.write(line)

if __name__ == "__main__":
    main()