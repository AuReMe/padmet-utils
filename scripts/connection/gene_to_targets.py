# -*- coding: utf-8 -*-
"""
Description:
    From a list of genes, get from the linked reactions the list of products.

    R1 is linked to G1, R1 produces M1 and M2.  output: M1,M2. Takes into account reversibility

::

    usage:
        gene_to_targets.py --padmetSpec=FILE --genes=FILE --output=FILE [-v]
    
    option:
        -h --help     Show help
        --padmetSpec=FILE    path to the padmet file
        --genes=FILE   path to the file containing gene ids, one id by line
        --output=FILE    path to the output file containing all tagerts which can by produced by all reactions associated to the given genes
        -v   print info
"""
from padmet.classes import PadmetSpec
import docopt

def main():
    #recovering args
    args = docopt.docopt(__doc__)
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    genes_file = args["--genes"]
    output = args["--output"]
    verbose = args["-v"]
    
    with open(genes_file,'r') as f:
        all_genes = f.read().splitlines()
        nb_genes = len(all_genes)
    all_targets = set()
    count = 0
    for gene_id in all_genes:
        count += 1
        try:
            #check if gene id in padmet
            padmetSpec.dicOfNode[gene_id]
            if verbose: print("%s/%s %s" % (count, nb_genes, gene_id))
        except KeyError:
            if verbose: print("Gene %s not found in padmetSpec" %gene_id)
            continue

            #get all reactions linked to the gene
        try:
            reactions_linked = [padmetSpec.dicOfNode[rlt.id_in] 
            for rlt in padmetSpec.dicOfRelationOut[gene_id]
            if rlt.type == "is_linked_to"]
        except KeyError:
            if verbose: print("the gene %s is not linked to any reaction" %gene_id)
            continue
        if verbose: print("\t%s reactions linked" %len(reactions_linked))

        for reaction_node in reactions_linked:
            reaction_node_id = reaction_node.id
            if verbose: print("\t\t"+reaction_node_id)
            #all products
            products = set([rlt.id_out 
            for rlt in padmetSpec.dicOfRelationIn.get(reaction_node_id, None)
            if rlt.type == "produces"])
            #add the products as targets
            all_targets = all_targets.union(products)
            #if reaction is reversible (reversible or UNKNOWN direction) add reactants
            if reaction_node.misc["DIRECTION"][0] != "LEFT-TO-RIGHT":
                #all reactants
                reactants = set([rlt.id_out 
                for rlt in padmetSpec.dicOfRelationIn.get(reaction_node_id, None)
                if rlt.type == "consumes"])
                all_targets = all_targets.union(reactants)
    
    if len(all_targets) != 0:
        with open(output, 'w') as f:
            f.write("\n".join(all_targets))
        
                
if __name__ == "__main__":
    main()