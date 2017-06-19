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
From a list of genes, recovere from the linked reactions the list of products.
R1 is linked to G1, R1 produces M1 and M2.  output: M1,M2. Takes into account reversibility

usage:
    gene_to_targets.py --padmetSpec=FILE --genes=FILE --output=FILE [-v]

option:
    -h --help     Show help
    --padmetSpec=FILE    pathanme to the padmet file
    --genes=FILE   pathname to the file containing gene ids, one id by line
    --output=FILE    pathname to the output file containing all tagerts which can by produced by all reactions associated to the given genes
    -v   print info

"""
from padmet.padmetSpec import PadmetSpec
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
            print("Gene %s not found in padmetSpec" %gene_id)
            continue

            #get all reactions linked to the gene
        try:
            reactions_linked = [padmetSpec.dicOfNode[rlt.id_in] 
            for rlt in padmetSpec.dicOfRelationOut[gene_id]
            if rlt.type == "is_linked_to"]
        except KeyError:
            print("the gene %s is not linked to any reaction" %gene_id)
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