# -*- coding: utf-8 -*-
"""
Description:
    Before running orthology based reconstruction it is necessary to check if the metabolic network 
    and the proteom of the model organism use the same ids for genes (or at least more than a given cutoff). 
    To only check this. Use the 2nd usage.

    If the genes ids are not the same, it is necessary to use a dictionnary of genes ids associating
    the genes ids from the proteom to the genes ids from the metabolic network.

    To create the correct proteom from the dictionnnary, use the 3nd usage
    Finnaly by using the 1st usage, it is possible to:

        1/ Check model_faa and model_metabolic for a given cutoff

        2/ if under the cutoff, convert model_faa to the correct one with dict_ids_file

        3/ if still under, SystemExit()

::
   
    usage:
        check_orthology_input.py    --model_metabolic=FILE    --model_faa=FILE    [--cutoff=FLOAT] [--dict_ids_file=FILE] --output=FILE    [-v]
        check_orthology_input.py    --model_metabolic=FILE    --model_faa=FILE    [--cutoff=FLOAT] [-v]
        check_orthology_input.py    --model_faa=FILE    --dict_ids_file=FILE    --output=FILE [-v]
    
    option:
        -h --help    Show help.
        --model_metabolic=FILE    pathname to the metabolic network of the model (sbml).
        --model_faa=FILE    pathname to the proteom of the model (faa)
        --cutoff=FLOAT    cutoff [0:1] for comparing model_metabolic and model_faa. [default: 0.70]. 
        --dict_ids_file=FILE    pathname to the dict associating genes ids from the model_metabolic to the model_faa. line = gene_id_in_metabolic_network\tgene_id_in_faa
        --output=FILE    output of get_valid_faa (a faa) or get_dict_ids (a dictionnary of gene ids in tsv)
        -v   print info
"""
import docopt
from padmet.utils.connection import check_orthology_input

def main():
    args = docopt.docopt(__doc__)
    model_metabolic = args["--model_metabolic"]
    model_faa = args["--model_faa"]
    dict_ids_file = args["--dict_ids_file"]
    output = args["--output"]
    verbose = args["-v"]
    cutoff = float(args["--cutoff"])
    check_orthology_input.check_orthology_input(model_metabolic, model_faa, dict_ids_file, output, verbose, cutoff)



if __name__ == "__main__":
    main()

