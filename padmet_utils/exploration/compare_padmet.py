# -*- coding: utf-8 -*-
"""
Description:
    #Compare 1-n padmet and create a folder output with files:
    genes.csv:
        fieldnames = [gene, padmet_a, padmet_b, padmet_a_rxn_assoc, padmet_b_rxn_assoc]
        line = [gene-a, 'present' (if in padmet_a), 'present' (if in padmet_b), rxn-1;rxn-2 (names of reactions associated to gene-a in padmet_a), rxn-2]
    reactions.csv:
        fieldnames = [reaction, padmet_a, padmet_b, padmet_a_genes_assoc, padmet_b_genes_assoc, padmet_a_formula, padmet_b_formula]
        line = [rxn-1, 'present' (if in padmet_a), 'present' (if in padmet_b), 'gene-a;gene-b; gene-a, 'cpd-1 + cpd-2 => cpd-3', 'cpd-1 + cpd-2 => cpd-3']
    pathways.csv:
        fieldnames = [pathway, padmet_a_completion_rate, padmet_b_completion_rate, padmet_a_rxn_assoc, padmet_b_rxn_assoc]
        line = [pwy-a, 0.80, 0.30, rxn-a;rxn-b; rxn-a]
    compounds.csv:
        fieldnames = ['metabolite', padmet_a_rxn_consume, padmet_a_rxn_produce, padmet_b_rxn_consume, padmet_rxn_produce]
        line = [cpd-1, rxn-1,'',rxn-1,'']

::
    
    usage:
        compare_padmet.py --padmet=FILES/DIR --output=DIR [--padmetRef=FILE] [-v]
    
    option:
        -h --help    Show help.
        --padmet=FILES/DIR    pathname of the padmet files, sep all files by ',', ex: /path/padmet1.padmet;/path/padmet2.padmet OR a folder
        --output=DIR    pathname of the output folder
        --padmetRef=FILE    pathanme of the database ref in padmet
"""
from padmet.classes import PadmetRef
from padmet.utils.exploration import compare_padmet

import docopt

def main():
    args = docopt.docopt(__doc__)
    output = args["--output"]
    verbose = args["-v"]
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    padmet_path = args["--padmet"]
    compare_padmet.compare_padmet(padmet_path, output, padmetRef, verbose)

if __name__ == "__main__":
    main()