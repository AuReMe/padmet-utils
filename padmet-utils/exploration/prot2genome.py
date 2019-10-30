#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
usage:
    prot2genome --query_faa=FILE --query_ids=FILE/STR --subject_gbk=FILE --subject_fna=FILE --subject_faa=FILE --output_folder=FILE [--cpu=INT] [blastp] [tblastn] [debug]
    prot2genome --query_faa=FILE --query_ids=FILE/STR --subject_gbk=FILE --subject_fna=FILE --subject_faa=FILE --output_folder=FILE --exonerate=PATH  [--cpu=INT] [blastp] [tblastn] [debug]
    prot2genome --padmet=FOLDER --output=FOLDER
    prot2genome --studied_organisms=FOLDER --output=FOLDER
    prot2genome --run=FOLDER --exonerate=PATH  --padmetRef=FILE [--cpu=INT] [blastp] [tblastn] [debug]

    From aucome run fromAucome():
        -1. Extract specifique reactions in spec_reactions folder with extractReactions()
        -2. Extract genes from spec_reactions files with extractGenes()
        -3. Run tblastn + exonerate with runAllAnalysis()
    
options:
    --query_faa=FILE #TODO. 
    --query_ids=FILE/STR #TODO. 
    --subject_gbk=FILE #TODO. 
    --subject_fna=FILE #TODO. 
    --subject_faa=FILE #TODO. 
    --output_folder=FILE #TODO. 
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu). [default: 1]
    --exonerate=PATH #TODO. 
    blastp #TODO. 
    tblastn #TODO. 
    debug #TODO.
    
"""

from padmet.utils.exploration import prot2genome
import subprocess
import docopt

#Utilise gbk to faa, ajouter une option pour faire un fichier fasta une sequence.
#pident a ajouter
#TODO:
#Seq to gbk, match to faa


def main():
    """
    """
    global exonerate_path
    args = docopt.docopt(__doc__)
    exonerate_path = args["--exonerate"]
    if exonerate_path:
        print("Test running exonerate...")
        subprocess.run([exonerate_path], shell = True)
        exonerate = True

    run_folder = args["--run"]
    cpu =  int(args["--cpu"])
    padmetRef = args["--padmetRef"]
    blastp = args["blastp"]
    tblastn = args["tblastn"]
    debug = args["debug"]

    if run_folder:
        prot2genome.fromAucome(run_folder, cpu, padmetRef, blastp, tblastn, exonerate, debug)

if __name__ == "__main__":
    main()    

