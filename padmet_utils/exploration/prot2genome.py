#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
usage:
    prot2genome --query_faa=FILE --query_ids=FILE/STR --subject_gbk=FILE --subject_fna=FILE --subject_faa=FILE --output_folder=FILE [--cpu=INT] [blastp] [tblastn] [exonerate]  [debug]

options:
    --query_faa=FILE #TODO. 
    --query_ids=FILE/STR #TODO. 
    --subject_gbk=FILE #TODO. 
    --subject_fna=FILE #TODO. 
    --subject_faa=FILE #TODO. 
    --output_folder=FILE #TODO. 
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu). [default: 1]
    blastp #TODO. 
    tblastn #TODO. 
    exonerate #TODO. 
    debug #TODO. 
    
"""

from Bio.Blast.Applications import NcbiblastpCommandline, NcbitblastnCommandline
from Bio import SearchIO
from Bio import SeqIO
import os
import csv
import subprocess
from multiprocessing import Pool
import docopt

#Utilise gbk to faa, ajouter une option pour faire un fichier fasta une sequence.
#pident a ajouter
#TODO:
#Seq to gbk, match to faa


def main():
    """
    """
    args = docopt.docopt(__doc__)
    query_faa = args["--query_faa"]
    query_ids = args["--query_ids"]
    subject_gbk = args["--subject_gbk"]
    subject_fna = args["--subject_fna"]
    subject_faa = args["--subject_faa"]
    output_folder = args["--output_folder"]
    cpu =  int(args["--cpu"])
    blastp = args["blastp"]
    tblastn = args["tblastn"]
    exonerate = args["exonerate"]
    debug = args["debug"]

    #if query_ids is a file, extract all query_ids (1 by line)
    #else, query_ids must represent query_ids sep by ';'
    if os.path.isfile(query_ids):
        with open(query_ids, 'r') as f:
            all_query_seq_ids = f.read().splitlines()
    else:
        all_query_seq_ids = query_ids.split(";")

    print("%s query ids to search" %len(all_query_seq_ids))
    #If GBK given and no fna, create fna from genbank
    if not os.path.exists(subject_fna):
        SeqIO.convert(subject_gbk, "genbank", subject_fna, "fasta")
    #Do the same for gbk to fasta ? (must be with padmet-utils for isoforms)
    #TODO

    #output file: header: see analysis_header from analysisOutput function
    analysis_output = os.path.join(output_folder, "blast_analysis.csv")

    pool = Pool(cpu)
    #list of dict to give to dictWritter
    all_analysis_result = []
    #Create a list of dict, each dict is the arg given to runAllAnalysis fct
    all_dict_args = []
    for query_seq_id in all_query_seq_ids:
        dict_args = {"query_seq_id": query_seq_id, "query_faa": query_faa, "subject_faa": subject_faa, "subject_fna": subject_fna, "output_folder": output_folder, "blastp": blastp, "tblastn": tblastn, "exonerate": exonerate, "debug": debug}
        all_dict_args.append(dict_args)
    #Run runAllAnalysis in multiproccess.
    mp_results = pool.map(runAllAnalysis, all_dict_args)
    for _list in mp_results:
        all_analysis_result += _list
    print("Creating output analysis: %s" %analysis_output)
    #Create output file
    analysisOutput(all_analysis_result, analysis_output)


def runBlastp(query_seq_faa, subject_faa, header=["sseqid", "evalue", "bitscore"], debug=False):
    """
    """
    print("\tRunning Blastp %s vs %s" %(os.path.basename(query_seq_faa), os.path.basename(subject_faa)))
    outfmt_arg = '"%s %s"'%(6, " ".join(header))
    output = NcbiblastpCommandline(query=query_seq_faa, subject=subject_faa, evalue=1e-10, outfmt=outfmt_arg)()[0]
    output = [line.split("\t") for line in output.splitlines()]
    blastp_result = {}
    count = 0
    for line in output:
        count += 1
        blastp_result[count] = {}
        for (index, h) in enumerate(header):
            blastp_result[count]["blastp_"+h] = line[index]
    try:
        max_bitscore = max([float(hsp["blastp_bitscore"]) for hsp in blastp_result.values()])
    except ValueError:
        max_bitscore = "N.A"
    if debug:
        print("\t\tAll hsp from blastp:")
        if blastp_result.values():
            for hsp in blastp_result.values():
                print(hsp)
            print("\t\tMax biscore: %s"%max_bitscore)
        else:
            print("\t\tNo HPS")
    try:
        result = [hsp for hsp in blastp_result.values() if float(hsp["blastp_bitscore"]) == max_bitscore][0]
    except IndexError:
        result = {}
    return result

def runTblastn(query_seq_faa, subject_fna, header=["sseqid", "evalue", "bitscore"], debug=False):
    """
    """    
    print("\tRunning tBlastn %s vs %s" %(os.path.basename(query_seq_faa), os.path.basename(subject_fna)))
    outfmt_arg = '"%s %s"'%(6, " ".join(header))
    output = NcbitblastnCommandline(query=query_seq_faa, subject=subject_fna, evalue=1e-10, outfmt=outfmt_arg)()[0]
    output = [line.split("\t") for line in output.splitlines()]
    tblastn_result = {}
    count = 0
    for line in output:
        count += 1
        tblastn_result[count] = {}
        for (index, h) in enumerate(header):
            tblastn_result[count]["tblastn_"+h] = line[index]
    try:
        max_bitscore = max([float(hsp["tblastn_bitscore"]) for hsp in tblastn_result.values()])
    except ValueError:
        max_bitscore = "N.A"
    if debug:
        print("\t\tAll hsp from tblastn:")
        for hsp in tblastn_result.values():
            print(hsp)
        print("\t\tMax biscore: %s"%max_bitscore)
    try:           
        result = [hsp for hsp in tblastn_result.values() if float(hsp["tblastn_bitscore"]) == max_bitscore][0]
    except IndexError:
        result = {}
    return result

def runExonerate(query_seq_faa, sseq_seq_faa, output, debug=False):
    """
    """
    print("\tRunning Exonerate %s vs %s" %(os.path.basename(query_seq_faa), os.path.basename(sseq_seq_faa)))
    exonerate_result = {}
    exonerate_path = "/home/maite/exonerate-2.2.0-x86_64/bin/exonerate"
    cmd_args = '{0} --model protein2genome {1} {2} --score 500 --showquerygff True  \
    --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" >> {3}'.format(exonerate_path, query_seq_faa, sseq_seq_faa, output)
    subprocess.run([cmd_args], shell = True)
    try:
        exonerate_raw_output = list(SearchIO.parse(output, 'exonerate-text'))
       
        best_hsp = exonerate_raw_output[0][0][0]
        for qresult in exonerate_raw_output:
            for hit in qresult:
                for hsp in hit:
                    if hsp.score > best_hsp.score:
                        best_hsp = hsp
        if debug:
            print("\t\tAll hsp from exonerate:")
            for qresult in exonerate_raw_output:
                for hit in qresult:
                    for hsp in hit:
                        print(hsp)
            print("\t\tMax score: %s"%best_hsp.score)
        exonerate_result = {"exonerate_score": best_hsp.score, "exonerate_hit_range": best_hsp.hit_range}
    except RuntimeError:
        if debug:
            print("\t\tNo HSP")
    return exonerate_result


def runAllAnalysis(dict_args):
    """
    """
    query_seq_id = dict_args["query_seq_id"]
    query_faa = dict_args["query_faa"]
    subject_faa = dict_args["subject_faa"]
    subject_fna = dict_args["subject_fna"]
    output_folder = dict_args["output_folder"]
    blastp = dict_args["blastp"]
    tblastn = dict_args["tblastn"]
    exonerate = dict_args["exonerate"]
    debug = dict_args["debug"]
        
    
    with open(query_faa, "r") as faa:
        query_seqs = [seq_record for seq_record in SeqIO.parse(faa, "fasta") if seq_record.id.startswith(query_seq_id+"_isoform") or seq_record.id == query_seq_id]
    if len(query_seqs) > 1:
        print("/!\ Isoforms found for %s: %s"%(query_seq_id, [i.name for i in query_seqs]))
    analysis_result = list()
    for query_seq in query_seqs:
        query_seq_id = query_seq.id
        query_seq_faa = os.path.join(output_folder,query_seq_id+".faa")
        if not os.path.exists(query_seq_faa):
            SeqIO.write(query_seq, query_seq_faa, "fasta")
        current_result ={"query_seq_id": query_seq_id,}
    
        # Run BLASTP and parse the output
        if blastp:
            blastp_result = runBlastp(query_seq_faa, subject_faa, debug=debug)
            if blastp_result:
                current_result.update(blastp_result)
        # Run TBLASTN and parse the output
        if tblastn:
            tblastn_result = runTblastn(query_seq_faa, subject_fna, debug=debug)
            if tblastn_result:
                current_result.update(tblastn_result)
        
        # Run exonerate, parse output
        if exonerate:
            if not tblastn:
                IOError("must run tblastn to be able to run exonerate")
            _tblastn_hit = True
            try:
                exonerate_target_id = current_result["tblastn_sseqid"]
            except KeyError:
                _tblastn_hit = False
                if debug:
                    print("No hit from tblastn, can't run exonerate")
            if _tblastn_hit:
                sseq_seq_faa = os.path.join(output_folder,exonerate_target_id+".fna")
                if not os.path.exists(sseq_seq_faa):
                    with open(subject_fna, "r") as fna:
                        sseq_seq = [seq_record for seq_record in SeqIO.parse(fna, "fasta") if seq_record.id == exonerate_target_id][0]
                        SeqIO.write(sseq_seq, sseq_seq_faa, "fasta")
                exonerate_output = os.path.join(output_folder, "exonerate_output_%s_vs_%s.txt"%(query_seq_id, exonerate_target_id))
                exonerate_result = runExonerate(query_seq_faa, sseq_seq_faa, exonerate_output, debug=debug)
                current_result.update(exonerate_result)
        analysis_result.append(current_result)

    return analysis_result

def analysisOutput(analysis_result, analysis_output):
    """
    """
    analysis_header = ["query_seq_id", "blastp_sseqid", "blastp_evalue", "blastp_bitscore",  "tblastn_sseqid", "tblastn_evalue", "tblastn_bitscore", "exonerate_score", "exonerate_hit_range"]
    with open(analysis_output,"w") as csvfile:
        dict_writer = csv.DictWriter(csvfile, fieldnames=analysis_header, delimiter="\t")
        dict_writer.writeheader()            
        dict_writer.writerows(analysis_result)    

if __name__ == "__main__":
    main()    

