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
1./ extract all data from sgs output
2./ If needed, convert reaction id to metacyc id

usage:
    enhanced_sgs_output.py --sgs_output=FILE --padmetRef=FILE --output=FILE [--db=ID] [--mnx=FILE] [-v]

options:
    -h --help     Show help.
    --sgs_output=FILE    pathname of a sgs run' result
    --db=ID    database origin of the reactions. will be converted to metacyc
    --padmetRef=FILE    pathanme to padmet file corresponding to the database of reference (the repair network)
    --mnx=FILE    pathanme to metanetx file for reactions (reac_xref.tsv)
    --output=FILE    pathname to tsv output file
"""
import csv
import re
import docopt
from itertools import izip_longest

def main():
    args = docopt.docopt(__doc__)
    
    sgs_output = args["--sgs_output"]
    output = args["--output"]
    verbose = args["-v"]
    output = args["--output"]
    
    verbose = True
    sgs_output = "/home/maite/Documents/data/E.faecalis_data/data_askomics/efaecalis_sgs_final.txt"
    output = "/home/maite/Documents/data/E.faecalis_data/data_askomics/sgs.txt"
    sgs_dict = {}
    #1: extract data from shogen output.
    count = 1
    if verbose: print("Reading sgs output")
    with open(sgs_output, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            genes, reactions = row["gene_set"], row["reaction_set"]
            genes = re.sub("ac:|nc:|na:|\"","",genes).split(" ")
            reactions = re.sub("nc:|\"|R_","",reactions).split(" ")
            sgs_id = "sgs-"+str(count)+"_"+str(len(genes))
            count += 1
            sgs_dict[sgs_id] = {"genes": genes, "reactions": reactions,\
            "start_position": row["start_position"], "end_position": row["end_position"],\
            "length": row["length"], "density": str(round(float(row["density"]),2)),"chr_id": row["chromosome_id"]}
    
    
    if verbose: print("creating output %s" %output)
    with open(output, 'w') as f:
        header = ["SGS","contains@gene","concerns@reaction","start_position","end_position","length","density","chromosome_id"]
        f.write("\t".join(header)+"\n")
        for sgs_id, data in sgs_dict.items():
            all_data_ziped = izip_longest(data["genes"],data["reactions"],[data["start_position"]],[data["end_position"]],[data["length"]],[data["density"]],[data["chr_id"]])
            for line_data in all_data_ziped:
                line = [sgs_id]
                line += [x if x else "" for x in line_data]
                f.write("\t".join(line)+"\n")
                        
if __name__ == "__main__":
    main()                  
