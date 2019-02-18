# -*- coding: utf-8 -*-
"""
Description:
    convert GBK to FAA with Bio package

::
    
    usage:
        gbk_to_faa.py    --gbk=FILE --output=FILE [--qualifier=STR] [-v]
    
    option:
        -h --help    Show help.
        --gbk=FILE    path to the gbk file.
        --output=FILE    path to the output, a FAA file.
        --qualifier=STR    the qualifier of the gene id [default: locus_tag].
        -v   print info
"""
import docopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import protein

def main():
    args = docopt.docopt(__doc__)
    gbk_file = args["--gbk"]
    faa_file = args["--output"]
    qualifier = args["--qualifier"]
    verbose = args["-v"]
    fasta_records = []
    with open(gbk_file, "r") as gbk:
        for seq_record in SeqIO.parse(gbk, "genbank"):
            seq_feature_cds = (seq_feature for seq_feature in seq_record.features if seq_feature.type == "CDS")
            for seq_feature in seq_feature_cds:
                try:
                    fasta_records.append(SeqRecord(Seq(seq_feature.qualifiers['translation'][0], protein), id=seq_feature.qualifiers[qualifier][0], description=""))
                except KeyError:
                    if verbose:
                        try:
                            print("locus without Translation: "+seq_feature.qualifiers['locus_tag'][0])
                        except KeyError:
                            print("locus without Translation: "+seq_feature.qualifiers.get('gene',["Unknown"])[0])

    SeqIO.write(fasta_records, faa_file, "fasta")

if __name__ == "__main__":
    main()