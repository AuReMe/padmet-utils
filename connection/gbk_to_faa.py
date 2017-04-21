# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
convert GBK to FAA with Bio package

usage:
    gbk_to_faa.py    --gbk=FILE --output=FILE [--qualifier=STR] [-v]

option:
    -h --help    Show help.
    --gbk=FILE    pathname to the gbk file.
    --output=FILE    pathename to the output, a FAA file.
    --qualifier=STR    the qualifier of the gene id [default: locus_tag].
    -v   print info


"""
from Bio import SeqIO
import docopt

def main():
    args = docopt.docopt(__doc__)
    gbk_file = args["--gbk"]
    faa_file = args["--output"]
    qualifier = args["--qualifier"]
    verbose = args["-v"]
    dict_faa = {}
    with open(gbk_file, "rU") as gbk:
        for seq_record in SeqIO.parse(gbk, "genbank"):
            seq_feature_cds = (seq_feature for seq_feature in seq_record.features if seq_feature.type == "CDS")
            for seq_feature in seq_feature_cds:
                try:
                    dict_faa[seq_feature.qualifiers[qualifier][0]] = seq_feature.qualifiers['translation'][0]
                except KeyError:
                    if verbose:
                        try:
                            print ("locus without Translation: "+seq_feature.qualifiers['locus_tag'][0])
                        except KeyError:
                            print ("locus without Translation: "+seq_feature.qualifiers.get('gene',["Unknown"])[0])

    with open(faa_file,'w') as f:
        for locus_tag, seq in dict_faa.iteritems():
            f.write(">%s\n%s\n" % (locus_tag,seq))

if __name__ == "__main__":
    main()