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
    with open(gbk_file, "r") as gbk:
        seqs = [seq for seq in SeqIO.parse(gbk, "genbank")]
        for seq_record in SeqIO.parse(gbk, "genbank"):
            seq_feature_cds = (seq_feature for seq_feature in seq_record.features if seq_feature.type == "CDS")
            for seq_feature in seq_feature_cds:
                try:
                    dict_faa[seq_feature.qualifiers[qualifier][0]] = seq_feature.qualifiers['translation'][0]
                except KeyError:
                    if verbose:
                        try:
                            print("locus without Translation: "+seq_feature.qualifiers['locus_tag'][0])
                        except KeyError:
                            print("locus without Translation: "+seq_feature.qualifiers.get('gene',["Unknown"])[0])

    with open(faa_file,'w') as f:
        for locus_tag, seq in dict_faa.items():
            f.write(">%s\n%s\n" % (locus_tag,seq))

if __name__ == "__main__":
    main()