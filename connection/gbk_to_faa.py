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
    with open(gbk_file, "rU") as gbk:
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