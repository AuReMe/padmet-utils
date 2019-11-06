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

from padmet.utils.connection import gbk_to_faa


def main():
    args = docopt.docopt(__doc__)
    gbk_file = args["--gbk"]
    output = args["--output"]
    qualifier = args["--qualifier"]
    verbose = args["-v"]
    gbk_to_faa.gbk_to_faa(gbk_file, output, qualifier, verbose)


if __name__ == "__main__":
    main()