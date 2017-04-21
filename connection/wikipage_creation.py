# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Create the files containing the wikicode for each page of the wiki relative to a
padmetSpec based on a database padmetRef. All the pages are stored in the folder
output_dir. The padmetRef is required to extract the information about all the reactions
of a given pathway. In the case of BIGG db, it's possible to use the padmetSpec as
padmetRef because there is no pathways.

usage:
    wikiPage_creation.py --padmetSpec=FILE --output_dir=DIR [--padmetRef=FILE] [-v]

option:
    -h --help    Show help.
    --padmetRef=FILE    pathname of the padmet of reference.
    --padmetSpec=FILE    pathname of the padmet of interest.
    --output_dir=DIR    pathname of the directory were to stock all the wikiPages.
    -v   print info
"""
import padmet.wikiGenerator as wg
from time import time
import docopt

def main():
    args = docopt.docopt(__doc__)
    
    padmetRef_file = args["--padmetRef"]
    padmetSpec_file = args["--padmetSpec"]
    output_dir = args["--output_dir"]
    verbose = args["-v"]

    chronoDepart = time()
    wg.create_all_wikiPages(padmetSpec_file, output_dir, padmetRef_file, verbose)
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print "done in: ", chrono, "s !"

if __name__ == "__main__":
    main()

