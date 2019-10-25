#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Contains all necessary functions to generate wikiPages from a padmet file and update 
    a wiki online. Require WikiManager module (with wikiMate,Vendor)

::

    usage:
        wikiGenerator.py --padmet=FILE/DIR --output=DIR --wiki_id=STR [--database=STR] [--padmetRef=FILE] [--log_file=FILE] [-v]
        wikiGenerator.py --aureme_run=DIR --padmetSpec=ID -v
    
    options:
        -h --help     Show help.
        --padmet=FILE    path to padmet file.
        --output=DIR    path to folder to create with all wikipages in subdir.
        --wiki_id=STR    id of the wiki.
        --padmetRef=FILE    path to padmet of reference, ex: metacyc_xx.padmet, if given, able to calcul pathway rate completion.
        --log_file=FILE    log file from an aureme run, use this file to create a wikipage with all the command used during the aureme run.
        --aureme_run=DIR    can use an aureme run as input, will use from config file information for model_id and log_file and padmetRef.
        -v    print info.
"""
from padmet.classes import PadmetRef
from padmet.utils.connection import wikiGenerator
import docopt

def main(): 
    #files to upload: folder genomic_data, all sbml in output ortho, annot, external, seeds, targets
    args = docopt.docopt(__doc__)
    padmet = args["--padmet"]
    verbose = args["-v"]
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    wiki_id = args["--wiki_id"]
    output = args["--output"]
    log_file = args["--log_file"]
    database = args["--database"]
    
    wikiGenerator.wikiGenerator(padmet, output, wiki_id, padmetRef, database, log_file, verbose)

if __name__ == "__main__":
    main()