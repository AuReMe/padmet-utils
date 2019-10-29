#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Create sbml from sbml. Use it to change sbml level.

::
    usage:
        sbml_to_sbml.py --input=FILE/FOLDER --output=FILE/FOLDER --new_sbml_lvl=STR [--cpu=INT]
    
    option:
        -h --help    Show help.
        --input=FILE    path of the sbml file/folder to convert into sbml
        --output=FILE    path of the sbml file/folder to generate.
        --new_sbml_lvl=STR    level of the new sbml.
        --cpu=FILE    number of cpu used [default: 1].
        -v   print info.
"""

import docopt
from padmet.utils.connection import sbml_to_sbml

def main():
    args = docopt.docopt(__doc__)
    input_sbml = args["--input"]
    output_sbml = args["--output"]
    new_sbml_level = args["--new_sbml_lvl"]
    cpu = args['--cpu']

    sbml_to_sbml.from_sbml_to_sbml(input_sbml, output_sbml, new_sbml_level, cpu)

if __name__ == "__main__":
    main()