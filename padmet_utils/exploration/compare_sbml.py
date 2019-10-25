# -*- coding: utf-8 -*-
"""
Description:
    compare reactions in two sbml.

    Returns if a reaction is missing

    And if a reaction with the same id is using different species or different reversibility

::

    usage:
        compare_sbml.py --sbml1=FILE --sbml2=FILE
    
    option:
        -h --help    Show help.
        --sbml1=FILE    path of the first sbml file
        --sbml2=FILE    path of the second sbml file
"""
from padmet.utils.exploration import compare_sbml
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml1_path = args["--sbml1"]
    sbml2_path = args["--sbml2"]
    compare_sbml.compare_sbml(sbml1_path, sbml2_path)

if __name__ == "__main__":
    main()