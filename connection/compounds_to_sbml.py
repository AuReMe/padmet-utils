# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:
Convert a list of compounds to an sbml file. Often in sbml, the compartment of the
compounds is concatenate to the id. If wanted, possible to add a compart each compounds with compart_name
If verbose, will also check if each compound is in a padmetRef and or a padmetSpec
In sbml: id = id encoded, name = original id
ex: id = cpd-a, id encoded = M_cpd__45__a_'value of compart_name'
if verbose True, will print if each compound was found in padmetRef and or padmetSpec

usage:
    compounds_to_sbml.py --compounds=FILE --output=FILE [--padmetRef=FILE] [--padmetSpec=FILE] [-v]

option:
    -h --help     Show help.
    --compounds=FILE    the pathname to the file containing the compounds ids and the compart, line = cpd-id\tcompart.
    --output=FILE    pathname to the sbml output.
    --padmetRef=FILE   pathname of the padmet used as reference (ex: metacyc_18.5.padmet)
    --padmetSpec=FILE   pathname of the padmet corresponding to the species studied
    -v    print info and check if ids in padmetRef and or padmet
"""
from padmet.sbmlGenerator import compounds_to_sbml
import docopt

def main():
    args = docopt.docopt(__doc__)
    compounds = args["--compounds"]     
    padmetRef = args["--padmetRef"]
    padmetSpec = args["--padmetSpec"]
    output = args["--output"]
    verbose = args["-v"]

    compounds_to_sbml(compounds, output, padmetRef, padmetSpec, verbose = verbose)

if __name__ == "__main__":
    main()