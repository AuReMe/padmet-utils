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
Use the script padmet_to_sbml to convert padmet to a sbml file
give only the id of the reaction to test (obj_coef)

usage:
    padmet_to_sbml.py --padmet=FILE --output=FILE [--obj_fct=ID] [--sbml_lvl=x] [--sbml_version=y] [-v]

option:
    -h --help    Show help.
    --padmet=FILE    pathname of the padmet file to convert into sbml
    --output=FILE    pathanme of the sbml file to generate.
    --obj_fct=ID    id of the reaction objective.
    --sbml_lvl=x    sbml level 2 is sufficient for FBA [default: 2].
    --sbml_version=y    sbml version 1 is sufficient for FBA [default: 1].
    -v   print info.
"""
from padmet.sbmlGenerator import padmet_to_sbml
from time import time
import docopt


def main():
    args = docopt.docopt(__doc__)
    padmet_file = args["--padmet"]
    output = args["--output"]
    sbml_lvl = int(args["--sbml_lvl"])
    sbml_version = int(args["--sbml_version"])
    obj_fct = args["--obj_fct"]

    verbose = args["-v"]
    chronoDepart = time()
        
    padmet_to_sbml(padmet_file, output, obj_fct, sbml_lvl, sbml_version, verbose)
    
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print "done in: ", chrono, "s !"

if __name__ == "__main__":
    main()