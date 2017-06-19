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
Convert a list of reactions to SBML file.

usage:
    reactions_to_sbml.py --padmetRef=FILE --reactions=FILE --output=FILE [-v]

option:
    -h --help     Show help.
    --padmetRef=FILE   pathname to the padmet used as reference (ex: metacyc.padmet)   
    --reactions=FILE    pathname to the file of reaction id. one id by line.
    --output=FILE    pathname to the sbml output.
    -v    print info
"""
from padmet.sbmlGenerator import reactions_to_SBML
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmetRef_file = args["--padmetRef"]
    reactions_file = args["--reactions"]     
    output = args["--output"]
    verbose = args["-v"]

    reactions_to_SBML(reactions_file, output, padmetRef_file, verbose=verbose)
            
if __name__ == "__main__":
    main()