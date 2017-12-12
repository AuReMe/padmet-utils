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

usage:
    change_compart.py --padmet=FILE --old=STR --new=STR [--output=FILE]

options:
    -h --help     Show help.
    --padmet=FILE    pathname of the padmet file
    --old=STR    compartment id to remove
    --new=STR    new compartment id
    --output=FILE    new padmet pathname, if none, overwritting the original padmet
"""
from padmet.padmetSpec import PadmetSpec
import docopt

def main():
    #parsing args
    args = docopt.docopt(__doc__)
    padmet_file = args["--padmet"]
    old_compart = args["--old"]
    new_compart = args["--new"]
    new_padmet = args["--output"]
    if new_padmet is None:
        new_padmet = args["--padmet"]
    padmetSpec = PadmetSpec(padmet_file)

    #create the file with compounds ids and the current compartment corresponding to 'compart'    
    padmetSpec.change_compart(old_compart, new_compart)
    padmetSpec.generateFile(new_padmet)            

if __name__ == "__main__":
    main()