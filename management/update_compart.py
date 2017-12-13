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
    update_compart.py --padmet=FILE
    update_compart.py --padmet=FILE --remove=STR [--output=FILE] [-v]
    update_compart.py --padmet=FILE --old=STR --new=STR [--output=FILE] [-v]

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
    to_remove = args["--remove"]
    verbose = args["-v"]
    if new_padmet is None:
        new_padmet = args["--padmet"]
    padmetSpec = PadmetSpec(padmet_file)

    if to_remove:
        if "," in to_remove:
            to_remove = to_remove.split(",")
        else:
            to_remove = [to_remove]
        for compart in to_remove:
            if verbose:
                print("removing reaction(s) in compartment %s" %compart)
            padmetSpec.delCompart(compart, verbose)
        padmetSpec.generateFile(new_padmet)            
    elif old_compart and new_compart:
        if verbose:
            print("Converting compartment %s to %s" %(old_compart, new_compart))
        padmetSpec.change_compart(old_compart, new_compart, verbose)
        padmetSpec.generateFile(new_padmet)            
    else:
        if verbose:
            print("List of compartments")
        print(list(padmetSpec.get_all_compart()))

if __name__ == "__main__":
    main()
