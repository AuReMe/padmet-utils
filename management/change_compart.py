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
Two disctinct usages:
    change_compart.py --padmetSpec=FILE --compart=STR --output=FILE [-v]
For a given compartment to check return a file with all metabolites in it
Then, to change the compartment of one or more metabolites, modify the precedent output
by changing the compartment.
ex: output line: M_icit_x	x
M_icit_x is in compart x. to change the compart to 'c' for example, edit 'x' to 'c'.
Finally use the second usages to update the model.
    change_compart.py --padmetSpec=FILE --updateFile=FILE [--new_padmet=FILE] [-v]
NB: the metabolite id also change: M_icit_x => M_icit_c

usage:
    change_compart.py --padmetSpec=FILE --compart=STR --output=FILE [-v]
    change_compart.py --padmetSpec=FILE --updateFile=FILE [--new_padmet=FILE] [-v]

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathname of the padmet file to update from updateFile
    --compart=STR    the compartment id to check
    --output=FILE    pathname of the tsv file with col 1 = compounds ids, col 2 = the current compartment
    --updateFile=FILE    pathname of the output after modifications of col 2
    --new_padmet=FILE   pathname of the new padmet. if None, will erase the current
    -v   print info
"""
from padmet.padmetSpec import PadmetSpec
import docopt

def main():
    #parsing args
    args = docopt.docopt(__doc__)
    padmetSpec_file = args["--padmetSpec"]
    compart = args["--compart"]
    updateFile = args["--updateFile"]
    output = args["--output"]
    verbose = args["-v"]
    new_padmet = args["--new_padmet"]
    if new_padmet is None:
        new_padmet = args["--padmetSpec"]
    padmetSpec = PadmetSpec(padmetSpec_file)

    #create the file with compounds ids and the current compartment corresponding to 'compart'    
    if updateFile is None:
        compounds = set([rlt.id_out for rlt in padmetSpec.getAllRelation()
        if rlt.type in ["consumes","produces"] and rlt.misc["COMPARTMENT"][0] == compart])
        fileInArray = [cpd+"\t"+compart for cpd in compounds]
        with open(output, 'w') as f:
            f.write("\n".join(fileInArray))        
        #if updateFile is not None, use it to update compartments
    else:
        with open(updateFile,'r') as f:
            data = [line.split("\t") for line in f.read().splitlines()]
        for compound_id, compart in data:
            for rlt in padmetSpec.getAllRelation():
                if rlt.type in ["consumes","produces"] and rlt.id_out == compound_id:
                    rlt.misc.update({"COMPARTMENT":compart})
        padmetSpec.generateFile(new_padmet)            

if __name__ == "__main__":
    main()