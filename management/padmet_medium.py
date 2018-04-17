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
For a given set of compounds representing the growth medium (or seeds). Create 2 reactions
for each compounds to maintain consistency of the network for flux analysis.
For each compounds create 
    -an exchange reaction: this reaction consumes the compound in the
compartment 'C-BOUNDARY' and produces the compound in the compartment 'e' extracellular
    -a transport reaction: this reaction consumes the compound in the compartment
'e' extracellular' and produces the compound in the compartment 'c' cytosol
ex: for seed 'cpd-a'
1/ check if cpd-a in padmetSpec, if not, copy from padmetRef.
2/ create exchange reaction: ExchangeSeed_cpd-a_b: 1 cpd-a (C-BOUNDARAY) <=> 1 cpd-a (e)
3/ create transport reaction: TransportSeed_cpd-a_e: 1 cpd-a (e) => 1 cpd-a (c)
4/ create a new file if output not None, or overwritte padmetSpec

usage:
    padmet_medium.py --padmetSpec=FILE
    padmet_medium.py --padmetSpec=FILE -r [--output=FILE] [-v]
    padmet_medium.py --padmetSpec=FILE --padmetRef=FILE --seeds=FILE [--output=FILE] [--b_compart=STR] [--e_compart=STR] [--c_compart=STR] [-v]

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathanme to the padmet file to update
    --padmetRef=FILE    pathname to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
    --seeds=FILE    the pathname to the file containing the compounds ids and the compart, line = cpd-id\tcompart.
    --output=FILE    If not None, pathname to the padmet file updated
    -v   print info

"""
from padmet.classes import PadmetSpec
from padmet.classes import PadmetRef
from padmet.classes import Node
from padmet.classes import Relation
import docopt

def main():
    args = docopt.docopt(__doc__)
    if args["--seeds"]:
        with open(args["--seeds"], 'r') as f:
            seeds = [line.split("\t")[0] for line in f.read().splitlines()]
    else:
        seeds = None
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    padmetRef  = PadmetRef(args["--padmetRef"])
    output = args["--output"]
    verbose = args["-v"]
    #b_compart = args["--b_compart"]
    #e_compart = args["--e_compart"]
    #c_compart = args["--c_compart"]
    remove = args["-r"]
    if output is None:
        output = args["--padmetSpec"]
    
    if not remove and not seeds:
        g_m = padmetSpec.get_growth_medium()
        print("List of growth medium:")
        if g_m:
            print(list(g_m))
        else:
            print("[]")

    elif seeds:
        padmetSpec.set_growth_medium(new_growth_medium = seeds, padmetRef = padmetRef, verbose = verbose)
        padmetSpec.generateFile(output)
    
    elif remove:
        padmetSpec.remove_growth_medium(verbose = verbose)
        padmetSpec = PadmetSpec(output)

if __name__ == "__main__":
    main()  

                
    
    