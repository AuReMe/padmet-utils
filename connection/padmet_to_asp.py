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
Convert PADMet to ASP following this predicats:
direction(reaction_id, reaction_direction). #LEFT-TO-RIGHT or REVERSIBLE
ec_number(reaction_id, ec(x,x,x)).
catalysed_by(reaction_id, enzyme_id).
uniprotID(enzyme_id, uniprot_id). #if has has_xref and db = "UNIPROT"
common_name({reaction_id or enzyme_id or pathway_id or compound_id} , common_name)

usage:
    padmet_to_ASP.py --padmet=FILE --output=FILE [-v]

options:
    -h --help     Show help.
    --padmet=FILE    pathname of the padmet file to convert
    --output=FILE    pathname of the lp file
    -v   print info
"""
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmet_file = args["--padmet"]
    output = args["--output"]
    verbose = args["-v"]

    padmet_to_asp(padmet_file, output, verbose)

if __name__ == "__main__":
    main()