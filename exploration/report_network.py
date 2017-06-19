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
Create reports of a padmet file. 
all_pathways.tsv: header = ["dbRef_id", "Common name", "Number of reaction found", 
            "Total number of reaction", "Ratio (Reaction found / Total)"]
all_reactions.tsv: header = ["dbRef_id", "Common name", "formula (with id)", 
            "formula (with common name)", "in pathways", "associated genes"]
all_metabolites.tsv: header = ["dbRef_id", "Common name", "Produced (p), Consumed (c), Both (cp)"]

usage:
    report_network.py --padmetSpec=FILE --output_dir=dir [--padmetRef=FILE] [-v]

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathname of the padmet file.
    --padmetRef=FILE    pathname of the padmet file used as database
    --output_dir=dir    directory for the results.
    -v   print info.
"""
from padmet.padmetSpec import PadmetSpec
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmetRef_file = args["--padmetRef"]
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    output_dir = args["--output_dir"]
    verbose = args["-v"]

    padmetSpec.network_report(output_dir, padmetRef_file, verbose)  

if __name__ == "__main__":
    main() 