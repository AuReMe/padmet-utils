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
Create the files containing the wikicode for each page of the wiki relative to a
padmetSpec based on a database padmetRef. All the pages are stored in the folder
output_dir. The padmetRef is required to extract the information about all the reactions
of a given pathway. In the case of BIGG db, it's possible to use the padmetSpec as
padmetRef because there is no pathways.

usage:
    wikiPage_creation.py --padmetSpec=FILE --output_dir=DIR [--padmetRef=FILE] [-v]

option:
    -h --help    Show help.
    --padmetRef=FILE    pathname of the padmet of reference.
    --padmetSpec=FILE    pathname of the padmet of interest.
    --output_dir=DIR    pathname of the directory were to stock all the wikiPages.
    -v   print info
"""
import padmet.wikiGenerator as wg
from time import time
import docopt

def main():
    args = docopt.docopt(__doc__)
    
    padmetRef_file = args["--padmetRef"]
    padmetSpec_file = args["--padmetSpec"]
    output_dir = args["--output_dir"]
    verbose = args["-v"]

    chronoDepart = time()
    wg.create_all_wikiPages(padmetSpec_file, output_dir, padmetRef_file, verbose)
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print "done in: ", chrono, "s !"

if __name__ == "__main__":
    main()

