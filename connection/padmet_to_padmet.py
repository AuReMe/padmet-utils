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
    padmet_to_padmet.py --init_padmet=FILE --padmetRef=FILE --to_add=FILE [--output=FILE] [-v]

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathname to the padmet file to update with the sbml. If there's no padmetSpec, just specify the output
    --padmetRef=FILE    pathanme to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
    --sbml=FILE    1 sbml file to convert into padmetSpec (ex: my_network.xml/sbml) OR a directory with n SBML
    --output=FILE   pathanme to the new padmet file
    -v   print info
"""
from padmet.padmetSpec import PadmetSpec
from padmet.padmetRef import PadmetRef
from time import time
from datetime import datetime
import os
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmetRef_file = args["--padmetRef"]
    output = args["--output"]
    verbose = args["-v"]
    padmetSpec_file = args["--init_padmet"]

    if os.path.isdir(args["--to_add"]):
        padmet_type = "path"
    elif os.path.isfile(args["--to_add"]):
        padmet_type = "file"
    else:
        raise TypeError("%s is not a dir or a file" %(args["--to_add"]))

    if padmetRef_file and os.path.isfile(padmetRef_file):
        padmetRef = PadmetRef(padmetRef_file)
    else:
        padmetRef = None

    if os.path.isfile(padmetSpec_file):
        padmetSpec = PadmetSpec(padmetSpec_file)
    else:
        now = datetime.now()
        today_date = now.strftime("%Y-%m-%d")

        padmetSpec = PadmetSpec()
        if padmetRef:
            padmetSpec.setInfo(padmetRef)
            padmetSpec.info["PADMET"]["creation"] = today_date
            padmetSpec.setPolicy(padmetRef)
        else:
            POLICY_IN_ARRAY = [['class','is_a_class','class'], ['class','has_name','name'], ['class','has_xref','xref'], ['class','has_suppData','suppData'],
                            ['compound','is_a_class','class'], ['compound','has_name','name'], ['compound','has_xref','xref'], ['compound','has_suppData','suppData'],
                            ['gene','is_a_class','class'], ['gene','has_name','name'], ['gene','has_xref','xref'], ['gene','has_suppData','suppData'], ['gene','codes_for','protein'],
                            ['pathway','is_a_class','class'], ['pathway','has_name','name'], ['pathway','has_xref','xref'], ['pathway','is_in_pathway','pathway'], 
                            ['protein','is_a_class','class'], ['protein','has_name','name'], ['protein','has_xref','xref'], ['protein','has_suppData','suppData'], ['protein','catalyses','reaction'],
                            ['protein','is_in_species','class'], 
                            ['reaction','is_a_class','class'], ['reaction','has_name','name'], ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','has_reconstructionData','reconstructionData'], ['reaction','is_in_pathway','pathway'],  
                            ['reaction','consumes','class','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','class','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                            ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                            ['reaction','consumes','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                            ['reaction','is_linked_to','gene','SOURCE:ASSIGNMENT','X:Y']]
            dbNotes = {"PADMET":{"creation":today_date,"version":"2.6"},"DB_info":{"DB":db,"version":version}}
            padmetSpec.setInfo(dbNotes)
            padmetSpec.setPolicy(POLICY_IN_ARRAY)

    #if sbml is a directory, recover all file path in a list. if no => only one file: create a list with only this file
    if padmet_type == "path":
        path = args["--to_add"]
        if not path.endswith("/"):
            path += "/"
        all_files = [i for i in os.walk(path).next()[2] if not i.startswith(".~lock")]
        padmetFiles = [path+i for i in all_files if i.endswith(".padmet")]
    else:
        padmetFiles = [args["--to_add"]]

    chronoDepart = time()
    #CORE
    for padmet_file in padmetFiles:
        if verbose:
            print("Updating %s from %s" %(os.path.basename(padmetSpec_file),os.path.basename(padmet_file)))
        padmet = PadmetSpec(padmet_file)
        padmetSpec.updateFromPadmet(padmet)

    if len(padmetFiles) == 0:
        if verbose: print("No padmet found in %s" %args["--to_add"])
    else:
        if not output:
            output = padmetSpec_file
        if verbose: print("Generated file: %s" %output)
        padmetSpec.generateFile(output)
        
    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print "done in: ", chrono, "s !"


if __name__ == "__main__":
    main()
