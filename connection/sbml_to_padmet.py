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
There are 3 cases of convertion sbml to padmet:

1./ Creation of a reference database in padmet format from sbml(s) (or updating one with new(s) sbml(s))
First usage, padmetRef is the padmetRef to create or to update. If it's an update case, the output
can be used to create a new padmet, if output None, will overwritte the input padmetRef.

2./ Creation of a padmet representing an organism in padmet format from sbml(s) (or updating one with new(s) sbml(s))
2.A/ Without a database of reference:
Second usage, padmetSpec is the padmetSpec to create or update. If it's an update case, the output
can be used to create a new padmet, if output None, will overwritte the input padmetSpec.

2.B/ With a database of refence:
Third usage, padmetSpec is the padmetSpec to create or update. If it's an update case, the output
can be used to create a new padmet, if output None, will overwritte the input padmetSpec.
padmetRef is the padmet representing the database of reference.

It is possible to define a specific policy and info for the padmet. To learn more about
policy and info check doc of lib.padmetRef/Spec.
if the ids of reactions/compounds are not the same between padmetRef and the sbml, it is possible to use
a dictionnary of association (sbml_id padmetRef_id)
with one line = 'id_sbml \t id_padmetRef'
Finally if a reaction from sbml is not in padmetRef, it is possible to force the copy and creating
a new reaction in padmetSpec with the arg -f

usage:
    sbml_to_padmet.py --sbml=FILE --padmetRef=FILE [--output=FILE] [--policy=LIST] [--info=DICT] [-v]
    sbml_to_padmet.py --sbml=FILE --padmetSpec=FILE [--output=FILE] [--policy=LIST] [--info=DICT] [-v]
    sbml_to_padmet.py --sbml=FILE --padmetRef=FILE --padmetSpec=FILE [--output=FILE] [--assocIdOriginRef=FILE] [-v] [-f]

options:
    -h --help     Show help.
    --padmetSpec=FILE    pathname to the padmet file to update with the sbml. If there's no padmetSpec, just specify the output
    --padmetRef=FILE    pathanme to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
    --sbml=FILE    1 sbml file to convert into padmetSpec (ex: my_network.xml/sbml) OR a directory with n SBML
    --output=FILE   pathanme to the new padmet file
    --assocIdOriginRef=FILE    dictionnary of association id_origin id_ref
    --policy=LIST    
    --info=DICT    
    -f   allows to create reactions not found in the padmetRef (force)
    -v   print info
"""
from padmet.padmetSpec import PadmetSpec
from padmet.padmetRef import PadmetRef
from time import time
from datetime import datetime
import os
import docopt

def main():
    default_policy = [['compound','has_name','name'], ['compound','has_xref','xref'], ['compound','has_suppData','suppData'],
                    ['gene','has_name','name'], ['gene','has_xref','xref'], ['gene','has_suppData','suppData'], ['gene','codes_for','protein'],
                    ['pathway','has_name','name'], ['pathway','has_xref','xref'], ['pathway','is_in_pathway','pathway'], 
                    ['protein','has_name','name'], ['protein','has_xref','xref'], ['protein','has_suppData','suppData'], ['protein','catalyses','reaction','ASSIGNMENT','X','SCORE','Y'],
                    ['reaction','has_name','name'], ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','is_in_pathway','pathway'],  
                    ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                    ['reaction','consumes','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], 
                    ['reaction','is_linked_to','gene','ASSIGNMENT','X']]
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")
    default_info = {"PADMET":{"creation":today_date,"version":"2.3"},"DB_info":{"DB":"NA","version":"NA"},"Model_info":{"Description":"NA"}}
    #recovering args
    args = docopt.docopt(__doc__)
    verbose = args["-v"]
    if args["--info"] is None:
        info = default_info
    else:
        info = dict(args["--info"])
    if args["--policy"] is None:
        policy = default_policy
    else:
        policy = list(args["--policy"])
    

    #if sbml is a directory, recover all file path in a list. if no => only one file: create a list with only this file
    sbmlFileName = []
    if os.path.isdir(args["--sbml"]):
        for abspath,y,filenames in os.walk(args["--sbml"]):
            sbmlFileName += [abspath+"/"+filename for filename in filenames if not filename.startswith(".") and filename.endswith((".sbml",".xml"))]
    else:
        sbmlFileName = [args["--sbml"]]
    if not sbmlFileName:
        if verbose: print("No sbml found in %s" %(args["--sbml"]))
        exit()

    if args["--padmetRef"] is not None and args["--padmetSpec"] is None:
        if args["--output"] is not None:
            output = args["--output"]
        else:
            output = args["--padmetRef"]
        if os.path.isfile(args["--padmetRef"]):
            padmetRef = PadmetRef(args["--padmetRef"])
        else:
            padmetRef = PadmetRef()
            padmetRef.setInfo(info)
            padmetRef.setPolicy(policy)
    
        chronoDepart = time()
        #CORE
        for sbml_file in sbmlFileName:
            padmetRef.initFromSbml(sbml_file, verbose)
        padmetRef.generateFile(output)
        if verbose: print("Generated file: %s" %output)
       
        chrono = (time() - chronoDepart)
        partie_entiere, partie_decimale = str(chrono).split('.')
        chrono = ".".join([partie_entiere, partie_decimale[:3]])
        if verbose: print "done in: ", chrono, "s !"
    
    elif args["--padmetSpec"] is not None and args["--padmetRef"] is None:
        if args["--output"] is not None:
            output = args["--output"]
        else:
            output = args["--padmetSpec"]
        if os.path.isfile(args["--padmetSpec"]):
            padmetSpec = PadmetSpec(args["--padmetSpec"])
        else:
            padmetSpec = PadmetSpec()
            padmetSpec.setInfo(info)
            padmetSpec.setPolicy(policy)

        chronoDepart = time()
        #CORE
        for sbml_file in sbmlFileName:
            padmetSpec.updateFromSbml(sbml_file, None, None, verbose)
        padmetSpec.generateFile(output)
        if verbose: print("Generated file: %s" %output)
        
        chrono = (time() - chronoDepart)
        partie_entiere, partie_decimale = str(chrono).split('.')
        chrono = ".".join([partie_entiere, partie_decimale[:3]])
        if verbose: print "done in: ", chrono, "s !"
    
    else:
        padmetRef = PadmetRef(args["--padmetRef"])
        assocIdOriginRef = args["--assocIdOriginRef"]
        force = args["-f"]
        if args["--output"] is not None:
            output = args["--output"]
        else:
            output = args["--padmetSpec"]
        if os.path.isfile(args["--padmetSpec"]):
            padmetSpec = PadmetSpec(args["--padmetSpec"])
        else:
            padmetSpec = PadmetSpec()
            padmetSpec.setInfo(padmetRef)
            padmetSpec.setPolicy(padmetRef)

        chronoDepart = time()
        #CORE
        for sbml_file in sbmlFileName:
            padmetSpec.updateFromSbml(sbml_file, padmetRef, assocIdOriginRef, verbose, force)
        padmetSpec.generateFile(output)
        if verbose: print("Generated file: %s" %output)
        
        chrono = (time() - chronoDepart)
        partie_entiere, partie_decimale = str(chrono).split('.')
        chrono = ".".join([partie_entiere, partie_decimale[:3]])
        if verbose: print "done in: ", chrono, "s !"


if __name__ == "__main__":
    main()
