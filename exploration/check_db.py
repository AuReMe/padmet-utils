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
For a given sbml file, return the DBs of each reaction, based on metanetx file.

usage:
    check_db.py --sbml=FILE --mnx=FILE

options:
    -h --help     Show help.
    --sbml=FILE     Path to the sbml file to be analysed.
    --mnx=FILE     path to the metatnetx file for reaction, reac_xref.tsv.
"""
from libsbml import *
from padmet.sbmlPlugin import convert_from_coded_id
import docopt

def main():
    args = docopt.docopt(__doc__)
    rxn_mnx = args["--mnx"]
    sbml_file = args["--sbml"]
    
    with open(rxn_mnx,'r') as f:
        db = dict([tuple(line.split("\t")[0].split(":")[::-1]) for line in f.read().splitlines() if not line.startswith("#")])

    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print (document.getError(i).getMessage())
    
    model = document.getModel()
    listOfSpecies = model.getListOfSpecies()
    listOfReactions = model.getListOfReactions()
    nbReactions = len(listOfReactions)
    nbSpecies = len(listOfSpecies)
    print("nb species: %s" %nbSpecies)
    print("nb reactions: %s" %nbReactions) 
    
    data = [convert_from_coded_id(r.id)[0] for r in listOfReactions]
    
    result_dict = {"None":[]}
    for rxn in data:
        if rxn in db.keys():
            try:
                result_dict[db[rxn]].append(rxn)
            except KeyError:
                result_dict[db[rxn]] = [rxn]
        else:
            result_dict["None"].append(rxn)
    for k,v in result_dict.iteritems():
        print("%s: %s reactions" %(k,len(v)))
        for r in v:
            print("\t%s" %r)

if __name__ == "__main__":
    main()
