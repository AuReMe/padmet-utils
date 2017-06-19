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
@author: Nicolas Guillaudeux nicolas.guillaudeux@inria.fr

Description:

usage:
    convert_sbml_db.py --mnx_rxn=FILE --mnx_cpd=FILE --sbml=FILE --output_dict=FILE --output_sbml=FILE --db_in=ID [-v]

options:
    -h --help     Show help.
"""
from libsbml import *
#from padmet.padmetRef import PadmetRef
import re
import docopt

def main():
    args = docopt.docopt(__doc__)
    mnx_rxn_file = args["--mnx_rxn"]
    mnx_cpd_file = args["--mnx_cpd"]
    db_in = args["--db_in"]
    db_out = "metacyc"
    sbml_file = args["--sbml"]
    output_sbml = args["--output_sbml"]
    output_dict = args["--output_dict"]
    metacyc_file = args["--output_dict"]
    verbose = args["-v"]

    #metacyc = PadmetRef(metacyc_file)

    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print (document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    listOfSpecies = model.getListOfSpecies()

    if verbose:
        print("nb reactions: %s" %len(listOfReactions)) 

    mnx_rxn_dict = mnx_reader(mnx_rxn_file)
    mnx_rxn_dict_filtred = dict([(k,v) for k,v in mnx_rxn_dict.items() if db_in in v.keys() and db_out in v.keys()])

    if mnx_cpd_file:
        mnx_cpd_dict = mnx_reader(mnx_cpd_file)
        mnx_cpd_dict_filtred = dict([(k,v) for k,v in mnx_cpd_dict.items() if db_in in v.keys() and db_out in v.keys()])
    
    #k: reaction orignial id, v = list of reactions from db target
    mapped_rxn = {}
    mapped_species = {}
    #total_mapped_rxn = 0
    for sbml_rxn in listOfReactions:
        rxn_id = sbml_rxn.id
        uncoded = convert_from_coded_id(rxn_id)[0]
        if uncoded in ["APor","R00904"]:
            mapped_rxn[rxn_id] = "RXN-6382"
        elif uncoded == "HAL":
            mapped_rxn[rxn_id] = "CARNOSINE-SYNTHASE-RXN"
        else:
            match_dicts = [v for v in mnx_rxn_dict_filtred.values() if uncoded in v[db_in]]
            
            if match_dicts:
                if len(match_dicts) > 1:
                    print("more than one mnx for:%s" %uncoded)
                else:
                    match_dict = match_dicts[0]
                    if len(match_dict[db_out]) > 1:
                        #print("more than one match for:%s" %uncoded)
                        #print(";".join(match_dict[db_out]))
                        pass
                    else:
                        mapped_rxn[rxn_id] = match_dict[db_out][0]
                        #print("%s matched with %s" %(rxn_id, match_dict[db_out][0]))
        """
        if uncoded in dict_data.keys():
            if rxn_id in mapped_rxn.keys():
                if dict_data[uncoded] == mapped_rxn[rxn_id]:
                    print("%s same match with dict_data" %uncoded)
                else:
                    print("/!\ %s matched with %s in dict vs %s in mnx" %(uncoded, dict_data[uncoded], mapped_rxn[rxn_id]))
            else:
                print("%s matched with %s only with dict" %(uncoded, dict_data[uncoded]))
        elif rxn_id in mapped_rxn.keys():
            print("%s matched with %s only with mnx data" %(uncoded, mapped_rxn[rxn_id]))
        """ 
    all_mapped_rxn = set(mapped_rxn.keys())
    print ("%s/%s reactions mapped" %(len(mapped_rxn),len(listOfReactions)))

    nb_new_mapped_rxn = 0
    for sbml_rxn in [sbml_rxn for sbml_rxn in listOfReactions if sbml_rxn.id not in mapped_rxn]:
        all_species_id = set([s.getSpecies() for s in sbml_rxn.getListOfReactants()] + [s.getSpecies() for s in sbml_rxn.getListOfProducts()])
        for species_id in all_species_id:
            uncoded = convert_from_coded_id(species_id)[0]
        
        if uncoded in ["h2o","C00001"]:
            mapped_species[species_id] = "WATER"
        elif uncoded == "h2s":
            mapped_species[species_id] = "HS"
        elif uncoded == "no3":
            mapped_species[species_id] = "NITRATE"
        elif uncoded == "heme0":
            mapped_species[species_id] = "HEME_O"
        elif uncoded in ["akg","C00026"]:
            mapped_species[species_id] = "2-KETOGLUTARATE"
        elif uncoded in ["pi","C00009"]:
            mapped_species[species_id] = "Pi"
        elif uncoded == "for":
            mapped_species[species_id] = "FORMATE"
        elif uncoded == "trnagln":
            mapped_species[species_id] = "GLN-tRNAs"
        elif uncoded == "coa":
            mapped_species[species_id] = "CO-A"
        elif uncoded == "adn":
            mapped_species[species_id] = "ADENOSINE"
        elif uncoded == "ade":
            mapped_species[species_id] = "ADENINE"
        elif uncoded == "cytd":
            mapped_species[species_id] = "CYTIDINE"
        elif uncoded == "dtdp":
            mapped_species[species_id] = "TDP"
        elif uncoded == "C00008":
            mapped_species[species_id] = "ADP"
        elif uncoded in ["oaa","C00036"]:
            mapped_species[species_id] = "OXALACETIC_ACID"
        elif uncoded == "dhnpt":
            mapped_species[species_id] = "DIHYDRO-NEO-PTERIN"
        elif uncoded == "suc6p":
            mapped_species[species_id] = "SUCROSE-6P"
        elif uncoded == "g3pe":
            mapped_species[species_id] = "ALKYL-SN-GLYCERO-PHOSPHOETHANOLAMINE"
        elif uncoded == "retinol":
            mapped_species[species_id] = "Retinols"
        elif uncoded == "phpyr":
            mapped_species[species_id] = "PHENYL-PYRUVATE"
        elif uncoded == "C03912":
            mapped_species[species_id] = "L-DELTA1-PYRROLINE_5-CARBOXYLATE"
        elif uncoded == "C00014":
            mapped_species[species_id] = "AMMONIA"
        elif uncoded == "C00013":
            mapped_species[species_id] = "PPI"
        elif uncoded == "bamppald":
            mapped_species[species_id] = "CPD-6082"
        elif uncoded == "C00006":
            mapped_species[species_id] = "NADP"
        elif uncoded == "C00005":
            mapped_species[species_id] = "NADPH"
        elif uncoded in ["nadph","nadp","nadh","nad","hco3","ACP","hso3",
        "gly","ppi","amp","cdp","cmp","udp","ump","gdp","adp"]:
            mapped_species[species_id] = uncoded.upper()
        else:
            match_dicts = [v for v in mnx_cpd_dict_filtred.values() if uncoded in v[db_in] and db_out in v.keys()]

            if match_dicts:
                if len(match_dicts) > 1:
                    print("More than one mnx for:%s" %uncoded)
                else:
                    match_dict = match_dicts[0]
                    if len(match_dict[db_out]) > 1:
                        print("more than one match for:%s" %uncoded)
                    else:
                        mapped_species[species_id] = match_dict[db_out][0]
        if not all_species_id.difference(set(mapped_species.keys())):
            nb_new_mapped_rxn += 1
            all_mapped_rxn.add(sbml_rxn.id)
    print ("reactions recovered from compounds: %s" %nb_new_mapped_rxn)
    print ("total reactions mapped %s/%s" %(len(mapped_rxn)+nb_new_mapped_rxn,len(listOfReactions)))

    for rxn_id in set([rxn.id for rxn in listOfReactions]).difference(all_mapped_rxn):
        listOfReactions.remove(rxn_id)
    
    species_in_rxn_temp = [[r.getSpecies() for r in rxn.getListOfReactants()]+[p.getSpecies() for p in rxn.getListOfProducts()] for rxn in listOfReactions]
    species_in_rxn = set()
    [species_in_rxn.update(set(x)) for x in species_in_rxn_temp]
    [mapped_species.pop(sId) for sId in set([x.id for x in listOfSpecies]).difference(species_in_rxn) if sId in mapped_species.keys()]    
    [listOfSpecies.remove(sId) for sId in set([x.id for x in listOfSpecies]).difference(species_in_rxn)]    
       
    
    with open(output_dict, 'w') as f:
        for k,v in mapped_rxn.items():
            f.write(k+"\t"+v+"\n")
        for k,v in mapped_species.items():
            f.write(k+"\t"+v+"\n")
    
    
    writeSBMLToFile(document, output_sbml) 

def convert_from_coded_id(coded):
    """
    convert an id from sbml format to the original id. try to extract the type of
    the id and the compart using strong regular expression
    @param coded: the encoded id
    @type coded: str
    @return: (the uncoded id, type=None, compart=None)
    @rtype: tuple
    """
    #replace DASH from very old sbmls
    coded = coded.replace('_DASH_', '_')
    #an original id starting with int will start with '_' in sbml
    if coded.startswith("_"):
        coded = coded[1:]
    #reg ex to find the ascii used to replace not allowed char
    codepat = re.compile('__(\d+)__')
    #replac ascii by the not allowed char of sbml
    coded = codepat.sub(ascii_replace, coded)
    
    reg_expr = re.compile('(?P<_type>^[MR]_)(?P<_id>.*)(?P<compart>_.*)')
    search_result = reg_expr.search(coded)
    if search_result is not None:
        compart = search_result.group('compart').replace("_","")    
        _type = search_result.group('_type').replace("_","")
        uncoded = search_result.group('_id')
    else:
        reg_expr = re.compile('(?P<_type>^[MR]_)(?P<_id>.*)')
        search_result = reg_expr.search(coded)
        if search_result is not None:
            compart = None
            _type = search_result.group('_type').replace("_","")
            uncoded = search_result.group('_id')
        else:
            reg_expr = re.compile('(?P<_id>.*)(?P<compart>_.*)')
            search_result = reg_expr.search(coded)
            if search_result is not None:
                _type = None
                compart = search_result.group('compart').replace("_","")    
                uncoded = search_result.group('_id')
            else:
                uncoded = coded
                _type = None
                compart = None
            
    return (uncoded, _type, compart)

def ascii_replace(match):
    """
    recover banned char from the integer ordinal in the reg.match
    """
    return chr(int(match.group(1)))

def mnx_reader(input_file):
    with open(input_file, "r") as source:
        lecture = source.read().splitlines()
        list_data = [line.split("\t")[:2] for line in lecture if not line.startswith("#")]
    #print list_data[:10]
    mnx_dict = dict()
    all_mnxs = set([k for v,k in list_data])
    mnx_dict = dict([(k, dict()) for k in all_mnxs])
    #print a[:10]
    #print mnx_dict.items()[:10]
    for v,k in list_data:
        try:
            db, _id = v.split(":")
        except ValueError:
            db = v.split(":")[0]
            _id = ":".join(v.split(":")[1:])
        try:
            mnx_dict[k][db].append(_id)
        except KeyError:
            mnx_dict[k][db] = [_id]
    return mnx_dict



if __name__ == "__main__":
    main()

