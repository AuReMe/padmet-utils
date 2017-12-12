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
    convert_sbml_db.py --mnx_rxn=FILE --mnx_cpd=FILE --sbml=FILE --output=FILE --db_out=ID [-v]

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
    db_out = args["--db_out"].upper()
    if db_out not in ["BIGG","METACYC","KEGG"]:
        raise ValueError('Please choose a database id in ["BIGG","METACYC","KEGG"]')
        exit()
    sbml_file = args["--sbml"]
    output_dict = args["--output"]
    verbose = args["-v"]

    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print (document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    listOfSpecies = model.getListOfSpecies()

    if verbose:
        print("nb reactions: %s" %len(listOfReactions)) 

    #For reactions: k = MNXid, v = {k=db_id,v=[list of ids]}
    mnx_rxn_dict = mnx_reader(mnx_rxn_file, db_out)

    #For species: k = MNXid, v = {k=db_id,v=[list of ids]} 
    mnx_cpd_dict = mnx_reader(mnx_cpd_file, db_out)

    #k: orignial id, v = ref id
    mapped_rxn = {}
    mapped_cpd = {}
    rxn_with_more_one_mapping = 0
    cpd_with_more_one_mapping = 0
    rxn_mapped_with_cpds = []
    for sbml_rxn in listOfReactions:
        rxn_id = sbml_rxn.id
        uncoded_rxn_id = convert_from_coded_id(rxn_id)[0]
        
        #first check in intern dict mapping
        match_id = intern_mapping(uncoded_rxn_id, db_out, "reaction")
        
        #check if in mnx_rxn_dict
        if match_id:
            mapped_rxn[rxn_id] = match_id
        else:
            for map_dict in mnx_rxn_dict.values():
                #print(mapp_dict)
                all_rxn_id = []
                [all_rxn_id.extend(i) for i in map_dict.values()]
                if uncoded_rxn_id in all_rxn_id:
                    matchs_rxns = map_dict[db_out]
                    if len(matchs_rxns) > 1: 
                        rxn_with_more_one_mapping += 1
                        if verbose:
                            print("More than one mapping for reaction %s:\t%s" %(rxn_id,matchs_rxns))
                    else:
                        mapped_rxn[rxn_id] = matchs_rxns[0]
                    break
    #for all non mapped rxn, check if able to map all speices
    for sbml_cpd in listOfSpecies:
        cpd_id = sbml_cpd.id
        uncoded_cpd_id = convert_from_coded_id(cpd_id)[0]

        match_id = intern_mapping(uncoded_cpd_id, db_out, "compound")

        #check if in mnx_rxn_dict
        if match_id:
            mapped_cpd[cpd_id] = match_id
        else:
            for map_dict in mnx_cpd_dict.values():
                #print(mapp_dict)
                all_cpd_id = []
                [all_cpd_id.extend(i) for i in map_dict.values()]
                if uncoded_cpd_id in all_cpd_id:
                    matchs_cpds = map_dict[db_out]
                    if len(matchs_cpds) > 1:
                        if db_out == "METACYC" and uncoded_cpd_id.upper() in matchs_cpds:
                            mapped_cpd[cpd_id] = uncoded_cpd_id.upper()
                        else: 
                            cpd_with_more_one_mapping += 1
                            if verbose:
                                print("More than one mapping for compound %s:\t%s" %(cpd_id,matchs_cpds))
                    else:
                        mapped_cpd[cpd_id] = matchs_cpds[0]
                    break

    for sbml_rxn in [i for i in listOfReactions if i.id not in mapped_rxn.keys()]:
        all_cpds = set([r.getSpecies() for r in sbml_rxn.getListOfReactants()] + [r.getSpecies() for r in sbml_rxn.getListOfProducts()])
        match_cpd_in_rxn = set([cpd_id for cpd_id in all_cpds if cpd_id in mapped_cpd.keys()])

        if len(match_cpd_in_rxn) == len(all_cpds):
            rxn_mapped_with_cpds.append(sbml_rxn.id)
    if verbose:
        print("#######")
        print("Mapped reactions: %s/%s" %(len(mapped_rxn.keys()),len(listOfReactions)))
        print("Reactions with more than one mapping: %s" %rxn_with_more_one_mapping)
        print("Mapped species: %s/%s" %(len(mapped_cpd.keys()),len(listOfSpecies)))            
        print("Species with more than one mapping: %s" %cpd_with_more_one_mapping)
        print("Mapped reactions from species: %s" %(len(rxn_mapped_with_cpds)))
        for i in rxn_mapped_with_cpds:
            print("\t%s" %i)
        print("Total reactions mapped:%s/%s" %(len(mapped_rxn.keys())+len(rxn_mapped_with_cpds),len(listOfReactions)))
        print("#######")

    with open(output_dict, 'w') as f:
        for k,v in mapped_rxn.items():
            f.write(k+"\t"+v+"\n")
        for k,v in mapped_cpd.items():
            f.write(k+"\t"+v+"\n")
    

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
    
    reg_expr = re.compile('(?P<_type>^[MRS]_)(?P<_id>.*)(?P<compart>_.*)')
    search_result = reg_expr.search(coded)
    if search_result is not None:
        compart = search_result.group('compart').replace("_","")    
        _type = search_result.group('_type').replace("_","")
        uncoded = search_result.group('_id')
    else:
        reg_expr = re.compile('(?P<_type>^[MRS]_)(?P<_id>.*)')
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

def mnx_reader(input_file, db_out):
    with open(input_file, "r") as f:
        dataInArray = [line.split("\t")[:2] for line in f.read().splitlines() if not line.startswith("#")]

    #print list_data[:10]
    mnx_dict = dict()
    all_mnxs = set([k for v,k in dataInArray])
    mnx_dict = dict([(k, dict()) for k in all_mnxs])
    #print a[:10]
    #print mnx_dict.items()[:10]
    for v,k in dataInArray:
        try:
            db, _id = v.split(":")
        except ValueError:
            db = v.split(":")[0]
            _id = ":".join(v.split(":")[1:])
        db = db.upper()
        try:
            mnx_dict[k][db].append(_id)
        except KeyError:
            mnx_dict[k][db] = [_id]
    #clean ids not in db_out
    for k,v in mnx_dict.items():
        if db_out not in v.keys() or len(v.keys()) == 1:
            mnx_dict.pop(k)
    """
    for v in mnx_dict.values():
        try:
            rxn_db_out = v[db_out]
            all_mapp_rxn = []
            [all_mapp_rxn.extend(i) for i in v.values()]
    """
    return mnx_dict


def intern_mapping(id_to_map, db_out, _type):
    
    intern_rxn_dict = {
    "UNIQ_ID_1":{"METACYC":["RXN-6382"],"KEGG":["R00904"],"UNKNOWN":["APor"]},
    "UNIQ_ID_2":{"METACYC":["CARNOSINE-SYNTHASE-RXN"],"UNKNOWN":["HAL"]},
    "UNIQ_ID_3":{"METACYC":["ALANINE-AMINOTRANSFERASE-RXN"],"KEGG":["R00258"]},
    "UNIQ_ID_4":{"METACYC":["ASPAMINOTRANS-RXN"],"KEGG":["R00355"],"BIGG":["ASPATh"]},
    "UNIQ_ID_5":{"METACYC":["PHEAMINOTRANS-RXN"],"KEGG":["R00694"],"BIGG":["POAT"]},
    "UNIQ_ID_6":{"METACYC":["ORNCARBAMTRANSFER-RXN"],"KEGG":["R01398"],"BIGG":["OCT"]},
    "UNIQ_ID_7":{"METACYC":["3PGAREARR-RXN"],"KEGG":["R01518"],"BIGG":["PGM"]},
    "UNIQ_ID_8":{"METACYC":["6PGLUCONDEHYDROG-RXN"],"KEGG":["R01528"],"BIGG":["PGDHh"]},
    "UNIQ_ID_9":{"METACYC":["4-HYDROXY-2-KETOPIMELATE-LYSIS-RXN"],"KEGG":["R01645"]},
    "UNIQ_ID_10":{"METACYC":["AMINEPHEN-RXN"],"KEGG":["R02613"]},
    "UNIQ_ID_11":{"METACYC":["SHIKIMATE-5-DEHYDROGENASE-RXN"],"KEGG":["R02413"]},
    "UNIQ_ID_12":{"METACYC":["CARBODEHYDRAT-RXN"],"BIGG":["HCO3E"]}
    }
    
    intern_cpd_dict = {
    "UNIQ_ID_1":{"METACYC":["WATER"],"KEGG":["C00001"],"BIGG":["h2o"],"UNKNOWN":["H2O"]},
    "UNIQ_ID_2":{"METACYC":["2-KETOGLUTARATE"],"BIGG":["akg"],"KEGG":["C00026"]},
    "UNIQ_ID_3":{"METACYC":["Pi"],"BIGG":["pi"],"KEGG":["C00009"]},
    "UNIQ_ID_4":{"METACYC":["OXALACETIC_ACID"],"BIGG":["oaa"],"KEGG":["C00036"]},

    "UNIQ_ID_5":{"METACYC":["HS"],"BIGG":["h2s"]},
    "UNIQ_ID_6":{"METACYC":["NITRATE"],"BIGG":["no3"]},
    "UNIQ_ID_7":{"METACYC":["HEME_O"],"BIGG":["heme0"]},
    "UNIQ_ID_8":{"METACYC":["FORMATE"],"BIGG":["for"]},
    "UNIQ_ID_9":{"METACYC":["GLN-tRNAs"],"BIGG":["trnagln"]},
    "UNIQ_ID_10":{"METACYC":["CO-A"],"BIGG":["coa"]},
    "UNIQ_ID_11":{"METACYC":["ADENOSINE"],"BIGG":["adn"]},
    "UNIQ_ID_12":{"METACYC":["ADENINE"],"BIGG":["ade"]},
    "UNIQ_ID_13":{"METACYC":["CYTIDINE"],"BIGG":["cytd"]},
    "UNIQ_ID_14":{"METACYC":["CYTIDINE"],"BIGG":["dtdp"]},
    "UNIQ_ID_15":{"METACYC":["DIHYDRO-NEO-PTERIN"],"BIGG":["dhnpt"]},
    "UNIQ_ID_16":{"METACYC":["SUCROSE-6P"],"BIGG":["suc6p"]},
    "UNIQ_ID_17":{"METACYC":["ALKYL-SN-GLYCERO-PHOSPHOETHANOLAMINE"],"BIGG":["g3pe"]},
    "UNIQ_ID_18":{"METACYC":["Retinols"],"BIGG":["retinol"]},
    "UNIQ_ID_19":{"METACYC":["PHENYL-PYRUVATE"],"BIGG":["phpyr"]},
    "UNIQ_ID_20":{"METACYC":["CPD-6082"],"BIGG":["bamppald"]},

    "UNIQ_ID_21":{"METACYC":["ADP"],"KEGG":["C00008"]},
    "UNIQ_ID_22":{"METACYC":["L-DELTA1-PYRROLINE_5-CARBOXYLATE"],"KEGG":["C03912"]},
    "UNIQ_ID_23":{"METACYC":["AMMONIA"],"KEGG":["C00014"]},
    "UNIQ_ID_24":{"METACYC":["PPI"],"KEGG":["C00013"]},
    "UNIQ_ID_25":{"METACYC":["NADP"],"KEGG":["C00006"]},
    "UNIQ_ID_26":{"METACYC":["NADPH"],"KEGG":["C00005"]},

    "UNIQ_ID_27":{"METACYC":["GLY"],"UNKNOWN":["Glycine"]},
    "UNIQ_ID_28":{"METACYC":["METOH"],"UNKNOWN":["Methanol"]},
    "UNIQ_ID_29":{"METACYC":["PPI"],"UNKNOWN":["Pyrophosphate"]},
    }


    if _type == "reaction":
        for mapp_dict in intern_rxn_dict.values():
            all_rxn_id = []
            [all_rxn_id.extend(i) for i in mapp_dict.values()]
            if id_to_map in all_rxn_id:
                return mapp_dict[db_out][0]

    elif _type == "compound":
        for mapp_dict in intern_cpd_dict.values():
            all_cpd_id = []
            [all_cpd_id.extend(i) for i in mapp_dict.values()]
            if id_to_map in all_cpd_id:
                return mapp_dict[db_out][0]

        if db_out == "METACYC":
            for mapp_dict in intern_cpd_dict.values():
                all_cpd_id = []
                [all_cpd_id.extend(i) for i in mapp_dict.values()]
                if id_to_map.upper() in all_cpd_id:
                    return mapp_dict[db_out][0]


    return None


if __name__ == "__main__":
    main()

