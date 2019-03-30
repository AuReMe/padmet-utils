# -*- coding: utf-8 -*-
"""
Description:
    For a given sbml using a specific database. Return a dictionnary of mapping from this database to a choosen
    the output is a file with line = reaction_id in origin database, reaction_id in db_out database
    if a reaction can not be mapped, try to map the compounds and if all the compounds are mapped. Insert in
    the output file the mapping of those compounds

::
    
    usage:
        convert_sbml_db.py --mnx_reac=FILE --mnx_chem=FILE --sbml=FILE --to-map=STR [-v]
        convert_sbml_db.py --mnx_folder=DIR --sbml=FILE --to-map=STR [-v]
        convert_sbml_db.py --mnx_folder=DIR --sbml=FILE --output=FILE --db_out=ID --to-map=STR [-v]
        convert_sbml_db.py --mnx_reac=FILE --mnx_chem=FILE --sbml=FILE --output=FILE --db_out=ID --to-map=STR [-v]
    
    options:
        -h --help     Show help.
        --mnx_reac=FILE     path to the MetaNetX file for reactions
        --mnx_chem=FILE     path to the MetaNetX file for compounds
        --sbml=FILE     path to the sbml file to convert
        --output=FILE     path to the file containing the mapping, sep = "\t"
        --db_out=FILE     id of the output database in ["BIGG","METACYC","KEGG"]
        -v     verbose.

"""
import libsbml
from padmet.utils.sbmlPlugin import get_all_decoded_version
import docopt
import os

def main():
    args = docopt.docopt(__doc__)
    verbose = args["-v"]
    to_map = args["--to-map"]
    if args["--mnx_folder"]:
        mnx_folder = args["--mnx_folder"]
        mnx_reac_file = os.path.join(mnx_folder, "reac_xref.tsv")
        mnx_chem_file = os.path.join(mnx_folder, "chem_xref.tsv")
    else:
        mnx_reac_file = args["--mnx_reac"]
        mnx_chem_file = args["--mnx_chem"]
    sbml_file = args["--sbml"]

    if args["--db_out"]:
        db_out = args["--db_out"].upper()
        output_dict = args["--output"]
        map_sbml(sbml_file, to_map, db_out, output_dict, verbose, mnx_reac_file, mnx_chem_file)
    else:
        db_select, db_found = check_sbml_db(sbml_file, to_map, verbose, mnx_reac_file, mnx_chem_file)
        print(db_select, db_found)

def check_sbml_db(sbml_file, to_map, verbose = False, mnx_reac_file = None, mnx_chem_file = None):
    """
    
    """
    if to_map not in ["all", "reaction", "species"]:
        raise ValueError("%s must be in [all, reaction, species]" %to_map)

    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    listOfSpecies = model.getListOfSpecies()

    unknown_db = "Unknown"
    db_found = {unknown_db: 0}
    if to_map == "all":
        db_found["total_reaction_species"] = 0
    if verbose:
        print("Check from which database is this sbml:")
    if to_map in ["all", "reaction"]:
        db_found['total_reaction'] = 0
        with open(mnx_reac_file, "r") as f:
            dict_reaction_id_db = dict([(line.split("\t")[0].split(":")[1], line.split("\t")[0].split(":")[0])  for line in f.read().splitlines()
            if not line.startswith("#") and ":" in line.split("\t")[0]])

        for reaction in listOfReactions:
            reaction_id = reaction.id
            db_found['total_reaction'] += 1
            if to_map == "all":
                db_found["total_reaction_species"] += 1
            all_reaction_id_decoded = get_all_decoded_version(reaction_id, "reaction")
   
            for reaction_id_decoded in all_reaction_id_decoded:
                db_match = dict_reaction_id_db.get(reaction_id_decoded, unknown_db)
                if db_match != unknown_db:
                    break
            if verbose:
                if db_match != unknown_db:
                    print("reaction: original id: %s; decoded id:%s; db_match:%s" %(reaction_id, reaction_id_decoded, db_match))
                else:
                    print("reaction: original id:%s; all decoded id:%s; no db found" %(reaction_id, all_reaction_id_decoded))
            try:
                db_found[db_match] += 1
            except KeyError:
                db_found[db_match] = 1

    if to_map in ["all", "species"]:
        db_found["total_species"] = 0
        with open(mnx_chem_file, "r") as f:
            dict_species_id_db = dict([(line.split("\t")[0].split(":")[1], line.split("\t")[0].split(":")[0])  for line in f.read().splitlines()
            if not line.startswith("#") and ":" in line.split("\t")[0]])

        for species in listOfSpecies:
            species_id = species.id
            db_found['total_species'] += 1
            if to_map == "all":
                db_found["total_reaction_species"] += 1

            all_species_id_decoded = get_all_decoded_version(species_id, "species")
    
            for species_id_decoded in all_species_id_decoded:
                db_match = dict_species_id_db.get(species_id_decoded, unknown_db)
                if db_match != unknown_db:
                    break
            if verbose:
                if db_match != unknown_db:
                    print("species: original id: %s; decoded id:%s; db_match:%s" %(species_id, species_id_decoded, db_match))
                else:
                    print("species: original id:%s; all decoded id:%s; no db found" %(species_id, all_species_id_decoded))
            try:
                db_found[db_match] += 1
            except KeyError:
                db_found[db_match] = 1

    if to_map == "all":
        db_select = [k for k, v in list(db_found.items())
                     if v == max([j for i,j in list(db_found.items()) if i != 'total_reaction_species'])][0]
    elif to_map == "reaction":
        db_select = [k for k, v in list(db_found.items())
                     if v == max([j for i,j in list(db_found.items()) if i != 'total_reaction'])][0]
    elif to_map == "species":
        db_select = [k for k, v in list(db_found.items())
                     if v == max([j for i,j in list(db_found.items()) if i != 'total_species'])][0]

    return (db_select, db_found)

def map_sbml(sbml_file, to_map, db_out, output_dict, verbose = False, mnx_reac_file = None, mnx_chem_file = None):
    """
    """
    map_from_cpd = False
    if to_map not in ["all", "reaction", "species"]:
        raise ValueError("%s must be in [all, reaction, species]" %to_map)

    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    listOfSpecies = model.getListOfSpecies()
                
    if db_out not in ["BIGG","METACYC","KEGG"]:
        raise ValueError('Please choose a database id in ["BIGG","METACYC","KEGG"]')


    #k: orignial id, v = ref id
    dict_sbml_id_mapped_id = {}
    if to_map in ["all", "reaction"]:
        #For reactions: k = MNXid, v = {k=db_id,v=[list of ids]}
        mnx_reac_dict = mnx_reader(mnx_reac_file, db_out)
        reaction_with_more_one_mapping = 0
        count_nb_stric_reac_mapped = 0

        for reaction in listOfReactions:
            reaction_id = reaction.id
            match_ids = None
            all_reaction_id_decoded = get_all_decoded_version(reaction_id, "reaction")
            for reaction_id_decoded in all_reaction_id_decoded:
                if not match_ids:
                    #first check in intern dict mapping
                    match_ids = intern_mapping(reaction_id_decoded, db_out, "reaction")
                    if match_ids:
                        dict_sbml_id_mapped_id[reaction_id] = match_ids
                        count_nb_stric_reac_mapped += 1
                        if verbose:
                            print("reaction: original id:%s; decoded id:%s; match_with:%s from intern mapping" %(reaction_id, reaction_id_decoded, match_ids))
                    #check if in mnx_reac_dict
                    else:
                        match_ids = get_from_mnx(mnx_reac_dict, reaction_id_decoded, db_out)
                        if match_ids:
                            if len(match_ids) > 1: 
                                reaction_with_more_one_mapping += 1
                                if verbose:
                                    print("reaction: original id:%s; decoded id:%s; More than one mapping:%s" %(reaction_id, reaction_id_decoded, match_ids))
                            else:
                                dict_sbml_id_mapped_id[reaction_id] = match_ids[0]
                                count_nb_stric_reac_mapped += 1
                                if verbose:
                                    print("reaction: original id:%s; decoded id:%s; match_with:%s from MNX mapping" %(reaction_id, reaction_id_decoded, match_ids[0]))
            if verbose and not match_ids:
                print("reaction: original id:%s; all decoded id:%s; not mapped" %(reaction_id, all_reaction_id_decoded))

    if to_map in ["all", "species"]:
        #For species: k = MNXid, v = {k=db_id,v=[list of ids]} 
        mnx_chem_dict = mnx_reader(mnx_chem_file, db_out)
        species_with_more_one_mapping = 0
        count_nb_stric_species_mapped = 0
  
        for species in listOfSpecies:
            species_id = species.id
            match_ids = None
            all_species_id_decoded = get_all_decoded_version(species_id, "species")
    
            for species_id_decoded in all_species_id_decoded:
                #first check in intern dict mapping
                match_ids = intern_mapping(species_id_decoded, db_out, "species")
        
                if match_ids:
                    dict_sbml_id_mapped_id[species_id] = match_ids
                    count_nb_stric_species_mapped += 1
                    if verbose:
                        print("species: original id:%s; decoded id:%s; match_with:%s from intern mapping" %(species_id, species_id_decoded, match_ids))
                    break
                #check if in mnx_chem_dict
                else:
                    match_ids = get_from_mnx(mnx_chem_dict, species_id_decoded, db_out)
                    if match_ids:
                        if len(match_ids) > 1:
                            species_with_more_one_mapping += 1
                            if verbose:
                                print("species: original id:%s; decoded id:%s; More than one mapping:%s" %(species_id, species_id_decoded, match_ids))
                        else: 
                            dict_sbml_id_mapped_id[species_id] = match_ids[0]
                            count_nb_stric_species_mapped += 1
                            if verbose:
                                print("species: original id:%s; decoded id:%s; match_with:%s from MNX mapping" %(species_id, species_id_decoded, match_ids))
            if verbose and not match_ids:
                print("species: original id:%s; all decoded id:%s; not mapped" %(species_id, all_species_id_decoded))

    if map_from_cpd:
        reaction_mapped_with_cpds = []

        #for all non mapped rxn, check if able to map all speices
        for sbml_rxn in [i for i in listOfReactions if i.id not in list(dict_sbml_id_mapped_id.keys())]:
            all_cpds = set([r.getSpecies() for r in sbml_rxn.getListOfReactants()] + [r.getSpecies() for r in sbml_rxn.getListOfProducts()])
            match_cpd_in_rxn = set([cpd_id for cpd_id in all_cpds if cpd_id in list(dict_sbml_id_mapped_id.keys())])
    
            if len(match_cpd_in_rxn) == len(all_cpds):
                reaction_mapped_with_cpds.append(sbml_rxn.id)

    if verbose:
        print("#######")
        if to_map in ["all", "reaction"]:
            print("Mapped reactions: %s/%s" %(count_nb_stric_reac_mapped,len(listOfReactions)))
            print("Reactions with more than one mapping: %s" %reaction_with_more_one_mapping)
        if to_map in ["all", "species"]:
            print("Mapped species: %s/%s" %(count_nb_stric_species_mapped,len(listOfSpecies)))
            print("Species with more than one mapping: %s" %species_with_more_one_mapping)

        if map_from_cpd:
            print("Mapped reactions from species: %s" %(len(reaction_mapped_with_cpds)))
            for i in reaction_mapped_with_cpds:
                print("\t%s" %i)
            print("Total reactions mapped:%s/%s" %(count_nb_stric_reac_mapped+len(reaction_mapped_with_cpds),len(listOfReactions)))
        print("#######")

    with open(output_dict, 'w') as f:
        for k,v in list(dict_sbml_id_mapped_id.items()):
            f.write(k+"\t"+v+"\n")

def get_from_mnx(mnx_dict, element_id, db_out):
    """
    """
    match_ids = None
    for map_dict in list(mnx_dict.values()):
        #print(map_dict)
        all_element_id = []
        [all_element_id.extend(i) for i in list(map_dict.values())]
        if element_id in all_element_id:
            match_ids = map_dict[db_out]
            break
    return match_ids
    

def mnx_reader(input_file, db_out):
    """
    """
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
    for k,v in list(mnx_dict.items()):
        if db_out not in list(v.keys()) or len(list(v.keys())) == 1:
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
    
    intern_reac_dict = {
    "UNIQ_ID_1":{"METACYC":["RXN-6382"],"KEGG":["R00904"],"UNKNOWN":["APor"]},
    "UNIQ_ID_2":{"METACYC":["CARNOSINE-SYNTHASE-RXN"],"UNKNOWN":["HAL"]},
    "UNIQ_ID_3":{"METACYC":["ALANINE-AMINOTRANSFERASE-RXN"],"KEGG":["R00258"]},
    "UNIQ_ID_4":{"METACYC":["ASPAMINOTRANS-RXN"],"KEGG":["R00355"],"BIGG":["ASPATh"]},
    "UNIQ_ID_5":{"METACYC":["PHEAMINOTRANS-RXN"],"KEGG":["R00694"],"BIGG":["POAT","POATh","POATm"]},
    "UNIQ_ID_6":{"METACYC":["ORNCARBAMTRANSFER-RXN"],"KEGG":["R01398"],"BIGG":["OCT","OCTh","OCTm"]},
    "UNIQ_ID_7":{"METACYC":["3PGAREARR-RXN"],"KEGG":["R01518"],"BIGG":["PGM","PGMf","PGMm"]},
    "UNIQ_ID_8":{"METACYC":["6PGLUCONDEHYDROG-RXN"],"KEGG":["R01528"],"BIGG":["PGDHh"]},
    "UNIQ_ID_9":{"METACYC":["4-HYDROXY-2-KETOPIMELATE-LYSIS-RXN"],"KEGG":["R01645"]},
    "UNIQ_ID_10":{"METACYC":["AMINEPHEN-RXN"],"KEGG":["R02613"]},
    "UNIQ_ID_11":{"METACYC":["SHIKIMATE-5-DEHYDROGENASE-RXN"],"KEGG":["R02413"]},
    "UNIQ_ID_12":{"METACYC":["CARBODEHYDRAT-RXN"],"BIGG":["HCO3E","HCO3Em","HCO3Ehi"]}
    }
    
    intern_chem_dict = {
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
    "UNIQ_ID_21":{"METACYC":["S-3-HYDROXYBUTANOYL-COA"],"BIGG":["3hbcoa"]},
    "UNIQ_ID_22":{"METACYC":["Thiopurines"],"BIGG":["6mpur"]},
    "UNIQ_ID_23":{"METACYC":["N-ACETYL-D-GLUCOSAMINE-6-P"],"BIGG":["acgam6p"]},
    "UNIQ_ID_24":{"METACYC":["CHLOROPHYLL-B"],"BIGG":["chlb"]},
    "UNIQ_ID_25":{"METACYC":["FORMALDEHYDE"],"BIGG":["fald"]},
    "UNIQ_ID_26":{"METACYC":["GLC-1-P"],"BIGG":["g1p"]},
    "UNIQ_ID_27":{"METACYC":["NA+"],"BIGG":["na1"]},

    "UNIQ_ID_28":{"METACYC":["ADP"],"KEGG":["C00008"]},
    "UNIQ_ID_29":{"METACYC":["L-DELTA1-PYRROLINE_5-CARBOXYLATE"],"KEGG":["C03912"]},
    "UNIQ_ID_30":{"METACYC":["AMMONIA"],"KEGG":["C00014"]},
    "UNIQ_ID_31":{"METACYC":["PPI"],"KEGG":["C00013"]},
    "UNIQ_ID_32":{"METACYC":["NADP"],"KEGG":["C00006"]},
    "UNIQ_ID_33":{"METACYC":["NADPH"],"KEGG":["C00005"]},

    "UNIQ_ID_34":{"METACYC":["GLY"],"UNKNOWN":["Glycine"]},
    "UNIQ_ID_35":{"METACYC":["METOH"],"UNKNOWN":["Methanol"]},
    "UNIQ_ID_36":{"METACYC":["PPI"],"UNKNOWN":["Pyrophosphate"]},
    }


    if _type == "reaction":
        for mapp_dict in list(intern_reac_dict.values()):
            all_reaction_id = []
            [all_reaction_id.extend(i) for i in list(mapp_dict.values())]
            if id_to_map in all_reaction_id:
                return mapp_dict.get(db_out, [None])[0]

    elif _type == "species":
        for mapp_dict in list(intern_chem_dict.values()):
            all_species_id = []
            [all_species_id.extend(i) for i in list(mapp_dict.values())]
            if id_to_map in all_species_id:
                return mapp_dict.get(db_out, [None])[0]

    return None


if __name__ == "__main__":
    main()

