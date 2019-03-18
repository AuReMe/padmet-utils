# -*- coding: utf-8 -*-
"""
Description:

    classes.dat:
    For each class:
    create new node / class = class
    UNIQUE-ID (1) => node._id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
    SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name synonyms)
     
    compounds.dat:
    for each compound:
    create new node / class = compound
    UNIQUE-ID (1) => node._id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    INCHI-KEY (0-1) {InChIKey=XXX} => node.misc['INCHI_KEY': XXX]
    MOLECULAR-WEIGHT (0-1) => node.misc()['MOLECULAR_WEIGHT'] = MOLECULAR-WEIGHT
    SMILES (0-1) => node.misc()['SMILES'] = SMILES
    TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
    SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name name)
    DBLINKS (0-n) {(db "id" ...)} => for each, create new node xref, create rlt (node has_xref xref)
    
    proteins.dat:
    for each protein:
    create new node / class = protein
    UNIQUE-ID (1) => node._id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    INCHI-KEY (0-1) {InChIKey=XXX} => node.misc['INCHI_KEY': XXX]
    MOLECULAR-WEIGHT (0-1) => node.misc()['MOLECULAR_WEIGHT'] = MOLECULAR-WEIGHT
    SMILES (0-1) => node.misc()['SMILES'] = SMILES
    TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
    SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name name)
    DBLINKS (0-n) {(db "id" ...)} => for each, create new node xref, create rlt (node has_xref xref)
    SPECIES (0-1) => for each, check or create new node class, create rlt (node is_in_species class)
    
    reactions.dat:
    for each reaction:
    create new node / class = reaction + node.misc()["DIRECTION"] = "UNKNOWN" by default
    UNIQUE-ID (1) => node._id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    EC-NUMBER (0-n) => node.Misc['EC-NUMBER'] = EC-NUMBER
    REACTION-DIRECTION (0-1) => node.Misc['DIRECTION'] = reaction-direction, if REVERSIBLE, else: LEFT-TO-RIGHT
    RXN-LOCATIONS (0,n) => node.misc['COMPARTMENT'] = rxn-location
    TYPES (0-n) => check or create new node class, create rlt (node.id is_a_class types's_node.id)
    DBLINKS (0-n) {(db "id" ...)} => create new node xref, create rlt (node has_xref xref's_node.id)
    SYNONYMS (0-n) => create new node name, create rlt (node has_name name's_node.id)
    --
    for LEFT and RIGHT, also check 2 next lines if info about 'coefficient' or 'compartment'
    defaut value: coefficient/stoichiometry = 1, compartment = unknown
    also check if the direction is 'RIGHT-TO-LEFT', if yes, inverse consumes and produces relations
    then change direction to 'LEFT-TO-RIGHT'
    LEFT (1-n) => create rlt (node.id consumes left's_node.id)
    RIGHT (1-n) => create rlt (node.id produces right's_node.id)
    
    enzrxns.dat:
    for each association enzyme/reaction:
    create new rlt / type = catalyses
    ENZYME (1) => stock enzyme as 'enzyme catalyses'
    REACTION (1-n) => for each reaction after, create relation 'enzyme catalyses reaction'
    
    pathways.dat:
    for each pathway:
    create new node / class = pathway
    UNIQUE-ID (1) => node._id = UNIQUE-ID
    TYPES (0-n) => check or create new node class, create rlt (node is_a_class types)
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    DBLINKS (0-n) {(db "id" ...)} => create new node xref, create rlt (node has_xref xref)
    SYNONYMS (0-n) => create new node name, create rlt (node has_name name)
    IN-PATHWAY (0-n) => check or create new node pathway, create rlt (node is_in_pathway name)
    REACTION-LIST (0-n) => check or create new node pathway, create rlt (node is_in_pathway name)
    
    metabolic-reaction.xml: optional
    for each reaction:

::

    usage:
        pgdb_to_padmet.py --output=FILE --directory=DIR --version=V --db=ID [--padmetRef=FILE] [--source=STR] [-v] [-g] [-m]
        pgdb_to_padmet.py --version=V --db=ID --output=FILE --classes_file=FILE --compounds_file=FILE --proteins_file=FILE --reactions_file=FILE --enzrxns_file=FILE --pathways_file=FILE [--genes_file=FILE] [--metabolic_reactions=FILE] [--source=STR] [-v]
    
    options:
        -h --help     Show help.
        --version=V    Xcyc version
        --db=ID    Biocyc database corresponding to the pgdb (metacyc, ecocyc, ...)
        --output=FILE    padmet file corresponding to the DB
        --directory=DIR    directory containg all the .dat files of metacyc (data)
        --padmetRef=FILE    padmet of reference
        --classes_file=FILE   
        --compounds_file=FILE   
        --proteins_file=FILE 
        --reactions_file=FILE 
        --enzrxns_file=FILE 
        --pathways_file=FILE
        --genes_file=FILE
        --metabolic_reactions=FILE
        --source=STR    Tag associated to the source of the reactions, used to ensure traceability
        -m    use the metabolic-reactions.xml file to enhance the database
        -g    use the genes_file (use if its a specie's pgdb, if metacyc, do not use)
        -v   print info

"""
import re
from datetime import datetime
from time import time
from padmet.classes import PadmetRef, PadmetSpec, Node, Relation
import padmet.utils.sbmlPlugin as sbmlPlugin
import libsbml
import docopt


#1 Parsing .dat
#2 creating nodes and relations
def main():
    chronoDepart = time()

    global regex_purge, regex_xref, list_of_relation, verbose, def_compart_in, def_compart_out, with_genes, source
    regex_purge = re.compile("<.*?>|\|")
    regex_xref = re.compile('^\((?P<DB>\S*)\s*"(?P<ID>\S*)"')
    list_of_relation = []
    def_compart_in = "c"
    def_compart_out = "e"
    #parsing args
    args = docopt.docopt(__doc__)
    version = args["--version"]
    db = args["--db"]
    output = args["--output"]
    path = args["--directory"]
    enhanced_db = args["-m"]
    with_genes = args["-g"]
    source = args["--source"]
    if not source: source = 'genome'
    padmetRef_file = args["--padmetRef"]
    if source : source = source.upper()
    
    if path is not None:
        if not path.endswith("/"): path += "/"

        classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file = \
        [(path + _file) for _file in ["classes.dat", "compounds.dat", "proteins.dat", "reactions.dat", "enzrxns.dat", "pathways.dat"]]
        if enhanced_db:
            metabolic_reactions = path + "metabolic-reactions.xml"
        else:
            metabolic_reactions = None
        if with_genes:
            genes_file = path + "genes.dat"
        else:
            genes_file = None

    else:
        classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file, genes_file, metabolic_reactions = \
        args["--classes_file"], args["--compounds_file"], args["--proteins_file"], args["--reactions_file"], args["--enzrxns_file"], args["--pathways_file"], args["--genes_file"], args["--metabolic_reactions"]
    verbose = args["-v"]
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")
    if padmetRef_file:
        padmet = PadmetSpec()
        padmetRef = PadmetRef(padmetRef_file)
        version = padmetRef.info["DB_info"]["version"]
        db = padmetRef.info["DB_info"]["DB"]
        dbNotes = {"PADMET":{"creation":today_date,"version":"2.6"},"DB_info":{"DB":db,"version":version}}
        padmet.setInfo(dbNotes)
        padmet.setPolicy(padmetRef)
        with open(reactions_file, 'rU') as f:
            rxns_id = [line.split(" - ")[1] for line in f.read().splitlines() if line.startswith("UNIQUE-ID")]
        count = 0
        for rxn_id in rxns_id:
            count += 1
            if verbose: print("%s/%s Copy %s" %(count, len(rxns_id), rxn_id))
            try:
                padmet.copyNode(padmetRef, rxn_id)
                reconstructionData_id = rxn_id+"_reconstructionData_"+source
                if reconstructionData_id in list(padmet.dicOfNode.keys()) and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source %s" %(rxn_id, source))
                reconstructionData = {"SOURCE":[source],"TOOL":["PATHWAYTOOLS"],"CATEGORY":["ANNOTATION"]}
                reconstructionData_rlt = Relation(rxn_id,"has_reconstructionData",reconstructionData_id)
                padmet.dicOfNode[reconstructionData_id] = Node("reconstructionData", reconstructionData_id, reconstructionData)
                padmet._addRelation(reconstructionData_rlt)

            except TypeError:
                if verbose: print("%s not in padmetRef" %(rxn_id))


        if verbose: print("parsing genes")
        dict_protein_gene_id = genes_parser(genes_file, padmet)
        if verbose: print("parsing association enzrxns")
        enzrxns_parser(enzrxns_file, padmet, dict_protein_gene_id)

    else:
        #print(verbose,today_date,version, output, classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file)
        #return    
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
        padmet = PadmetRef()
        if verbose: print("setting policy")
        padmet.setPolicy(POLICY_IN_ARRAY)
        if verbose: print("setting dbInfo")
        padmet.setInfo(dbNotes)
    
    
        if verbose: print("parsing classes")
        classes_parser(classes_file, padmet)
    
        if verbose: print("parsing compounds")
        compounds_parser(compounds_file, padmet)
    
        if verbose: print("parsing reactions")
        reactions_parser(reactions_file, padmet)
    
        if verbose: print("parsing pathways")
        pathways_parser(pathways_file, padmet)
    
        if with_genes:
            if verbose: print("parsing genes")
            dict_protein_gene_id = genes_parser(genes_file, padmet)
            if verbose: print("parsing association enzrxns")
            enzrxns_parser(enzrxns_file, padmet, dict_protein_gene_id)
    
        if metabolic_reactions is not None:
            if verbose: print("enhancing db from metabolic-reactions.xml")
            enhance_db(metabolic_reactions, padmet, with_genes)
    
    for rlt in list_of_relation:
        try:
            padmet.dicOfRelationIn[rlt.id_in].append(rlt)
        except KeyError:
            padmet.dicOfRelationIn[rlt.id_in] = [rlt]
        try:
            padmet.dicOfRelationOut[rlt.id_out].append(rlt)
        except KeyError:
            padmet.dicOfRelationOut[rlt.id_out] = [rlt]

    if with_genes:
        #temp, removing or not reactions without gene assoc:
        all_reactions = [node for node in list(padmet.dicOfNode.values()) if node.type == "reaction"]
        rxn_to_del = [r for r in all_reactions if not any([rlt for rlt in padmet.dicOfRelationIn[r.id] if rlt.type == "is_linked_to"])]
        #[padmet.delNode(node.id) for node in rxn_to_del]
        for rxn in rxn_to_del:
            padmet.delNode(rxn.id)
        if verbose: print("%s/%s reactions without gene association deleted" %(len(rxn_to_del), len(all_reactions)))
        all_genes_linked = set([rlt.id_out for rlt in padmet.getAllRelation() if rlt.type == "is_linked_to"])
        all_genes = set([node.id for node in list(padmet.dicOfNode.values()) if node.type == "gene"])
        count = 0
        for gene_id in [g for g in all_genes if g not in all_genes_linked]:
            count += 1
            #if verbose: print("Removing gene without gene assoc %s" %gene_id)
            padmet.dicOfNode.pop(gene_id)
        if verbose: print("%s/%s genes not linked to any reactions deleted" %(count, len(all_genes)))
    else:
        rxns = [node.id for node in list(padmet.dicOfNode.values()) if node.type == "reaction"]
        for rxn_id in rxns:
            cp_rlts = set([rlt.type for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]])
            if len(cp_rlts) == 1:
                padmet.delNode(rxn_id)

    padmet.generateFile(output)

    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print("done in: ", chrono, "s !")

def from_pgdb_to_padmet(pgdb_folder, db, version, arg_verbose, arg_with_genes, arg_source, enhanced_db, padmetRef_file):
    chronoDepart = time()

    global regex_purge, regex_xref, list_of_relation, verbose, def_compart_in, def_compart_out, with_genes, source
    regex_purge = re.compile("<.*?>|\|")
    regex_xref = re.compile('^\((?P<DB>\S*)\s*"(?P<ID>\S*)"')
    list_of_relation = []
    def_compart_in = "c"
    def_compart_out = "e"
    #parsing args
    verbose = arg_verbose
    with_genes = arg_with_genes
    source = arg_source
    if not source: source = 'genome'

    if source : source = source.upper()

    if pgdb_folder is not None:
        if not pgdb_folder.endswith("/"): pgdb_folder += "/"

        classes_file, compounds_file, proteins_file, reactions_file, enzrxns_file, pathways_file = \
        [(pgdb_folder + _file) for _file in ["classes.dat", "compounds.dat", "proteins.dat", "reactions.dat", "enzrxns.dat", "pathways.dat"]]
        if enhanced_db:
            metabolic_reactions = pgdb_folder + "metabolic-reactions.xml"
        else:
            metabolic_reactions = None
        if with_genes:
            genes_file = pgdb_folder + "genes.dat"
        else:
            genes_file = None

    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")
    if padmetRef_file:
        padmet = PadmetSpec()
        padmetRef = PadmetRef(padmetRef_file)
        version = padmetRef.info["DB_info"]["version"]
        db = padmetRef.info["DB_info"]["DB"]
        dbNotes = {"PADMET":{"creation":today_date,"version":"2.6"},"DB_info":{"DB":db,"version":version}}
        padmet.setInfo(dbNotes)
        padmet.setPolicy(padmetRef)
        with open(reactions_file, 'rU') as f:
            rxns_id = [line.split(" - ")[1] for line in f.read().splitlines() if line.startswith("UNIQUE-ID")]
        count = 0
        for rxn_id in rxns_id:
            count += 1
            if verbose: print("%s/%s Copy %s" %(count, len(rxns_id), rxn_id))
            try:
                padmet.copyNode(padmetRef, rxn_id)
                reconstructionData_id = rxn_id+"_reconstructionData_"+source
                if reconstructionData_id in list(padmet.dicOfNode.keys()) and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source %s" %(rxn_id, source))
                reconstructionData = {"SOURCE":[source],"TOOL":["PATHWAYTOOLS"],"CATEGORY":["ANNOTATION"]}
                reconstructionData_rlt = Relation(rxn_id,"has_reconstructionData",reconstructionData_id)
                padmet.dicOfNode[reconstructionData_id] = Node("reconstructionData", reconstructionData_id, reconstructionData)
                padmet._addRelation(reconstructionData_rlt)

            except TypeError:
                if verbose: print("%s not in padmetRef" %(rxn_id))


        if verbose: print("parsing genes")
        dict_protein_gene_id = genes_parser(genes_file, padmet)
        if verbose: print("parsing association enzrxns")
        enzrxns_parser(enzrxns_file, padmet, dict_protein_gene_id)

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
        padmet = PadmetRef()
        if verbose: print("setting policy")
        padmet.setPolicy(POLICY_IN_ARRAY)
        if verbose: print("setting dbInfo")
        padmet.setInfo(dbNotes)
    
    
        if verbose: print("parsing classes")
        classes_parser(classes_file, padmet)
    
        if verbose: print("parsing compounds")
        compounds_parser(compounds_file, padmet)
    
        if verbose: print("parsing reactions")
        reactions_parser(reactions_file, padmet)
    
        if verbose: print("parsing pathways")
        pathways_parser(pathways_file, padmet)
    
        if with_genes:
            if verbose: print("parsing genes")
            dict_protein_gene_id = genes_parser(genes_file, padmet)
            if verbose: print("parsing association enzrxns")
            enzrxns_parser(enzrxns_file, padmet, dict_protein_gene_id)
    
        if metabolic_reactions is not None:
            if verbose: print("enhancing db from metabolic-reactions.xml")
            enhance_db(metabolic_reactions, padmet, with_genes)
    
    for rlt in list_of_relation:
        try:
            padmet.dicOfRelationIn[rlt.id_in].append(rlt)
        except KeyError:
            padmet.dicOfRelationIn[rlt.id_in] = [rlt]
        try:
            padmet.dicOfRelationOut[rlt.id_out].append(rlt)
        except KeyError:
            padmet.dicOfRelationOut[rlt.id_out] = [rlt]

    if with_genes:
        #temp, removing or not reactions without gene assoc:
        all_reactions = [node for node in list(padmet.dicOfNode.values()) if node.type == "reaction"]
        rxn_to_del = [r for r in all_reactions if not any([rlt for rlt in padmet.dicOfRelationIn[r.id] if rlt.type == "is_linked_to"])]
        #[padmet.delNode(node.id) for node in rxn_to_del]
        for rxn in rxn_to_del:
            padmet.delNode(rxn.id)
        if verbose: print("%s/%s reactions without gene association deleted" %(len(rxn_to_del), len(all_reactions)))
        all_genes_linked = set([rlt.id_out for rlt in padmet.getAllRelation() if rlt.type == "is_linked_to"])
        all_genes = set([node.id for node in list(padmet.dicOfNode.values()) if node.type == "gene"])
        count = 0
        for gene_id in [g for g in all_genes if g not in all_genes_linked]:
            count += 1
            #if verbose: print("Removing gene without gene assoc %s" %gene_id)
            padmet.dicOfNode.pop(gene_id)
        if verbose: print("%s/%s genes not linked to any reactions deleted" %(count, len(all_genes)))
    else:
        rxns = [node.id for node in list(padmet.dicOfNode.values()) if node.type == "reaction"]
        for rxn_id in rxns:
            cp_rlts = set([rlt.type for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type in ["consumes","produces"]])
            if len(cp_rlts) == 1:
                padmet.delNode(rxn_id)

    chrono = (time() - chronoDepart)
    partie_entiere, partie_decimale = str(chrono).split('.')
    chrono = ".".join([partie_entiere, partie_decimale[:3]])
    if verbose: print("done in: ", chrono, "s !")

    return padmet

def classes_parser(filePath, padmet):
    """
    from class.dat: get for each class, the UNIQUE-ID, COMMON-NAME, TYPES, SYNONYMS, DBLINKS
    Create a class node with node.id = UNIQUE-ID,  node.misc = {COMMON-NAME:[COMMON-NAMES]}
    - For each types:
    A type is in fact a class. this information is stocked in padmet as: is_a_class relation btw a node and a class_node
    check if the type is already in the padmet
    if not create a new class_node (var: subClass) with subClass_node.id = type
    Create a relation current node is_a_class type
    - For each Synonyms:
    this information is stocked in padmet as: has_name relation btw a node and a name_node
    create a new name_node with name_node.id = class_id+"_names" and name_node.misc = {LABEL:[synonyms]}
    Create a relation current node has_name name_node.id
    - For each DBLINKS:
    DBLINKS is parsed with regex_xref to get the db and the id
    this information is stocked in padmet as: has_xref relation btw a node and a xref_node
    create a new xref_node with xref_node.id = class_id+"_xrefs" and xref_node.misc = {db:[id]}
    Create a relation current node has_xref xref_node.id
    """
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            try:
                #if len of value is 0 then ValueError raised
                attrib, value = line.split(" - ")
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME",\
                "TYPES", "SYNONYMS", "DBLINKS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]
            except ValueError:
                pass
    count = 0
    nb_classes = str(len(list(dict_data.keys())))
    for class_id, dict_values in dict_data.items():
        count += 1
        #if verbose: print("%s/%s\t%s" %(count, nb_classes, class_id))
        class_node = Node("class", class_id)
        padmet.dicOfNode[class_id] = class_node
        try:
            class_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
        except KeyError:
            pass
        try:
            types = dict_values["TYPES"]
            _setType(types, class_id, padmet)
        except KeyError:
            pass
        try:
            syns = dict_values["SYNONYMS"]
            _setSyns(syns, class_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["DBLINKS"]
            _setXrefs(xrefs, class_id, padmet)
        except KeyError:
            pass


def reactions_parser(filePath, padmet):
    """
    from reaction.dat: get for each reaction, the UNIQUE-ID, COMMON-NAME, TYPES, SYNONYMS, DBLINKS
    Create a reaction node with node.id = UNIQUE-ID,  node.misc = {COMMON-NAME:[COMMON-NAMES]}
    - For each types:
    A type is in fact a class. this information is stocked in padmet as: is_a_class relation btw a node and a class_node
    check if the type is already in the padmet
    if not create a new class_node (var: subClass) with subClass_node.id = type
    Create a relation current node is_a_class type
    - For each Synonyms:
    this information is stocked in padmet as: has_name relation btw a node and a name_node
    create a new name_node with name_node.id = reaction_id+"_names" name_node.misc = {LABEL:[synonyms]}
    Create a relation current node has_name name_node.id
    - For each DBLINKS:
    DBLINKS is parsed with regex_xref to get the db and the id
    this information is stocked in padmet as: has_xref relation btw a node and a xref_node
    create a new xref_node with xref_node.id = reaction_id+"_xrefs" and xref_node.misc = {db:[id]}
    Create a relation current node has_xref xref_node.id
    """
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = [line for line in f.read().splitlines() if not line.startswith("#") and not line == "//"]
        index = -1        
        for line in data:
            index += 1
            try:
                #if len of value is 0 then ValueError raised
                attrib, value = line.split(" - ")
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME", "EC-NUMBER", "REACTION-DIRECTION",\
                "TYPES", "SYNONYMS", "DBLINKS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]

                elif attrib in ["LEFT", "RIGHT"]:
                    #set default values
                    compartment = def_compart_in
                    stoichiometry = "1"
                    #check if information about stoechiometry and compartment in line + 1 and line + 2
                    try:
                        first_next_line = data[index+1]
                    #last line of the file, just add the LEFT/RIGHT
                    except IndexError:
                        first_next_line = ""
                    try:
                        second_next_line = data[index+2]
                    except IndexError:
                        second_next_line = ""

                    if first_next_line.startswith("^COEFFICIENT"):
                        stoichiometry = first_next_line.split(" - ")[1]
                        if second_next_line.startswith("^COMPARTMENT"):
                            compartment = second_next_line.split(" - ")[1]
                    if first_next_line.startswith("^COMPARTMENT"):
                        compartment = first_next_line.split(" - ")[1]
                        if second_next_line.startswith("^COEFFICIENT"):
                            stoichiometry = second_next_line.split(" - ")[1]                            
                    
                    #delete all tags
                    compartment = re.sub(regex_purge, "", compartment)
                    stoichiometry = re.sub(regex_purge, "", stoichiometry)
                    try:
                        dict_data[current_id][attrib].append((value, stoichiometry, compartment))
                    except KeyError:
                        dict_data[current_id][attrib] = [(value, stoichiometry, compartment)]

            except ValueError:
                pass

    count = 0
    nb_rxn = str(len(list(dict_data.keys())))
    for rxn_id, dict_values in dict_data.items():
        if "LEFT" in list(dict_values.keys()) or "RIGHT" in list(dict_values.keys()):
            count += 1
            if verbose: print("%s/%s\t%s" %(count, nb_rxn, rxn_id))
            rxn_node = Node("reaction", rxn_id)
            padmet.dicOfNode[rxn_id] = rxn_node
            try:
                rxn_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
            except KeyError:
                pass
            try:
                rxn_node.misc["EC-NUMBER"] = dict_values["EC-NUMBER"]
            except KeyError:
                pass
            try:
                rxn_dir = dict_values["REACTION-DIRECTION"][0]
                if rxn_dir == "REVERSIBLE":
                    rxn_node.misc["DIRECTION"] = ["REVERSIBLE"]
                elif "LEFT-TO-RIGHT" in rxn_dir:
                    #if:LEFT-TO-RIGHT, IRREVERSIBLE-LEFT-TO-RIGHT, PHYSIOL-RIGHT-TO-LEFT
                    rxn_node.misc["DIRECTION"] = ["LEFT-TO-RIGHT"]
                elif "RIGHT-TO-LEFT" in rxn_dir:
                    #Temporarily set direaction as RIGHT-TO-LEFT
                    #then, RIGHT' metabolites will be LEFT and LEFT -> RIGHT
                    #To finish set back DIRECTION to LEFT-TO-RIGHT
                    rxn_node.misc["DIRECTION"] = ["RIGHT-TO-LEFT"] 
            except KeyError:
                rxn_node.misc["DIRECTION"] = ["REVERSIBLE"]
            """
            try:
                rxn_node.misc["COMPARTMENT"] = dict_values["RXN-LOCATIONS"]
            except KeyError:
                pass
            """
            try:
                types = dict_values["TYPES"]
                _setType(types, rxn_id, padmet)
            except KeyError:
                pass
            try:
                syns = dict_values["SYNONYMS"]
                _setSyns(syns, rxn_id, padmet)
            except KeyError:
                pass
            try:
                xrefs = dict_values["DBLINKS"]
                _setXrefs(xrefs, rxn_id, padmet)
            except KeyError:
                pass
            if with_genes:
                reconstructionData_id = rxn_id+"_reconstructionData_"+source
                if reconstructionData_id in list(padmet.dicOfNode.keys()) and verbose:
                    print("Warning: The reaction %s seems to be already added from the same source %s" %(rxn_id, source))
                reconstructionData = {"SOURCE":[source],"TOOL":["PATHWAYTOOLS"],"CATEGORY":["ANNOTATION"]}
                reconstructionData_rlt = Relation(rxn_id,"has_reconstructionData",reconstructionData_id)
                padmet.dicOfNode[reconstructionData_id] = Node("reconstructionData", reconstructionData_id, reconstructionData)
                list_of_relation.append(reconstructionData_rlt)
            try:
                reactants_data = dict_values["LEFT"]
                for reactant_id, stoichiometry, compartment in reactants_data:
                    if compartment == "CCO-OUT":
                        compartment = def_compart_out
                    else:
                        compartment = def_compart_in
                    try:
                        reactant_node = padmet.dicOfNode[reactant_id]
                    except KeyError:
                        reactant_node = Node("compound", reactant_id)
                        padmet.dicOfNode[reactant_id] = reactant_node
    
                    #if the reaction direction was set to RIGHT-TO-LEFT, then this compound is in fact a product
                    if rxn_node.misc["DIRECTION"][0] == "RIGHT-TO-LEFT":
                        produces_rlt = Relation(rxn_id, "produces", reactant_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(produces_rlt)
                    else:
                        consumes_rlt = Relation(rxn_id, "consumes", reactant_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(consumes_rlt)
            except KeyError:
                pass
            try:
                products_data = dict_values["RIGHT"]
                for product_id, stoichiometry, compartment in products_data:
                    if compartment == "CCO-OUT":
                        compartment = def_compart_out
                    else:
                        compartment = def_compart_in
                    try:
                        product_node = padmet.dicOfNode[product_id]
                    except KeyError:
                        product_node = Node("compound", product_id)
                        padmet.dicOfNode[product_id] = product_node

                    #if the reaction direction was set to RIGHT-TO-LEFT, then this compound is in fact a reactant
                    if rxn_node.misc["DIRECTION"][0] == "RIGHT-TO-LEFT":
                        consumes_rlt = Relation(rxn_id, "consumes", product_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(consumes_rlt)
                    else:
                        produces_rlt = Relation(rxn_id, "produces", product_id, {"STOICHIOMETRY": [stoichiometry], "COMPARTMENT": [compartment]})
                        list_of_relation.append(produces_rlt)
            except KeyError:
                pass
            if rxn_node.misc["DIRECTION"][0] == "RIGHT-TO-LEFT":
                rxn_node.misc["DIRECTION"] = ["LEFT-TO-RIGHT"]

def pathways_parser(filePath, padmet):
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            try:
                #if len of value is 0 then ValueError raised
                attrib, value = line.split(" - ")
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME", "TAXONOMIC-RANGE",\
                "TYPES", "SYNONYMS", "DBLINKS", "IN-PATHWAY", "REACTION-LIST"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]
            except ValueError:
                pass
    
    count = 0
    nb_pathways = str(len(list(dict_data.keys())))
    for pathway_id, dict_values in dict_data.items():
        count += 1
        #if verbose: print("%s/%s\t%s" %(count, nb_pathways, pathway_id))
        pathway_node = Node("pathway", pathway_id)
        padmet.dicOfNode[pathway_id] = pathway_node
        try:
            pathway_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
        except KeyError:
            pass
        try:
            pathway_node.misc["TAXONOMIC-RANGE"] = dict_values["TAXONOMIC-RANGE"]
        except KeyError:
            pass
        try:
            types = dict_values["TYPES"]
            _setType(types, pathway_id, padmet)
        except KeyError:
            pass
        try:
            syns = dict_values["SYNONYMS"]
            _setSyns(syns, pathway_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["DBLINKS"]
            _setXrefs(xrefs, pathway_id, padmet)
        except KeyError:
            pass
        try:
            subPathways = dict_values["IN-PATHWAY"]
            for subPathway in subPathways:
                #add the hierachization info, current pathway is_in_pathway subpathway
                is_in_pathway_rlt = Relation(pathway_id, "is_in_pathway", subPathway)
                list_of_relation.append(is_in_pathway_rlt)
        except KeyError:
            pass
        try:
            subNodes = dict_values["REACTION-LIST"]
            for subNode in subNodes:
                #add the hierachization info, Reaction/pathway is_in_pathway current pathway
                is_in_pathway_rlt = Relation(subNode, "is_in_pathway", pathway_id)
                list_of_relation.append(is_in_pathway_rlt)
        except KeyError:
            pass

def compounds_parser(filePath, padmet):
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            try:
                #if len of value is 0 then ValueError raised
                attrib, value = line.split(" - ")
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","INCHI-KEY","MOLECULAR-WEIGHT","SMILES",\
                "TYPES", "SYNONYMS", "DBLINKS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]
            except ValueError:
                pass

    count = 0
    nb_cpds = str(len(list(dict_data.keys())))
    for compound_id, dict_values in dict_data.items():
        count += 1
        #if verbose: print("%s/%s\t%s" %(count, nb_cpds, compound_id))
        compound_node = Node("compound", compound_id)
        padmet.dicOfNode[compound_id] = compound_node
        try:
            compound_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
        except KeyError:
            pass
        try:
            compound_node.misc["INCHI-KEY"] = dict_values["INCHI-KEY"]
        except KeyError:
            pass
        try:
            compound_node.misc["MOLECULAR-WEIGHT"] = dict_values["MOLECULAR-WEIGHT"]
        except KeyError:
            pass
        try:
            compound_node.misc["SMILES"] = dict_values["SMILES"]
        except KeyError:
            pass
        try:
            types = dict_values["TYPES"]
            _setType(types,compound_id, padmet)
        except KeyError:
            pass
        try:
            syns = dict_values["SYNONYMS"]
            _setSyns(syns, compound_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["DBLINKS"]
            _setXrefs(xrefs, compound_id, padmet)
        except KeyError:
            pass

def genes_parser(filePath, padmet):
    dict_data = {}
    #k='ACCESSION-1', v ='PRODUCT'
    dict_protein_gene_id = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            try:
                #if len of value is 0 then ValueError raised
                attrib, value = line.split(" - ")
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","ACCESSION-1","CENTISOME-POSITION","LEFT-END-POSITION","RIGHT-END-POSITION","SYNONYMS","TRANSCRIPTION-DIRECTION","PRODUCT"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]
            except ValueError:
                pass

    count = 0
    nb_genes = str(len(list(dict_data.keys())))
    for current_id, dict_values in dict_data.items():
        count += 1
        try:
            gene_id = dict_values.get("ACCESSION-1",[current_id])[0]
            enzyme_id = dict_values["PRODUCT"][0]
            dict_protein_gene_id[enzyme_id] = gene_id
            gene_node = Node("gene", gene_id)
            padmet.dicOfNode[gene_id] = gene_node
            try:
                if dict_values["COMMON-NAME"][0] != gene_id:
                    gene_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
            except KeyError:
                pass
            try:
                if dict_values["TRANSCRIPTION-DIRECTION"][0] == "-":
                    gene_node.misc["TRANSCRIPTION-DIRECTION"] = ["NEGATIVE"]
                elif dict_values["TRANSCRIPTION-DIRECTION"][0] == "+":
                    gene_node.misc["TRANSCRIPTION-DIRECTION"] = ["POSITIVE"]
            except KeyError:
                pass
            try:
                gene_node.misc["CENTISOME-POSITION"] = dict_values["CENTISOME-POSITION"]
            except KeyError:
                pass
            try:
                gene_node.misc["LEFT-END-POSITION"] = dict_values["LEFT-END-POSITION"]
            except KeyError:
                pass
            try:
                gene_node.misc["RIGHT-END-POSITION"] = dict_values["RIGHT-END-POSITION"]
            except KeyError:
                pass
            try:
                gene_node.misc["RIGHT-END-POSITION"] = dict_values["RIGHT-END-POSITION"]
            except KeyError:
                pass
            try:
                syns = dict_values["SYNONYMS"]
                _setSyns(syns, gene_id, padmet)
            except KeyError:
                pass
        except KeyError:
            pass
        if verbose: print("%s/%s\t%s/%s" %(count, nb_genes, current_id, gene_id))
    return dict_protein_gene_id

def proteins_parser(filePath, padmet, dict_gene_unique_id_real_id = None):
    dict_data = {}
    dict_id_protein_gene_real = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            try:
                #if len of value is 0 then ValueError raised
                attrib, value = line.split(" - ")
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","INCHI-KEY","MOLECULAR-WEIGHT","SMILES",\
                "TYPES", "SPECIES", "SYNONYMS", "DBLINKS", "GENE", "GO-TERMS"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]
            except ValueError:
                pass
            

    count = 0
    nb_proteins = str(len(list(dict_data.keys())))
    for protein_id, dict_values in dict_data.items():
        count += 1
        if verbose: print("%s/%s\t%s" %(count, nb_proteins, protein_id))
        protein_node = Node("protein", protein_id)
        padmet.dicOfNode[protein_id] = protein_node
        try:
            protein_node.misc["COMMON-NAME"] = dict_values["COMMON-NAME"]
        except KeyError:
            pass
        try:
            protein_node.misc["INCHI-KEY"] = dict_values["INCHI-KEY"]
        except KeyError:
            pass
        try:
            protein_node.misc["MOLECULAR-WEIGHT"] = dict_values["MOLECULAR-WEIGHT"]
        except KeyError:
            pass
        try:
            protein_node.misc["SMILES"] = dict_values["SMILES"]
        except KeyError:
            pass
        try:
            types = dict_values["TYPES"]
            _setType(types, protein_id, padmet)
        except KeyError:
            pass
        try:
            species = dict_values["SPECIES"]
            _setType(species, protein_id, padmet)
        except KeyError:
            pass
        try:
            syns = dict_values["SYNONYMS"]
            _setSyns(syns, protein_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["DBLINKS"]
            _setXrefs(xrefs, protein_id, padmet)
        except KeyError:
            pass
        try:
            xrefs = dict_values["GO-TERMS"]
            _setXrefs(xrefs, protein_id, padmet)
        except KeyError:
            pass
        if dict_gene_unique_id_real_id is not None:
            try:
                genes = dict_values["GENE"]
                for gene in genes:
                    try:
                        gene_id = dict_gene_unique_id_real_id[gene]
                    except KeyError:
                        gene_id = gene
                        gene_node = Node("gene", gene_id)                        
                        padmet.dicOfNode[gene_id] = gene_node
                    try:
                        dict_id_protein_gene_real[protein_id].append(gene_id)
                    except KeyError:
                        dict_id_protein_gene_real[protein_id] = [gene_id]
    
                    codes_for_rlt = Relation(gene_id, "codes_for", protein_id)
                    list_of_relation.append(codes_for_rlt)
            except KeyError:
                pass
    return dict_id_protein_gene_real

def enzrxns_parser(filePath, padmet, dict_protein_gene_id = None):
    dict_data = {}
    with open(filePath, 'r', encoding='windows-1252') as f:
        data = (line for line in f.read().splitlines() if not line.startswith("#") and not line == "//")
        for line in data:
            try:
                #if len of value is 0 then ValueError raised
                attrib, value = line.split(" - ")
                #delete all tags
                value = re.sub(regex_purge,"",value)
                if attrib == "UNIQUE-ID":
                    current_id = value
                    dict_data[current_id] = {}
                if attrib in ["COMMON-NAME","ENZYME","REACTION","BASIS-FOR-ASSIGNMENT"]:
                    try:
                        dict_data[current_id][attrib].append(value)
                    except KeyError:
                        dict_data[current_id][attrib] = [value]
            except ValueError:
                pass
    for k,v in list(dict_data.items()):
        if len(v["ENZYME"]) > 1:
            print(k)
    count = 0
    nb_enzrxns = str(len(list(dict_data.keys())))
    for current_id, dict_values in dict_data.items():
        count += 1
        #if verbose: print("%s/%s\t%s" %(count, nb_enzrxns, current_id))
        rxn_id = dict_values["REACTION"][0]
        names = dict_values.get("COMMON-NAME",[])
        for name in names: 
            if name.endswith("_"):
                names[names.index(name)] = name[:-1] 
        protein = dict_values["ENZYME"][0]
        try:
            rxn_node = padmet.dicOfNode[rxn_id]
            try:
                [rxn_node.misc["COMMON-NAME"].append(name) for name in names if name not in rxn_node.misc["COMMON-NAME"]]
            except KeyError:
                rxn_node.misc["COMMON-NAME"] = names
            if dict_protein_gene_id is not None:
                try:
                    gene_id = dict_protein_gene_id[protein]
                    try:
                        assignment = dict_values["BASIS-FOR-ASSIGNMENT"][0]
                        if assignment.startswith(":"): assignment = assignment[1:]
                    except KeyError:
                        assignment = "NA"
                    is_linked_rlt = Relation(rxn_id, "is_linked_to", gene_id, {"SOURCE:ASSIGNMENT":[source+":"+assignment]})
                    list_of_relation.append(is_linked_rlt)
                except KeyError:
                    pass
        except KeyError:
            pass


def _setType(types, current_id, padmet):
    for subClass_id in types:
        #Type allow the hierachization of the current
        #XX is_a_class type
        #if type is not already in padmet, create a new class node (subClass_node)
        #subClass_node id == type
        #create a relation xx is_a_class type
        try:
            subClass_node = padmet.dicOfNode[subClass_id]
        except KeyError:
            subClass_node = Node("class", subClass_id)
            padmet.dicOfNode[subClass_id] = subClass_node
        is_a_class_rlt = Relation(current_id, "is_a_class", subClass_id)
        list_of_relation.append(is_a_class_rlt)
    
def _setSyns(syns, current_id, padmet):
    name_id = current_id+"_names"
    try:
        name_node = padmet.dicOfNode[name_id]
    except KeyError:
        #create node name
        name_node = Node("name", name_id, {"LABEL":[]})
        padmet.dicOfNode[name_id] = name_node
        has_name_rlt = Relation(current_id, "has_name", name_id)
        list_of_relation.append(has_name_rlt)
    [name_node.misc["LABEL"].append(syn) for syn in syns 
    if syn not in name_node.misc["LABEL"]]

def _setXrefs(xrefs, current_id, padmet):
    xref_id = current_id+"_xrefs"
    try:
        xref_node = padmet.dicOfNode[xref_id]
    except KeyError:
        #create node xref
        xref_node = Node("xref", xref_id)
        padmet.dicOfNode[xref_id] = xref_node
        has_xref_rlt = Relation(current_id, "has_xref", xref_id)
        list_of_relation.append(has_xref_rlt)

    for xref in xrefs:
        #an xref is like: (REFSEQ "NP_417401" NIL NIL NIL NIL NIL)
        #in this example DB = REFSEQ and ID = NP_417401
        #update node xref, with in misc k = DB and v = [ID]
        #node id is created by incrementing meta_max_id
        xref_search = regex_xref.search(xref)
        if xref_search is not None:
            xref_dict = xref_search.groupdict()
            db = xref_dict["DB"]
            _id = xref_dict["ID"]
        else:
            db = "GO-TERMS"
            _id = xref
        if db in list(xref_node.misc.keys()) and _id not in xref_node.misc[db]:
            xref_node.misc[db].append(_id)
        else:
            xref_node.misc[db] = [_id]

def enhance_db(metabolic_reactions, padmet, with_genes=False):
    print("loading sbml file: %s" %metabolic_reactions)
    reader = libsbml.SBMLReader()
    document = reader.readSBML(metabolic_reactions)
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    #recovere the reactions that are not in the basic metacyc but in the sbml file
    #use the reactions_name instead of ids because the ids are encoded, the name is the non-encoded version of the id
    padmet_reactions_id = set([node.id for node in list(padmet.dicOfNode.values()) if node.type == "reaction"])
    reaction_to_add = [reaction for reaction in listOfReactions 
    if reaction.getName() not in padmet_reactions_id]
    count = 0
    if verbose: print(str(len(reaction_to_add))+" reactions to add")
    for reactionSBML in reaction_to_add:
        count += 1
        reaction_id = reactionSBML.getName()
        if verbose: print(str(count)+"/"+str(len(reaction_to_add))+"\t"+reaction_id)
        if reactionSBML.getReversible():
            reaction_dir = "REVERSIBLE"
        else:
            reaction_dir = "LEFT-TO-RIGHT"
        try:
            reaction_node = padmet.dicOfNode[reaction_id]
        except KeyError:
            reaction_node = Node("reaction", reaction_id, {"DIRECTION": [reaction_dir]})
            padmet.dicOfNode[reaction_id] = reaction_node
        reactants = reactionSBML.getListOfReactants()
        for reactant in reactants: #convert ids
            reactant_id, _type, reactant_compart = sbmlPlugin.convert_from_coded_id(reactant.getSpecies())
            if reactant_id not in list(padmet.dicOfNode.keys()):
                reactant_node = Node("compound",reactant_id)
                padmet.dicOfNode[reaction_id] = reactant_node
            reactant_stoich = reactant.getStoichiometry()
            consumes_rlt = Relation(reaction_id,"consumes",reactant_id, {"STOICHIOMETRY":[reactant_stoich], "COMPARTMENT": [reactant_compart]})
            list_of_relation.append(consumes_rlt)

        products = reactionSBML.getListOfProducts()
        for product in products:
            product_id, _type, product_compart = sbmlPlugin.convert_from_coded_id(product.getSpecies())
            if product_id not in list(padmet.dicOfNode.keys()):
                product_node = Node("compound",product_id)
                padmet.dicOfNode[product_id] = product_node
            product_stoich = product.getStoichiometry()
            produces_rlt = Relation(reaction_id,"produces",product_id,{"STOICHIOMETRY": [product_stoich], "COMPARTMENT": [product_compart]})
            list_of_relation.append(produces_rlt)
        
        if with_genes:
            notes = sbmlPlugin.parseNotes(reactionSBML)
            if "GENE_ASSOCIATION" in list(notes.keys()):
                #Using sbmlPlugin to recover all genes associated to the reaction
                listOfGenes = sbmlPlugin.parseGeneAssoc(notes["GENE_ASSOCIATION"][0])
                if len(listOfGenes) != 0:
                    for gene in listOfGenes:
                        try:
                            #check if gene already in the padmet
                            padmet.dicOfNode[gene]
                        except TypeError:
                            gene_node = Node("gene",gene)
                            padmet.dicOfNode[gene] = gene_node
                        is_linked_rlt = Relation(reaction_id, "is_linked_to", gene)
                        list_of_relation.append(is_linked_rlt)



if __name__ == "__main__":
    main()      
