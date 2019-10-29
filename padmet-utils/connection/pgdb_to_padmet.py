# -*- coding: utf-8 -*-
"""
Description:

    Read a PGDB folder (from BIOCYC/PATHWAYTOOLS) and create a padmet.
    1./ To create a padmet without any genes information extracted use the first usage with:
        pgdb: path to pgdb folder
        output: path to the padmet to create
        version: to specify the version of the pgdb (20.0, 22.0)
        db: to sepcify the name of the database (METACYC, ECOCYC, ...)
        enhance: to also read the file metabolic-reaction.xml and add the to the padmet
    2./ To create a padmet and add only reactions from pgdb if they are in padmetRef specifie.
        Copy information of the reaction not from the pgdb but from the padmetRef.
        This allow to uniform reaction to the same version of metacyc represented in the padmetRef
        For example, in some case 2 pgdb from different version can contain different information for a same reaction,pathway...
        In this case use:
            padmetRef: path to the padmet of reference
    3./ To create a padmet wth genes information extracted use:
        extract-gene
    3.1/ To remove from the final padmet all reactions without genes associated use:
        no-orphan
    4./ To read the metabolic-reaction.xml file, a sbml with some missing reactions in PGDB use:
        enhance

    For more information of the parsing process read information below.

    
    
    classes.dat:
    For each class:
    create new node / class = class
    UNIQUE-ID (1) => node.id = UNIQUE-ID
    COMMON-NAME (0-n) => node.Misc['COMMON-NAME'] = COMMON-NAME
    TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
    SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name synonyms)
     
    compounds.dat:
    for each compound:
    create new node / class = compound
    UNIQUE-ID (1) => node.id = UNIQUE-ID
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
    UNIQUE-ID (1) => node.id = UNIQUE-ID
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
    UNIQUE-ID (1) => node.id = UNIQUE-ID
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
    

::

    usage:
        pgdb_to_padmet.py --pgdb=DIR --output=FILE [--version=V] [--db=ID] [--padmetRef=FILE] [--source=STR] [-v] [--enhance]
        pgdb_to_padmet.py --pgdb=DIR --output=FILE --extract-gene [--no-orphan]  [--version=V] [--db=ID] [--padmetRef=FILE] [--source=STR] [-v] [--enhance]
    
    options:
        -h --help     Show help.
        --version=V    Xcyc version [default: N.A].
        --db=ID    Biocyc database corresponding to the pgdb (metacyc, ecocyc, ...) [default: N.A].
        --output=FILE    padmet file corresponding to the DB.
        --pgdb=DIR    directory containg all the .dat files of metacyc (data).
        --padmetRef=FILE    padmet of reference.
        --source=STR    Tag associated to the source of the reactions, used to ensure traceability [default: GENOME].
        --enhance    use the metabolic-reactions.xml file to enhance the database.
        --extract-gene    use the genes_file (use if its a specie's pgdb, if metacyc, do not use).
        --no-orhpan    use the genes_file (use if its a specie's pgdb, if metacyc, do not use).
        -v   print info.

"""
from padmet.utils.connection import pgdb_to_padmet
import docopt

#1 Parsing .dat
#2 creating nodes and relations
def main():
    #parsing args
    args = docopt.docopt(__doc__)
    version = args["--version"]
    db = args["--db"]
    source = args["--source"]
    output = args["--output"]
    pgdb_folder = args["--pgdb"]
    enhanced_db = args["--enhance"]
    extract_gene = args["--extract-gene"]
    no_orphan = args["--no-orphan"]
    padmetRef_file = args["--padmetRef"]
    verbose = args["-v"]
    padmet = pgdb_to_padmet.from_pgdb_to_padmet(pgdb_folder, db , version, source, extract_gene, no_orphan, enhanced_db, padmetRef_file, verbose)
    padmet.generateFile(output)

if __name__ == "__main__":
    main()      
