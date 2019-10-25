# -*- coding: utf-8 -*-
"""
Description:
    This tool is use the MetaNetX database to check or convert a sbml. Flat files
    from MetaNetx are required to run this tool. They can be found in the aureme workflow
    or from the MetaNetx website.
    To use the tool set:
        mnx_folder= the path to a folder containing MetaNetx flat files.
        the files must be named as 'reac_xref.tsv' and 'chem_xref.tsv'
        or set manually the different path of the flat files with:
            mnx_reac= path to the flat file for reactions
            
            mnx_chem= path to the flat file for chemical compounds (species)

    To check the database used in a sbml:
        to check all element of sbml (reaction and species) set:
            to--map=all
        to check only reaction of sbml set:
            to--map=reaction
        to check only species of sbml set:
            to--map=species

    To map a sbml and obtain a file of mapping ids to a given database set:
        to-map:
            as previously explained
        db_out:
            the name of the database target: ['metacyc', 'bigg', 'kegg'] only
        output:
            the path to the output file
        
        For a given sbml using a specific database.
        
        Return a dictionnary of mapping.
        
        the output is a file with line = reaction_id/or species in sbml, reaction_id/species in db_out database
        
        ex:
            For a sbml based on kegg database, db_out=metacyc: the output file will contains for ex:
        R02283 ACETYLORNTRANSAM-RXN
    
::
    
    usage:
        convert_sbml_db.py --mnx_reac=FILE --mnx_chem=FILE --sbml=FILE --to-map=STR [-v]
        convert_sbml_db.py --mnx_folder=DIR --sbml=FILE --to-map=STR [-v]
        convert_sbml_db.py --mnx_folder=DIR --sbml=FILE --output=FILE --db_out=ID --to-map=STR [-v]
        convert_sbml_db.py --mnx_reac=FILE --mnx_chem=FILE --sbml=FILE --output=FILE --db_out=ID --to-map=STR [-v]
    
    options:
        -h --help     Show help.
        --to-map=STR     select the part of the sbml to check or convert, must be in ['all', 'reaction', 'species']
        --mnx_reac=FILE     path to the MetaNetX file for reactions
        --mnx_chem=FILE     path to the MetaNetX file for compounds
        --sbml=FILE     path to the sbml file to convert
        --output=FILE     path to the file containing the mapping, sep = "\t"
        --db_out=FILE     id of the output database in ["BIGG","METACYC","KEGG"]
        -v     verbose.

"""
from padmet.utils.exploration import convert_sbml_db
import docopt

def main():
    args = docopt.docopt(__doc__)
    verbose = args["-v"]
    to_map = args["--to-map"]
    mnx_folder = args["--mnx_folder"]
    mnx_reac_file = args["--mnx_reac"]
    mnx_chem_file = args["--mnx_chem"]
    sbml_file = args["--sbml"]

    if args["--db_out"]:
        db_out = args["--db_out"].upper()
        output = args["--output"]
        convert_sbml_db.map_sbml(sbml_file, to_map, db_out, output, verbose, mnx_reac_file, mnx_chem_file, mnx_folder)
    else:
        db_select, db_found = convert_sbml_db.check_sbml_db(sbml_file, to_map, verbose, mnx_reac_file, mnx_chem_file, mnx_folder)
        print("Best matching database: %s" %db_select)
        print(db_found)


if __name__ == "__main__":
    main()

