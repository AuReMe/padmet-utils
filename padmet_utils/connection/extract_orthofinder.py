# -*- coding: utf-8 -*-
"""
Description:
    After running orthofinder on n fasta file, read the output file 'Orthogroups.csv'
    
    Require a folder 'orthology_based_folder' with this archi:
    
        |-- model_a
            -- model_a.sbml
        |-- model_b
            --model_b.sbml

    And the name of the studied organism 'study_id'

    1. Read the orthogroups file, extract orthogroups in dict 'all_orthogroups', and all org names

    2. In orthology folder search for sbml files 'extension = .sbml'

    3. For each models regroup all information in a dict dict_data:
        
        {'study_id': study_id,
        'model_id' : model_id,
        'sbml_template': path to sbml of model',
        'output': path to the output sbml,
        'verbose': bool, if true print information
        }

        The output is by default:
            \output_orthofinder_from_'model_id'.sbml

    4. Store all previous dict_data in a list all_dict_data

    5. iter on dict from all_dict_data and use function dict_data_to_sbml

    Use a dict of data dict_data and dict of orthogroups dict_orthogroup to create sbml files.
    
    dict_data and dict_orthogroup are obtained with fun orthofinder_to_sbml
    
    6./ Read dict_orthogroups and check if model associated to dict_data and study org share orthologue
    
    7./ Read sbml of model, parse all reactions and get genes associated to reaction.
    
    8./ For each reactions:
        
        Parse genes associated to sub part (ex: (gene-a and gene-b) or gene-c) = [(gene-a,gene-b), gene-c]
        
        Check if study org have orthologue with at least one sub part (gene-a, gene-b) or gene-c
        
        if yes: add the reaction to the new sbml and change genes ids by study org genes ids
        
        Create the new sbml file.
    
::
   
    usage:
        extract_orthofinder --sbml=FILE/DIR --orthologues=DIR --study_id=STR --output=DIR [--workflow=STR] [-v]
        extract_orthofinder --sbml=DIR --orthogroups=FILE --study_id=STR --output=DIR [--workflow=STR] [-v]

    option:
        -h --help    Show help.
        --sbml=DIR   Folder with sub folder named as models name within sbml file name as model_name.sbml
        --orthogroups=FILE   Output file of Orthofinder run Orthogroups.tsv
        --orthologues=DIR   Output directory of Orthofinder run Orthologues
        --study_id=ID   name of the studied organism
        --workflow=ID   worklow id in ['aureme','aucome']. specific run architecture where to search sbml files
       --output=DIR   folder where to create all sbml output files
        -v   print info

"""
import docopt
from padmet.utils.connection import extract_orthofinder


def main():
    args = docopt.docopt(__doc__)
    verbose = args["-v"]
    sbml = args["--sbml"]
    orthogroups_file = args["--orthogroups"]
    orthologue_folder = args["--orthologues"]
    output_folder = args["--output"]
    study_id = args["--study_id"]
    workflow = args["--workflow"]
    all_model_sbml = extract_orthofinder.get_sbml_files(sbml, workflow, verbose)
    if orthogroups_file:
        extract_orthofinder.orthogroups_to_sbml(orthogroups_file, all_model_sbml, output_folder, study_id, verbose)
    elif orthologue_folder:
        extract_orthofinder.orthologue_to_sbml(orthologue_folder, all_model_sbml, output_folder, study_id, verbose)

if __name__ == "__main__":
    main()

