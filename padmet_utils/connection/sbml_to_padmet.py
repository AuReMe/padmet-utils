# -*- coding: utf-8 -*-
"""
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

::
    
    usage:
        sbml_to_padmet.py --sbml=DIR/FILE --padmetRef=FILE [--output=FILE] [--db=STR] [--version=STR] [-v]
        sbml_to_padmet.py --sbml=DIR/FILE --padmetSpec=FILE [--output=FILE] [--source_tool=STR] [--source_category=STR] [--source_id=STR] [-v]
        sbml_to_padmet.py --sbml=DIR/FILE --padmetSpec=FILE  --padmetRef=FILE  [--mapping=DIR/FILE] [--mapping_tag=STR] [--output=FILE] [--source_tool=STR] [--source_category=STR] [--source_id=STR] [-v] [-f]
    
    options:
        -h --help     Show help.
        --padmetSpec=FILE    path to the padmet file to update with the sbml. If there's no padmetSpec, just specify the output
        --padmetRef=FILE    path to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
        --sbml=FILE    1 sbml file to convert into padmetSpec (ex: my_network.xml/sbml) OR a directory with n SBML
        --output=FILE   pathanme to the new padmet file
        --mapping=FILE    dictionnary of association id_origin id_ref
        --mapping_tag=STR    if sbml is a folder, use a tag to define mapping files ex: org1.sbml and org1_dict.csv, '_dict.csv' will be the mapping tag. [default: _dict.csv]
        --db=STR    database name
        --version=STR    database version
        -v   print info
"""
from padmet.utils.connection import sbml_to_padmet
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmetRef_file = args["--padmetRef"]
    sbml = args["--sbml"]
    output = args["--output"]
    verbose = args["-v"]
    db = args["--db"]
    if not db: db = "NA"
    version = args["--version"]
    if not version: version = "NA"
    padmetSpec_file = args["--padmetSpec"]
    source_tool = args["--source_tool"]
    source_category = args["--source_category"]
    mapping = args["--mapping"]
    mapping_tag = args["--mapping_tag"]


    if padmetSpec_file:
        sbml_to_padmet.sbml_to_padmetSpec(sbml, padmetSpec_file, padmetRef_file=padmetRef_file, output=output, mapping=mapping, mapping_tag=mapping_tag, source_tool=source_tool, source_category=source_category, db=db, version=version, verbose=verbose)
    else:
        sbml_to_padmet.sbml_to_padmetRef(sbml, padmetRef_file, output, db, version, verbose)

if __name__ == "__main__":
    main()
