# -*- coding: utf-8 -*-
"""
Description:
    Update a padmetSpec by filling specific forms.

    1./ Create new reaction(s) to padmet file. \n
    - Get the template form with --template_new_rxn
    - Fill the template
    - set --data as path to the filled template

    2./ Add reaction(s) from padmetRef or remove reactions(s). \n
    - Get the template form with --template_add_delete_rxn
    - Fill the template
    - set --date as path to the filled template

    Update padmetSpec and create a new padmet (new_padmet) or overwrite the input

::

    usage:
        manual_curation.py --padmetSpec=FILE --data=FILE [--padmetRef=FILE] [--output=FILE] [--tool=STR] [--category=STR] [-v]
        manual_curation.py --template_new_rxn=FILE
        manual_curation.py --template_add_delete_rxn=FILE

    option:
        -h --help    Show help.
        --padmetSpec=FILE    path to the padmet to update
        --padmetRef=FILE    path of the padmet representing the reference database
        --data=FILE    path to the form with data for curation
        --output=FILE    path to the output. if None. Overwriting padmetSpec
        --tool=STR    specification of the tool used to allow this curation: ex a tool of gapfilling (meneco)
        --category=STR    specification of the category of curation: ex if a reaction is added based on annotation info, use 'annotation'
        --template_new_rxn=FILE    create a form used to create new reaction, use this form as input for 'data' option
        --template_add_delete_rxn=FILE    create a form used to add or delete reaction, use this form as input for 'data' option
        -v    print info
"""
import os
import docopt
from padmet.classes import PadmetSpec
from padmet.classes import PadmetRef
from padmet.utils.management import manual_curation

def main():
    args = docopt.docopt(__doc__)
    data_file = args["--data"]
    output = args["--output"]
    verbose = args["-v"]
    if data_file:
        filename = os.path.splitext(os.path.basename(data_file))[0]
        source = filename

    category = args["--category"]
    tool = args["--tool"]
    if args["--template_new_rxn"]:
        output = args["--template_new_rxn"]
        manual_curation.template_new_rxn(output)
    elif args["--template_add_delete_rxn"]:
        output = args["--template_add_delete_rxn"]
        manual_curation.template_add_delete(output)
    else:
        padmetSpec = PadmetSpec(args["--padmetSpec"])
        if not output:
            output = args["--padmetSpec"]
        if args["--padmetRef"]:
            padmetRef = PadmetRef(args["--padmetRef"])
        else:
            padmetRef = None
        to_do = manual_curation.sniff_datafile(data_file)

        if to_do == "rxn_creator":
            manual_curation.rxn_creator(data_file, padmetSpec, output, padmetRef, source, tool, category, verbose)
        elif to_do == "add_delete_rxn":
            manual_curation.add_delete_rxn(data_file, padmetSpec, output, padmetRef, source, tool, category, verbose)

if __name__ == "__main__":
    main()
