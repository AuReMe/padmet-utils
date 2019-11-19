#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    usage:
        relation_curation.py    --padmet=FILE --id_in=STR [--type=STR] [-v]
        relation_curation.py    --padmet=FILE --id_out=STR [--type=STR] [-v]
        relation_curation.py    --padmet=FILE --id_in=STR [--type=STR] --to-remove==STR --output=FILE [-v]
        relation_curation.py    --padmet=FILE --id_in=STR --id_out=STR [--type=STR] --to-remove==STR --output=FILE [-v]
    
    option:
        -h --help    Show help.
        --model_metabolic=FILE    pathname to the metabolic network of the model (sbml).
        --study_metabolic=FILE    ****.
        --inp=FILE    ****.
        --omcl=FILE    ****.
        --output=FILE    ****.
        -v   print info.

"""
import docopt
from padmet.classes import PadmetSpec
from padmet.utils.management import relation_curation

def main():
    args = docopt.docopt(__doc__)
    padmet_file = args["--padmet"]
    id_in = args["--id_in"]
    id_out = args["--id_out"]
    _type = args["--type"]
    output = args["--output"]
    verbose = args["-v"]
    to_remove = args["--to-remove"]
    padmet = PadmetSpec(padmet_file)
    relation_curation.get_relations(padmet=padmet, id_in=id_in, id_out=id_out, _type=_type, verbose=verbose)

    if to_remove:
        if to_remove != "all":
            to_remove = to_remove.split(";")
        relation_curation.get_relations(padmet=padmet, id_in=id_in, id_out=id_out, _type=_type, to_remove=to_remove, output=output, verbose=False)


if __name__ == "__main__":
    main()
