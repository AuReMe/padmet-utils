# -*- coding: utf-8 -*-
"""
Description:
    The standard output of meneco return ids of reactions corresponding to the solution for gapfilling.

    The ids are those from the sbml and so they are encoded.

    This script extract the solution corresponding to the union of reactions
    "Computing union of reactions from all completion"
    Based on padmetRef return a file with more information for each reaction.

    ex: RXN__45__5

    RXN-5, common_name, ec-number, Formula (with id),Formula (with cname),Action,Comment
    Also, the output can be used as input of the script update_padmetSpec.py
    In the column Action: 'add' => To add the reaction, '' => to do nothing

    Comment: the reason of adding the reaction (ex: added for gap-filling by meneco)

::
    
    usage:
        enhanced_meneco_output.py --meneco_output=FILE --padmetRef=FILE --output=FILE [-v]
    
    options:
        -h --help     Show help.
        --meneco_output=FILE    pathname of a meneco run' result
        --padmetRef=FILE    path to padmet file corresponding to the database of reference (the repair network)
        --output=FILE    path to tsv output file
"""
from padmet.classes import PadmetRef
from padmet.utils.connection import enhanced_meneco_output
import docopt


def main():
    args = docopt.docopt(__doc__)
    
    meneco_output_file = args["--meneco_output"]
    output = args["--output"]
    verbose = args["-v"]
    padmetRef = PadmetRef(args["--padmetRef"])
    enhanced_meneco_output.enhanced_meneco_output(meneco_output_file, padmetRef, output, verbose)

if __name__ == "__main__":
    main()           