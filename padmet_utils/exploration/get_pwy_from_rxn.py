# -*- coding: utf-8 -*-
"""
Description:
    From a file containing a list of reaction, return the pathways where these reactions 
    are involved.
    ex: if rxn-a in pwy-x => return, pwy-x; all rxn ids in pwy-x; all rxn ids in pwy-x FROM the list; ratio

::
    
    usage:
        get_pwy_from_rxn.py --reaction_file=FILE --padmetRef=FILE  --output=FILE
    
    options:
        -h --help     Show help.
        --reaction_file=FILE    pathname of the file containing the reactions id, 1/line 
        --padmetRef=FILE    pathname of the padmet representing the database.
        --output=FILE    pathname of the file with line = pathway id, all reactions id, reactions ids from reaction file, ratio. sep = "\t"
"""
from padmet.classes import PadmetSpec
from padmet.utils.exploration import get_pwy_from_rxn
import docopt

def main():
    args = docopt.docopt(__doc__)
    reaction_file = args["--reaction_file"]
    with open(reaction_file, 'r') as f:
        reactions = set(f.read().splitlines())
    padmet_file = args["--padmetRef"]
    output = args["--output"]
    padmet = PadmetSpec(padmet_file)
    dict_pwy = get_pwy_from_rxn.extract_pwys(padmet, reactions)
    get_pwy_from_rxn.dict_pwys_to_file(dict_pwy, output)

if __name__ == "__main__":
    main()    
