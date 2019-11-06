# -*- coding: utf-8 -*-
"""
Description:
    compare reactions in sbml and padmet file

::
    
    usage:
        compare_sbml_padmet.py --padmet=FILE --sbml=FILE
    
    option:
        -h --help    Show help.
        --padmet=FILE    path of the padmet file
        --sbml=FILE    path of the sbml file
"""
from padmet.classes import PadmetSpec
from padmet.utils.exploration import compare_sbml_padmet
import libsbml
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML(sbml_file)
    for i in range(sbml_document.getNumErrors()):
        print(sbml_document.getError(i).getMessage())

    padmet = PadmetSpec(args["--padmet"])
    compare_sbml_padmet.compare_sbml_padmet(sbml_document, padmet)
    
if __name__ == "__main__":
    main()