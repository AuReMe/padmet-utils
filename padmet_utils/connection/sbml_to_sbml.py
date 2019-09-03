#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Create sbml from sbml. Use it to change sbml level.

::
    usage:
        sbml_to_sbml.py --input=FILE/FOLDER --output=FILE/FOLDER --new_sbml_lvl=STR [--cpu=INT]
    
    option:
        -h --help    Show help.
        --input=FILE    path of the sbml file/folder to convert into sbml
        --output=FILE    path of the sbml file/folder to generate.
        --new_sbml_lvl=STR    level of the new sbml.
        --cpu=FILE    number of cpu used [default: 1].
        -v   print info.
"""

import docopt
import os
import sys

from multiprocessing import Pool
from padmet_utils.connection.sbmlGenerator import padmet_to_sbml
from padmet_utils.connection.sbml_to_padmet import from_sbml_to_padmet

def main():
    args = docopt.docopt(__doc__)
    input_sbml = args["--input"]
    output_sbml = args["--output"]
    new_sbml_level = args["--new_sbml_lvl"]
    cpu = args['--cpu']

    from_sbml_to_sbml(input_sbml, output_sbml, new_sbml_level, cpu)


def run_sbml_to_sbml(multiprocess_data):
    """Turn sbml to sbml.

    Parameters
    ----------
    multiprocess_data: dict
        pathname to species sbml file, pathname to output sbml file, new sbml level

    Returns
    -------
    bool:
        True if sbml file exists
    """
    padmet = from_sbml_to_padmet(sbml=multiprocess_data['sbml_file'], db=None, version=None,
                                padmetSpec_file=None, source_tool=None, source_category=None, source_id=None, padmetRef_file=None, mapping=None, verbose=None)
    padmet_to_sbml(padmet, multiprocess_data['sbml_output_file'], sbml_lvl=multiprocess_data['new_sbml_level'], verbose=False)

    if multiprocess_data['sbml_output_file'] and not os.access(multiprocess_data['sbml_output_file'], os.W_OK):
        try:
            open(multiprocess_data['sbml_output_file'], 'w').close()
            os.unlink(multiprocess_data['sbml_output_file'])
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True


def from_sbml_to_sbml(input_sbml, output_sbml, new_sbml_level, cpu):
    """Turn sbml to sbml.

    Parameters
    ----------
    input_sbml: str
        pathname to species sbml file/folder
    output_sbml: str
        pathname to output sbml file/folder
    new_sbml_level: int
        new sbml level
    cpu: int
        number of cpu

    Returns
    -------
    str:
        pathname to output sbml file/folder
    """
    try:
        new_sbml_level = int(new_sbml_level)
    except:
        print('SBML level must be an int.')

    if new_sbml_level not in [2, 3]:
        sys.exit('New SBML level must be 2 or 3.')

    sbml_to_sbml_pool = Pool(processes=cpu)

    multiprocess_datas = []
    if os.path.isdir(input_sbml):
        if not os.path.isdir(output_sbml):
            try:
                os.makedirs(output_sbml)
            except OSError:
                print('Impossible to create output folder.')
        for sbml_file in os.listdir(input_sbml):
            multiprocess_data = {}
            multiprocess_data['sbml_file'] = input_sbml + '/' + sbml_file
            multiprocess_data['sbml_output_file'] = output_sbml + '/' + sbml_file
            multiprocess_data['new_sbml_level'] = new_sbml_level
            multiprocess_datas.append(multiprocess_data)

    elif os.path.isfile(input_sbml):
        multiprocess_data = {}
        multiprocess_data['sbml_file'] = input_sbml
        multiprocess_data['sbml_output_file'] = output_sbml
        multiprocess_data['new_sbml_level'] = new_sbml_level
        multiprocess_datas.append(multiprocess_data)

    sbml_checks = sbml_to_sbml_pool.map(run_sbml_to_sbml, multiprocess_datas)

    sbml_to_sbml_pool.close()
    sbml_to_sbml_pool.join()

    if all(sbml_checks):
        return output_sbml

    else:
        print("Error during sbml creation.")

if __name__ == "__main__":
    main()