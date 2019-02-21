====================
Scripts: Exploration
====================
Description:

#TODO

compare\_padmet
==========================================

.. automodule:: scripts.exploration.compare_padmet
    :members:
    :undoc-members:
    :show-inheritance:

compare\_sbml
========================================
Description:
    compare reactions in two sbml

::

    usage:
        compare_sbml.py --sbml1=FILE --sbml2=FILE
    
    option:
        -h --help    Show help.
        --sbml1=FILE    path of the first sbml file
        --sbml2=FILE    path of the second sbml file

compare\_sbml\_padmet
================================================

.. automodule:: scripts.exploration.compare_sbml_padmet
    :members:
    :undoc-members:
    :show-inheritance:

convert\_sbml\_db
============================================

.. automodule:: scripts.exploration.convert_sbml_db
    :members:
    :undoc-members:
    :show-inheritance:

flux\_analysis
=========================================
Description:
    Run flux balance analyse with cobra package. If the flux is >0. Run also FVA
    and return result in standard output

::
    
    usage:
        flux_analysis.py --sbml=FILE
        flux_analysis.py --sbml=FILE --seeds=FILE --targets=FILE
        flux_analysis.py --sbml=FILE --all_species
        
    
    option:
        -h --help    Show help.
        --sbml=FILE    pathname to the sbml file to test for fba and fva.
        --seeds=FILE    pathname to the sbml file containing the seeds (medium).
        --targets=FILE    pathname to the sbml file containing the targets.
        --all_species    allow to make FBA on all the metabolites of the given model.

get\_pwy\_from\_rxn
==============================================

.. automodule:: scripts.exploration.get_pwy_from_rxn
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_stats
========================================
#TODO

report\_network
==========================================

.. automodule:: scripts.exploration.report_network
    :members:
    :undoc-members:
    :show-inheritance:

visu\_path
=====================================

.. automodule:: scripts.exploration.visu_path
    :members:
    :undoc-members:
    :show-inheritance:

