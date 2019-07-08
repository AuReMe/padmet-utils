====================
Scripts: Exploration
====================
Description:

#TODO

compare\_padmet
==========================================

.. automodule:: padmet_utils.exploration.compare_padmet
    :members:
    :undoc-members:
    :show-inheritance:

compare\_sbml
========================================

Description:
    compare reactions in two sbml.

    Returns if a reaction is missing

    And if a reaction with the same id is using different species or different reversibility

::

    usage:
        compare_sbml.py --sbml1=FILE --sbml2=FILE
    
    option:
        -h --help    Show help.
        --sbml1=FILE    path of the first sbml file
        --sbml2=FILE    path of the second sbml file

compare\_sbml\_padmet
================================================

.. automodule:: padmet_utils.exploration.compare_sbml_padmet
    :members:
    :undoc-members:
    :show-inheritance:

convert\_sbml\_db
============================================

.. automodule:: padmet_utils.exploration.convert_sbml_db
    :members:
    :undoc-members:
    :show-inheritance:

dendrogram\_reactions\_distance
============================================

Description:
    Use reactions.csv file from compare_padmet.py to create a dendrogram using a Jaccard distance.
    
    From the matrix absence/presence of reactions in different species computes a Jaccard distance between these species.
    Apply a hierarchical clustering on these data with a complete linkage. Then create a dendrogram.
    Apply also intervene to create an upset graph on the data.


::

    usage:
        dendrogram_reactions_distance.py --reactions=FILE --output=FILE [--padmetRef=STR] [--pvclust] [--upset=INT] [-v]
    
    option:
        -h --help    Show help.
        -r --reactions=FILE    pathname of the file containing reactions in each species of the comparison.
        -o --output=FOLDER    path to the output folder.
        --pvclust    launch pvclust dendrogram using R
        --padmetRef=STR    path to the padmet Ref file
        -u --upset=INT    number of cluster in the upset graph.
        -v    verbose mode.

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

.. automodule:: padmet_utils.exploration.get_pwy_from_rxn
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_stats
========================================

.. automodule:: padmet_utils.exploration.get_pwy_from_rxn
    :members:
    :undoc-members:
    :show-inheritance:

report\_network
==========================================

.. automodule:: padmet_utils.exploration.report_network
    :members:
    :undoc-members:
    :show-inheritance:

visu\_path
=====================================

.. automodule:: padmet_utils.exploration.visu_path
    :members:
    :undoc-members:
    :show-inheritance:

