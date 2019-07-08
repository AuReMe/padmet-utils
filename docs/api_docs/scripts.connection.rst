===================
Scripts: Connection
===================
Description:

#TODO

biggAPI\_to\_padmet
======================================

.. automodule:: padmet_utils.connection.biggAPI_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

check\_orthology\_input
======================================

.. automodule:: padmet_utils.connection.check_orthology_input
    :members:
    :undoc-members:
    :show-inheritance:

enhanced\_meneco\_output
======================================

.. automodule:: padmet_utils.connection.enhanced_meneco_output
    :members:
    :undoc-members:
    :show-inheritance:

extract\_orthofinder
======================================

.. automodule:: padmet_utils.connection.extract_orthofinder
    :members:
    :undoc-members:
    :show-inheritance:

extract\_rxn\_with\_gene\_assoc
======================================

.. automodule:: padmet_utils.connection.extract_rxn_with_gene_assoc
    :members:
    :undoc-members:
    :show-inheritance:

gbk\_to\_faa
======================================

.. automodule:: padmet_utils.connection.gbk_to_faa
    :members:
    :undoc-members:
    :show-inheritance:

gene\_to\_targets
======================================

.. automodule:: padmet_utils.connection.gene_to_targets
    :members:
    :undoc-members:
    :show-inheritance:

modelSeed\_to\_padmet
======================================

.. automodule:: padmet_utils.connection.modelSeed_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_asp
======================================

.. automodule:: padmet_utils.connection.padmet_to_asp
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_matrix
======================================

.. automodule:: padmet_utils.connection.padmet_to_matrix
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_padmet
======================================

.. automodule:: padmet_utils.connection.padmet_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_tsv
======================================

.. automodule:: padmet_utils.connection.padmet_to_tsv
    :members:
    :undoc-members:
    :show-inheritance:

pgdb\_to\_padmet
======================================

.. automodule:: padmet_utils.connection.pgdb_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

sbmlGenerator
======================================

.. automodule:: padmet_utils.connection.sbmlGenerator
    :members:
    :undoc-members:
    :show-inheritance:

sbml\_to\_curation\_form
======================================
Description:
    extract all reactions from a sbml file to the form used in aureme for curation.

::

    usage:
        sbml_to_curation_form.py --sbml=FILE --output=FILE --comment=STR [--rxn_id=ID]
    
    options:
        -h --help     Show help.
        --sbml=FILE    path of the sbml.
        --output=FILE    form containing the reaction extracted, form used for manual curation in aureme.
        --rxn_id=FILE    id of one reaction or n reactions sep by ';', if None try to extract the reaction with objective coefficient == 1.
        --comment=STR    comment associated to the reactions in the form. Used to track sources of curation in aureme.

sbml\_to\_padmet
======================================

.. automodule:: padmet_utils.connection.sbml_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

wikiGenerator
======================================

.. automodule:: padmet_utils.connection.wikiGenerator
    :members:
    :undoc-members:
    :show-inheritance:

