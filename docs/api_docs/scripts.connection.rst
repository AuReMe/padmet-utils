===================
Scripts: Connection
===================
Description:

#TODO

biggAPI\_to\_padmet
======================================

.. automodule:: scripts.connection.biggAPI_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

enhanced\_meneco\_output
======================================

.. automodule:: scripts.connection.enhanced_meneco_output
    :members:
    :undoc-members:
    :show-inheritance:

enhanced\_sgs\_output
======================================

.. automodule:: scripts.connection.enhanced_sgs_output
    :members:
    :undoc-members:
    :show-inheritance:

extract\_rxn\_with\_gene\_assoc
======================================

.. automodule:: scripts.connection.extract_rxn_with_gene_assoc
    :members:
    :undoc-members:
    :show-inheritance:

gbk\_to\_faa
======================================

.. automodule:: scripts.connection.gbk_to_faa
    :members:
    :undoc-members:
    :show-inheritance:

gene\_to\_targets
======================================

.. automodule:: scripts.connection.gene_to_targets
    :members:
    :undoc-members:
    :show-inheritance:

modelSeed\_to\_padmet
======================================

.. automodule:: scripts.connection.modelSeed_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_asp
======================================

.. automodule:: scripts.connection.padmet_to_asp
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_matrix
======================================

.. automodule:: scripts.connection.padmet_to_matrix
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_padmet
======================================

.. automodule:: scripts.connection.padmet_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

padmet\_to\_tsv
======================================

.. automodule:: scripts.connection.padmet_to_tsv
    :members:
    :undoc-members:
    :show-inheritance:

pgdb\_to\_padmet
======================================

.. automodule:: scripts.connection.pgdb_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

post\_pantograph\_gbr
======================================

.. automodule:: scripts.connection.post_pantograph_gbr
    :members:
    :undoc-members:
    :show-inheritance:

pre\_pantograph
======================================

.. automodule:: scripts.connection.pre_pantograph
    :members:
    :undoc-members:
    :show-inheritance:

sbmlGenerator
======================================

.. automodule:: scripts.connection.sbmlGenerator
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

.. automodule:: scripts.connection.sbml_to_padmet
    :members:
    :undoc-members:
    :show-inheritance:

wikiGenerator
======================================

.. automodule:: scripts.connection.wikiGenerator
    :members:
    :undoc-members:
    :show-inheritance:

