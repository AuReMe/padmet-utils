# -*- coding: utf-8 -*-
"""
Description:
    convert a padmet representing a database (padmetRef) and/or a padmet representing a model (padmetSpec)
    to tsv files for askomics.

    1./ Folder creation
    given the output directory. Create this directory if required and create a folder
    padmetRef filename and/or padmetSpec filename
    
    2./

    2.1/ For padmetRef:

        2.1.a/ Nodes
            get all reactions nodes => extract data from misc with extract_nodes(rxn_nodes, "reaction", "../rxn.tsv")

            get all compounds nodes => extract data from misc with extract_nodes(cpd_nodes, "compounds", "../cpd.tsv")

            get all pathways nodes => extract data from misc with extract_nodes(pwy_nodes, "pathway", "../pwy.tsv")

            get all xrefs nodes => extract data from misc with extract_nodes(xref_nodes, "xref", "../xref.tsv")
    
        2.1.b/ Relations
            for each rxn in rxn nodes:

                get all rlt consumes/produces => create list of data with extract_rxn_cpd(rxn_cpd_rlt)
                    fieldnames = "rxn_cpd","concerns@reaction","consumes@compound","produces@compound","stoichiometry","compartment"

                get all rlt is_in_pathway => create list of data with extract_rxn_pwy(rxn_pwy_rlt)
                    fieldnames = "rxn_pwy","concerns@reaction","in_pwy@pathway"

                get all rlt has_xref => create list of data with extract_entity_xref(rxn_xref_rlt)

            for each cpd in cpd nodes:

                get all rlt has_xref => update previous list of data with extract_entity_xref(cpd_xref_rlt)
                    fieldnames = "entity_xref","concerns@reaction","concerns@compound","has_xref@xref"

::
    
    usage:
        padmet_to_tsv.py --padmetSpec=FILE [--padmetRef=FILE] --output_dir=DIR [-v]
        padmet_to_tsv.py --padmetRef=FILE [--padmetSpec=FILE] --output_dir=DIR [-v]
    
    options:
        -h --help     Show help.
        --padmetSpec=FILE    path of the padmet representing the network to convert
        --padmetRef=FILE    path of the padmet representing the database
        --output_dir=DIR
        -v
"""
from padmet.utils.connection import padmet_to_tsv
import docopt

def main():
    args = docopt.docopt(__doc__)
    padmetSpec_file = args["--padmetSpec"]
    padmetRef_file = args["--padmetRef"]
    output_dir = args["--output_dir"]
    verbose = args["-v"]
    padmet_to_tsv.padmet_to_tsv(padmetSpec_file, padmetRef_file, output_dir, verbose)


if __name__ == "__main__":
    main()
