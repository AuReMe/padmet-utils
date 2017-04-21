---
title:  padmet-utils - Documentation
author: Meziane AITE
date: 2017-04-20
version: 2.4
geometry: margin=2cm
---

                  ___                      _
                 |  _`\                   ( )
                 | (_) )   __     _ _    _| |  ___ ___     __
                 | ,  /  /'__`\ /'_` ) /'_` |/' _ ` _ `\ /'__`\
                 | |\ \ (  ___/( (_| |( (_| || ( ) ( ) |(  ___/
                 (_) (_)`\____)`\__,_)`\__,_)(_) (_) (_)`\____)

                        @author: Meziane AITE

################################################################################

##Descirption
The main concept underlying padmet-utils is to provide solutions that ensure the consistency, the internal standardization and the reconciliation of the information used within any workflow that combines several tools involving metabolic networks reconstruction or analysis. The PADMet package is at the core of the AuReMe workflow, dedicated to the primary reconstruction of genome-scale metabolic networks from raw data. It allows the study of organisms for which few experimental data are available. Its main feature is to undergo the reconstruction of the metabolic network by combining several
heterogeneous knowledge and data sources, including the information reported by several scaffold metabolic networks for cousin species.

##Documentation

.aureme.add_seeds_reactions.py
	Description:
	For a given set of compounds representing the growth medium (or seeds). Create 2 reactions
	for each compounds to maintain consistency of the network for flux analysis.
	For each compounds create 
	    -an exchange reaction: this reaction consumes the compound in the
	compartment 'C-BOUNDARY' and produces the compound in the compartment 'e' extracellular
	    -a transport reaction: this reaction consumes the compound in the compartment
	'e' extracellular' and produces the compound in the compartment 'c' cytosol
	ex: for seed 'cpd-a'
	1/ check if cpd-a in padmetSpec, if not, copy from padmetRef.
	2/ create exchange reaction: ExchangeSeed_cpd-a_b: 1 cpd-a (C-BOUNDARAY) => 1 cpd-a (e)
	3/ create transport reaction: TransportSeed_cpd-a_e: 1 cpd-a (e) => 1 cpd-a (c)
	4/ create a new file if output not None, or overwritte padmetSpec

	usage:
	    add_seeds_reactions.py --padmetSpec=FILE --padmetRef=FILE --seeds=FILE [--output=FILE] [-v]

	options:
	    -h --help     Show help.
	    --padmetSpec=FILE    pathanme to the padmet file to update
	    --padmetRef=FILE    pathname to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
	    --seeds=FILE    pathname to the file of seeds ids. 1/line
	    --output=FILE    If not None, pathname to the padmet file updated
	    -v   print info

.aureme.compounds_to_sbml.py
	Convert a list of compounds to an sbml file. Often in sbml, the compartment of the
	compounds is concatenate to the id. If wanted, possible to add a compart each compounds with compart_name
	If verbose, will also check if each compound is in a padmetRef and or a padmetSpec
	In sbml: id = id encoded, name = original id
	ex: id = cpd-a, id encoded = M_cpd__45__a_'value of compart_name'
	if verbose True, will print if each compound was found in padmetRef and or padmetSpec

	usage:
	    compounds_to_sbml.py --compounds=FILE --output=FILE [--padmetRef=FILE] [--padmetSpec=FILE] [--compart_name=STR] [-v]

	option:
	    -h --help     Show help.
	    --compounds=FILE    pathname of the file containing the compound's id. one by line.
	    --output=FILE    pathname to the sbml output.
	    --padmetRef=FILE   pathname of the padmet used as reference (ex: metacyc_18.5.padmet)
	    --padmetSpec=FILE   pathname of the padmet corresponding to the species studied
	    --compart_name=str   compart id to concatenate to compound id
	    -v    print info and check if ids in padmetRef and or padmet

.aureme.enhanced_meneco_output.py
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

	usage:
	    enhanced_meneco_output.py --meneco_output=FILE --padmetRef=FILE --output=FILE [-v]

	options:
	    -h --help     Show help.
	    --meneco_output=FILE    pathname of a meneco run' result
	    --padmetRef=FILE    pathanme to padmet file corresponding to the database of reference (the repair network)
	    --output=FILE    pathname to tsv file containing more informations about the reactions.

.aureme.extract_deadends.py
	Description:
	extract first solution of dead_ends on output of clingo deadend_encoding.lp

	usage:
	    extract_deadends.py --clasp_output=FILE --output=FILE

	option:
	    -h --help    Show help.
	    --clasp_output=FILE
	    --output=FILE    file with deadends id.

.aureme.extract_rxn_with_gene_assoc.py
	Description:
	From a given sbml file, create a sbml with only the reactions associated to a gene.
	Need for a reaction, in section 'note', 'GENE_ASSOCIATION': ....

	usage:
	    extract_rxn_with_gene_assoc.py --sbml=FILE --output=FILE [-v]

	options:
	    -h --help     Show help.
	    --sbml=FILE    pathname to the sbml file
	    --output=FILE    pathname to the sbml output (with only rxn with genes assoc)
	    -v   print info

.aureme.fba_test.py
	Description:
	Run flux balance analyse with cobra package. If the flux is >0. Run also FVA
	and return result in standard output

	usage:
	    fba_test.py --sbml=FILE

	option:
	    -h --help    Show help.
	    --sbml=FILE    pathname to the sbml file to test for fba and fva.

.aureme.gbk_to_faa.py
	Description:
	convert GBK to FAA with Bio package

	usage:
	    gbk_to_faa.py    --gbk=FILE --output=FILE [--qualifier=STR] [-v]

	option:
	    -h --help    Show help.
	    --gbk=FILE    pathname to the gbk file.
	    --output=FILE    pathename to the output, a FAA file.
	    --qualifier=STR    the qualifier of the gene id [default: locus_tag].
	    -v   print info

.aureme.gene_to_targets.py
	Description:
	From a list of genes, recovere from the linked reactions the list of products.
	R1 is linked to G1, R1 produces M1 and M2.  output: M1,M2. Takes into account reversibility

	usage:
	    gene_to_targets.py --padmetSpec=FILE --genes=FILE --output=FILE [-v]

	option:
	    -h --help     Show help
	    --padmetSpec=FILE    pathanme to the padmet file
	    --genes=FILE   pathname to the file containing gene ids, one id by line
	    --output=FILE    pathname to the output file containing all tagerts which can by produced by all reactions associated to the given genes
	    -v   print info

.aureme.manual_curation.py
	Description:
	This script allows to combine reaction_creator.py and update_padmetSpec.py.
	This script was created specially for AuReMe and the default metabolic network
	reconstruction workflow.
	If the file reaction_to_add_delete exist: calling update_padmetSpec
	If the file new_reaction_data exist: calling reaction_creator.py
	Update padmetSpec and creating a new padmet (new_padmet) or overwritte the input

	usage:
	    manual_curation.py --padmetSpec=FILE --padmetRef=FILE --reaction_to_add_delete=FILE  --new_reaction_data=FILE [--new_padmet=FILE] [-v]

	option:
	    -h --help    Show help.
	    --padmetSpec=FILE    pathname to the padmet to update
	    --padmetRef=FILE    pathname of the padmet representing the reference database
	    --reaction_to_add_delete=FILE    pathname to the file used for update_padmetSpec.py
	    --reaction_creator=FILE    pathname to the file used for reaction_creator.py
	    --new_padmet=FILE    pathname to the ouput. if None. Overwritting padmetSpec
	    -v    print info

.aureme.pre_pantograph.py
	Description:
	Before running pantograph it is necessary to check if the metabolic network 
	and the proteom of the model organism use the same ids for genes (or at least more than a given cutoff). 
	To only check this. Use the 2nd usage.
	If the genes ids are not the same, it is necessary to use a dictionnary of genes ids associating
	the genes ids from the proteom to the genes ids from the metabolic network.
	To create the correct proteom from the dictionnnary, use the 3nd usage
	Finnaly by using the 1st usage, it is possible to:
	    1/ Check model_faa and model_metabolic for a given cutoff
	    2/ if under the cutoff, convert model_faa to the correct one with dict_ids_file
	    3/ if still under, SystemExit()

	usage:
	    pre_pantograph.py    --model_metabolic=FILE    --model_faa=FILE    [--cutoff=FLOAT] [--dict_ids_file=FILE] --output=FILE    [-v]
	    pre_pantograph.py    --model_metabolic=FILE    --model_faa=FILE    [--cutoff=FLOAT] [-v]
	    pre_pantograph.py    --model_faa=FILE    --dict_ids_file=FILE    --output=FILE [-v]

	option:
	    -h --help    Show help.
	    --model_metabolic=FILE    pathname to the metabolic network of the model (sbml).
	    --model_faa=FILE    pathname to the proteom of the model (faa)
	    --cutoff=FLOAT    cutoff [0:1] for comparing model_metabolic and model_faa. [default: 0.70]. 
	    --dict_ids_file=FILE    pathname to the dict associating genes ids from the model_metabolic to the model_faa. line = 
	    --output=FILE    output of get_valid_faa (a faa) or get_dict_ids (a dictionnary of gene ids in tsv)
	    -v   print info

.aureme.reactions_to_sbml.py
	Description:
	Convert a list of reactions to SBML file.

	usage:
	    reactions_to_sbml.py --padmetRef=FILE --reactions=FILE --output=FILE [-v]

	option:
	    -h --help     Show help.
	    --padmetRef=FILE   pathname to the padmet used as reference (ex: metacyc.padmet)   
	    --reactions=FILE    pathname to the file of reaction id. one id by line.
	    --output=FILE    pathname to the sbml output.
	    -v    print info

.misc.bigg_to_padmet.py
	Description:
	Create a padmet file from 2 tsv files (reactions.tsv and metabolites.tsv)
	metabolites.tsv: col1 = metabolite ID, col2 = metabolite Name 
	ex: 10fthf_c 10-Formyltetrahydrofolate
	reactions.tsv: col1 = reaction ID, col2 = reaction Name, col3 = formula, col4 = reversibility
	ex: 10FTHF5GLUtl 10FTHF5GLUtl 1.0:10fthf5glu_c => 1.0:10fthf5glu_l False
	NB: the formula must be in this format: STOICHIOMETRY:METABOLITE_ID + ... => STOICHIOMETRY:METABOLITE_ID

	for reaction in reactions:
	    /!\ reactio id is modified with sbmlPlugin.converto_from_coded_id. delete the type information (R_)
	    create Node, class = reaction, id = id, misc.COMMON_NAME = name, misc.DIRECTION = reversible ('True' = 'REVERSIBLE', 'False' = 'LEFT-TO-RIGHT')
	    Parse formula to recovere listOfReactants, listOfProducts
	    for metabolite in listOfReactants:
		extract compart with convert_from_coded_id. reactant_id = id without compart (if compart not None)
		create relation, type = consumes, idIn = reaction_id, idOut = reactant_id, relation.misc["COMPARTMENT"] = [compart], relation.misc["STOICHIOMETRY"] = [stoic] 
	    for metabolite in listOfProducts:
		extract compart with convert_from_coded_id. product_id = id without compart (if compart not None)
		create relation, type = produces, idIn = reaction_id, idOut = product_id, relation.misc["COMPARTMENT"] = [compart/Unknown], relation.misc["STOICHIOMETRY"] = [stoic] 

	for metabolite in metabolites.tsv:
	    /!\ metabolite id is modified. => delete the type information (M_) and the compart info (_c/e/p..)
	    extract compart with convert_from_coded_id. product_id = id without compart (if compart not None)
	    add in node.misc["COMMON_NAME"] = metaboName

	usage:
	    bigg_to_padmet.py --version=STR --output=FILE --rxn_file=FILE [--met_file=FILE] [-v]

	options:
	    -h --help     Show help.
	    --version=STR    The tag associated to the origin of the information used
	    --output=FILE    pathname of the padmet file to create
	    --rxn_file=FILE    pathanme of the reactions file. col1 = reaction ID, col2 = reaction Name, col3 = formula, col4 = reversibility.
	    --met_file=FILE    pathanme of the metabolites file. col1 = metabolite ID, col2 = metabolite Name 
	    -v   print info

.misc.change_compart.py
	Description:
	Two disctinct usages:
	    change_compart.py --padmetSpec=FILE --compart=STR --output=FILE [-v]
	For a given compartment to check return a file with all metabolites in it
	Then, to change the compartment of one or more metabolites, modify the precedent output
	by changing the compartment.
	ex: output line: M_icit_x	x
	M_icit_x is in compart x. to change the compart to 'c' for example, edit 'x' to 'c'.
	Finally use the second usages to update the model.
	    change_compart.py --padmetSpec=FILE --updateFile=FILE [--new_padmet=FILE] [-v]
	NB: the metabolite id also change: M_icit_x => M_icit_c

	usage:
	    change_compart.py --padmetSpec=FILE --compart=STR --output=FILE [-v]
	    change_compart.py --padmetSpec=FILE --updateFile=FILE [--new_padmet=FILE] [-v]

	options:
	    -h --help     Show help.
	    --padmetSpec=FILE    pathname of the padmet file to update from updateFile
	    --compart=STR    the compartment id to check
	    --output=FILE    pathname of the tsv file with col 1 = compounds ids, col 2 = the current compartment
	    --updateFile=FILE    pathname of the output after modifications of col 2
	    --new_padmet=FILE   pathname of the new padmet. if None, will erase the current
	    -v   print info

.misc.compare_sbml_padmet.py
	Description:
	compare reactions in sbml and padmet files

	usage:
	    compare_sbml_padmet.py --padmet=FILE --sbml=FILE

	option:
	    -h --help    Show help.
	    --padmet=FILE    pathname of the padmet file
	    --sbml=FILE    pathanme of the sbml file

.misc.padmetSpec_to_ASP_for_deadends.py
	Description:
	use the script padmetSpec_to_asp_for_deadend to convert a padmet file to .lp for
	asp. the output is then used with clingo and deadend_encoding.lp to find deadends
	usage:
	    padmetSpec_to_ASP_for_deadends.py --padmetSpec=FILE --seeds=FILE [--targets=FILE] --output=FILE [-v]

	option:
	    -h --help    Show help.
	    --padmetSpec_file=FILE    padmetSpec file to convert into .lp for asp.
	    --seeds_file=FILE
	    --targets_file=FILE
	    --output=FILE    lp file to create.
	    -v   print info.

.misc.padmet_to_asp.py
	Description:
	Convert PADMet to ASP following this predicats:
	direction(reaction_id, reaction_direction). #LEFT-TO-RIGHT or REVERSIBLE
	ec_number(reaction_id, ec(x,x,x)).
	catalysed_by(reaction_id, enzyme_id).
	uniprotID(enzyme_id, uniprot_id). #if has has_xref and db = "UNIPROT"
	common_name({reaction_id or enzyme_id or pathway_id or compound_id} , common_name)

	usage:
	    padmet_to_ASP.py --padmet=FILE --output=FILE [-v]

	options:
	    -h --help     Show help.
	    --padmet=FILE    pathname of the padmet file to convert
	    --output=FILE    pathname of the lp file
	    -v   print info

.misc.padmet_to_sbml.py
	Description:
	Use the script padmet_to_sbml to convert padmet to a sbml file
	give only the id of the reaction to test (obj_coef)

	usage:
	    padmet_to_sbml.py --padmet=FILE --output=FILE [--obj_fct=ID] [--sbml_lvl=x] [--sbml_version=y] [-v]

	option:
	    -h --help    Show help.
	    --padmet=FILE    pathname of the padmet file to convert into sbml
	    --output=FILE    pathanme of the sbml file to generate.
	    --obj_fct=ID    id of the reaction objective.
	    --sbml_lvl=x    sbml level 2 is sufficient for FBA [default: 2].
	    --sbml_version=y    sbml version 1 is sufficient for FBA [default: 1].
	    -v   print info.

.misc.pgdb_to_padmet.py
	Description:

	classes.dat:
	For each class:
	create new node / class = class
	UNIQUE-ID (1) => node._id = UNIQUE-ID
	COMMON-NAME (0-n) => node.Misc['COMMON_NAME'] = COMMON-NAME
	TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
	SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name synonyms)
	 
	compounds.dat:
	for each compound:
	create new node / class = compound
	UNIQUE-ID (1) => node._id = UNIQUE-ID
	COMMON-NAME (0-n) => node.Misc['COMMON_NAME'] = COMMON-NAME
	INCHI-KEY (0-1) {InChIKey=XXX} => node.misc['INCHI_KEY': XXX]
	MOLECULAR-WEIGHT (0-1) => node.misc()['MOLECULAR_WEIGHT'] = MOLECULAR-WEIGHT
	SMILES (0-1) => node.misc()['SMILES'] = SMILES
	TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
	SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name name)
	DBLINKS (0-n) {(db "id" ...)} => for each, create new node xref, create rlt (node has_xref xref)

	proteins.dat:
	for each protein:
	create new node / class = protein
	UNIQUE-ID (1) => node._id = UNIQUE-ID
	COMMON-NAME (0-n) => node.Misc['COMMON_NAME'] = COMMON-NAME
	INCHI-KEY (0-1) {InChIKey=XXX} => node.misc['INCHI_KEY': XXX]
	MOLECULAR-WEIGHT (0-1) => node.misc()['MOLECULAR_WEIGHT'] = MOLECULAR-WEIGHT
	SMILES (0-1) => node.misc()['SMILES'] = SMILES
	TYPES (0-n) => for each, check or create new node class, create rlt (node is_a_class types)
	SYNONYMS (0-n) => for each, create new node name, create rlt (node has_name name)
	DBLINKS (0-n) {(db "id" ...)} => for each, create new node xref, create rlt (node has_xref xref)
	SPECIES (0-1) => for each, check or create new node class, create rlt (node is_in_species class)

	reactions.dat:
	for each reaction:
	create new node / class = reaction + node.misc()["DIRECTION"] = "UNKNOWN" by default
	UNIQUE-ID (1) => node._id = UNIQUE-ID
	COMMON-NAME (0-n) => node.Misc['COMMON_NAME'] = COMMON-NAME
	EC-NUMBER (0-n) => node.Misc['EC-NUMBER'] = EC-NUMBER
	REACTION-DIRECTION (0-1) => node.Misc['DIRECTION'] = reaction-direction, if REVERSIBLE, else: LEFT-TO-RIGHT
	RXN-LOCATIONS (0,n) => node.misc['COMPARTMENT'] = rxn-location
	TYPES (0-n) => check or create new node class, create rlt (node.id is_a_class types's_node.id)
	DBLINKS (0-n) {(db "id" ...)} => create new node xref, create rlt (node has_xref xref's_node.id)
	SYNONYMS (0-n) => create new node name, create rlt (node has_name name's_node.id)
	--
	for LEFT and RIGHT, also check 2 next lines if info about 'coefficient' or 'compartment'
	defaut value: coefficient/stoichiometry = 1, compartment = unknown
	also check if the direction is 'RIGHT-TO-LEFT', if yes, inverse consumes and produces relations
	then change direction to 'LEFT-TO-RIGHT'
	LEFT (1-n) => create rlt (node.id consumes left's_node.id)
	RIGHT (1-n) => create rlt (node.id produces right's_node.id)

	enzrxns.dat:
	for each association enzyme/reaction:
	create new rlt / type = catalyses
	ENZYME (1) => stock enzyme as 'enzyme catalyses'
	REACTION (1-n) => for each reaction after, create relation 'enzyme catalyses reaction'

	pathways.dat:
	for each pathway:
	create new node / class = pathway
	UNIQUE-ID (1) => node._id = UNIQUE-ID
	TYPES (0-n) => check or create new node class, create rlt (node is_a_class types)
	COMMON-NAME (0-n) => node.Misc['COMMON_NAME'] = COMMON-NAME
	DBLINKS (0-n) {(db "id" ...)} => create new node xref, create rlt (node has_xref xref)
	SYNONYMS (0-n) => create new node name, create rlt (node has_name name)
	IN-PATHWAY (0-n) => check or create new node pathway, create rlt (node is_in_pathway name)
	REACTION-LIST (0-n) => check or create new node pathway, create rlt (node is_in_pathway name)

	metabolic-reaction.xml: optional
	for each reaction:

	usage:
	    pgdb_to_padmet.py --version=V --db=ID --output=FILE --directory=DIR [-v] [-g] [-m]
	    pgdb_to_padmet.py --version=V --db=ID --output=FILE --classes_file=FILE --compounds_file=FILE --proteins_file=FILE --reactions_file=FILE --enzrxns_file=FILE --pathways_file=FILE [--genes_file=FILE] [--metabolic_reactions=FILE] [-v]

	options:
	    -h --help     Show help.
	    --version=V    metacyc version
	    --output=FILE    padmet file corresponding to the DB
	    --directory=DIR    directory containg all the .dat files of metacyc (data)
	    --classes_file=FILE   
	    --compounds_file=FILE   
	    --proteins_file=FILE 
	    --reactions_file=FILE 
	    --enzrxns_file=FILE 
	    --pathways_file=FILE
	    --genes_file=FILE
	    --metabolic_reactions=FILE
	    -m    use the metabolic-reactions.xml file to enhance the database
	    -g    use the genes_file (use if its a specie's pgdb, if metacyc, do not use)
	    -v   print info

.misc.reaction_creator.py
	Description:
	Allows to create news reactions that are not in the database of reference and add
	them in a padmet file. First, fill the template (to get the template use the 2nd usage).
	Then, use the 1st usage to add the reaction(s) in new_reaction_data to padmetSpec.
	The script will check if the choosen id are available and if the compounds used exist or not.
	If 1-n compounds don't exist, they will be created. Checking the existence of compounds
	in padmetSpec AND padmetRef. 

	usage:
	    reaction_creator.py --padmetSpec=FILE --padmetRef=FILE --new_reaction_data=FILE [--new_padmet=FILE] [-v]
	    reaction_creator.py --getTemplate --output=FILE

	options:
	    -h --help     Show help.
	    --padmetSpec=FILE    pathname of the padmet file to update
	    --padmetRef=FILE    pathname of the padmet file representing the database of reference
	    --reaction_creator=FILE    pathname to the file containing the data (template obtained from 2nd usage)
	    --new_padmet=FILE    new name of padmet after update
	    -v   print info

.misc.report_network.py
	Description:
	Create reports of a padmet file. 
	all_pathways.tsv: header = ["dbRef_id", "Common name", "Number of reaction found", 
		    "Total number of reaction", "Ratio (Reaction found / Total)"]
	all_reactions.tsv: header = ["dbRef_id", "Common name", "formula (with id)", 
		    "formula (with common name)", "in pathways", "associated genes"]
	all_metabolites.tsv: header = ["dbRef_id", "Common name", "Produced (p), Consumed (c), Both (cp)"]

	usage:
	    report_network.py --padmetSpec=FILE --padmetRef=FILE --output_dir=dir [-v]

	options:
	    -h --help     Show help.
	    --padmetSpec=FILE    pathname of the padmet file.
	    --padmetRef=FILE    pathname of the padmet file used as database
	    --output_dir=dir    directory for the results.
	    -v   print info.

.misc.sbml_to_padmet.py
	Description:
	To convert a sbml to padmet without a padmet of database reference use the 2nd usage.
	It is possible to define a specific policy and info for the padmet. To learn more about
	policy and info check doc of lib.padmetRef/Spec.
	To convert a sbml to padmet based on a padmet of database reference use the 1st usage.
	If it is to update a given padmetSpec, then specify the pathname to --padmetSpec.
	If it is to initialize a padmetSpec then directly give the output to create.
	if in the ids of reactions are not the same from padmetRef and the sbml, it is possible to use
	a dictionnary of association (sbml_id padmetRef_id)
	The association id-sbml id-padmetRef can be recovered from a file (assocIdOriginRef)
	with one line = 'id_sbml \t id_padmetRef'
	Finally if a reaction from sbml is not in padmetRef, it is possible to force the copy and creating
	a new reaction in padmetSpec with the arg -f

	usage:
	    sbml_to_padmet.py --padmetRef=FILE --sbml=FILE [--padmetSpec=FILE] [--output=FILE] [--assocIdOriginRef=FILE] [-v] [-f]
	    sbml_to_padmet.py --sbml=FILE --output=FILE [--policy=LIST] [--info=DICT] [-v]

	options:
	    -h --help     Show help.
	    --padmetSpec=FILE    pathname to the padmet file to update with the sbml. If there's no padmetSpec, just specify the output
	    --padmetRef=FILE    pathanme to the padmet file representing to the database of reference (ex: metacyc_18.5.padmet)
	    --sbml=FILE    1 sbml file to convert into padmetSpec (ex: my_network.xml/sbml) OR a directory with n SBML
	    --output=FILE   pathanme to the new padmet file
	    --assocIdOriginRef=FILE    dictionnary of association id_origin id_ref
	    -f   allows to create reactions not found in the padmetRef (force)
	    -v   print info

.misc.update_padmetSpec.py
	Description:
	Allows to add or delete reactions/metabolites/pathways based on the given id.
	From the updateFile (obtained with 2nd usage) the script extract the id (col1),
	why adding/deleting this element (comment col2) and the action (col3)
	THE ACTION IS IN: ['add','delete','']
	When add: the element is added from padmetRef, and the COMMENT and SOURCE is stored in the node of
	the element in padmetSpec.
	When delete: the element is deleted from padmetSpec.
	Whene '': Nothing to do, but allow to manage the output of enhancedMenecoOutput.py

	usage:
	    update_padmetSpec.py --padmetSpec=FILE --padmetRef=FILE --updateFile=FILE [--source=STR] [--new_padmet=FILE] [-v]
	    update_padmetSpec.py --getTemplate --output=FILE

	option:
	    -h --help     Show help.
	    --padmetSpec=FILE    pathanme of the padmet file to update
	    --padmetRef=FILE  pathname of the padmet used as database of reactions
	    --updateFile=FILE    pathname of the file containing elements ids to add or delete
	    --new_padmet=FILE    pathanme of the new padmetSpec after update [default: --padmetSpec]
	    --source=STR   source of the reactions to add [default: manual]
	    -v   print info

.misc.visu_path.py
	Description:
	Allows to visualize a pathway in padmet network.
	Color code: 
	reactions associated to the pathway, present in the network: lightGreen
	reactions associated to the pathway, not present in the network: red
	compounds: skyblue

	usage:
	    visu_path.py --padmetSpec=FILE --padmetRef=FILE --pathway=ID

	options:
	    -h --help     Show help.
	    --padmetSpec=FILE    pathname to the PADMet file of the network.
	    --padmetRef=FILE    pathname to the PADMet file of the db of reference.
	    --pathway=ID    pathway id to visualize.

.misc.wikipage_creation.py
	Description:
	Create the files containing the wikicode for each page of the wiki relative to a
	padmetSpec based on a database padmetRef. All the pages are stored in the folder
	output_dir. The padmetRef is required to extract the information about all the reactions
	of a given pathway. In the case of BIGG db, it's possible to use the padmetSpec as
	padmetRef because there is no pathways.

	usage:
	    wikiPage_creation.py --padmetRef=FILE --padmetSpec=FILE --output_dir=DIR [-v]

	option:
	    -h --help    Show help.
	    --padmetRef=FILE    pathname of the padmet of reference.
	    --padmetSpec=FILE    pathname of the padmet of interest.
	    --output_dir=DIR    pathname of the directory were to stock all the wikiPages.
	    -v   print info

III./ Wiki management:



