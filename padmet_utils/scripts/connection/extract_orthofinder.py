# -*- coding: utf-8 -*-
"""
Description:
    #todo
    
::
   
    usage:
        extract_orthofinder --orthology_based=DIR --orthgroups=FILE --study_id=STR --output=DIR [--mnx_rxn=FILE] [--mnx_cpd=FILE] [-v]

    option:
        -h --help    Show help.
        -v   print info
"""
import re
import csv
import libsbml
import docopt
from padmet.utils import sbmlPlugin as sp
import os
import gbr.gbr as gbr_compile


def main():
    args = docopt.docopt(__doc__)
    verbose = args["-v"]
    orthology_based_folder = args["--orthology_based"]
    orthogroups_file = args["--orthogroups"]
    output_folder = args["--output"]
    study_id = args["--study_id"]

    orthofinder_to_sbml(orthogroups_file, orthology_based_folder, output_folder, study_id, verbose)
 
    
def orthofinder_to_sbml(orthogroups_file, orthology_based_folder, output_folder, study_id, verbose = False):
    if verbose:
        print("Parsing Orthofinder output %s" %orthogroups_file)
    #k=orthogroup_id, v = {k = name, v = set of genes}
    dict_orthogroup = {}
    all_orgs = set()
    with open(orthogroups_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter = "\t")
        for row in reader:
            orth_id = row['']
            row.pop('')
            new_dict = dict([(name, set(genes.split(","))) for (name, genes) in list(row.items()) if genes])
            dict_orthogroup[orth_id] = new_dict
            all_orgs.update(new_dict.keys())
    #check all sbml
    all_model_sbml = {}
    for path, _, files in os.walk(orthology_based_folder):
        fold = os.path.basename(path)
        if fold in all_orgs:
            for _file in files:
                if _file == fold+".sbml":
                    all_model_sbml[fold] = os.path.join(path, _file)
    if verbose:
        print("Start sbml creation...")
    all_dict_data = []
    if not os.path.exists(output_folder):
        if verbose:
            print("\tCreating folder %s" %output_folder)
        os.makedirs(output_folder)
    all_models = all_orgs - set([study_id])
    for model_id in all_models:
        output_name = "output_orthofinder_from_{0}.sbml".format(model_id)
        output = os.path.join(output_folder, output_name)
        sbml_template = all_model_sbml[model_id]
        dict_data = {'study_id': study_id, 'model_id': model_id, 'sbml_template': sbml_template, 'output': output, 'verbose':verbose}
        all_dict_data.append(dict_data)

    for dict_data in all_dict_data:
        dict_data_to_sbml(dict_data, dict_orthogroup)


def dict_data_to_sbml(dict_data, dict_orthogroup):
    """
    dict_orthogroup: global var
    """
    #dict_data = {'study_name':'', 'o_compare_name': '', sbml_template':'', 'output':''}
    study_id = dict_data['study_id']
    model_id = dict_data['model_id']
    sbml_template = dict_data['sbml_template']
    output = dict_data['output']
    verbose = dict_data.get('verbose')
    if os.path.isfile(output):
        if verbose:
            print("*{0} is already created, skip".format(os.path.basename(output)))
        return
    if verbose:
        print("*Extracting orthology data to create sbml of {0} from {1}".format(study_id, model_id))

    #k = gene_id from to_compare, v = list of genes id of study
    sub_dict_orth = {}
    for k in list(dict_orthogroup.values()):
        try:
            all_to_compare_genes = k[model_id]
            all_study_genes = k[study_id]
            for to_compare_gene in all_to_compare_genes:
                try:
                    sub_dict_orth[to_compare_gene].update(all_study_genes)
                except KeyError:
                    sub_dict_orth[to_compare_gene] = set(all_study_genes)
        except KeyError:
            pass
                    
    if not sub_dict_orth:
        if verbose:
            print("\t{0} and {1} don't share any ortholgue".format(study_id, model_id))
        return

    reader = libsbml.SBMLReader()
    document_to_compare = reader.readSBML(sbml_template)
    for i in range(document_to_compare.getNumErrors()):
        print(document_to_compare.getError(i).getMessage())
    model_to_compare = document_to_compare.getModel()
    listOfReactions_with_genes = [rxn for rxn in model_to_compare.getListOfReactions()
                                  if sp.parseNotes(rxn).get("GENE_ASSOCIATION",[None])[0]]
    if verbose:
        print("\tSbml of {0} contains {1}/{2} reactions with genes assocation".format(model_id, len(listOfReactions_with_genes), len(model_to_compare.getListOfReactions())))
    dict_rxn_ga = {}
    for rxn in listOfReactions_with_genes:
        ga = sp.parseNotes(rxn)['GENE_ASSOCIATION'][0]
        ga_for_gbr = re.sub(r" or " , "|", ga)
        ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
        ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
        #ga_subsets = eval(subprocess.check_output("python3 grammar-boolean-rapsody.py "+ga_for_gbr, shell=True))
        if re.findall("\||&", ga_for_gbr):
            to_compare_ga_subsets = list(gbr_compile.compile_input(ga_for_gbr))
        else:
            to_compare_ga_subsets = ga_for_gbr
        study_ga_subsets = []
        """
        to_compare_ga_subsets = [('a','c','d'),('c',)]
        sub_dict_orth = {'a':['a_a'],'c':['c_c'], 'd':['d_d']}
        """
        for to_compare_subset in to_compare_ga_subsets:
            study_subset = set()
            for gene in to_compare_subset:
                if gene in list(sub_dict_orth.keys()):
                    study_subset.update(sub_dict_orth[gene])
                else:
                    study_subset = set()
                    break
            if study_subset:
                """
                if verbose:
                    print("\t\t{0} == {1}".format(tuple(to_compare_subset), tuple(study_subset)))
                """
                study_ga_subsets.append(study_subset)
        if study_ga_subsets:
            study_ga = " or ".join(["("+" and ".join(subset)+")" for subset in study_ga_subsets])
            if verbose:
                print("\t\tAdding %s" %rxn.id)
                print("\t\tGENE_ASSOCIATION: %s" %(study_ga))
            dict_rxn_ga[rxn.id] = study_ga
    if not dict_rxn_ga:
        if verbose:
            print("\tNo reaction added from {0} to {1} because of missing orthologues".format(model_id, study_id))
        return
    rxn_id_to_remove = set([rxn.id for rxn in model_to_compare.getListOfReactions()]).difference(list(dict_rxn_ga.keys()))
    if verbose:
        print("\tRemoving %s unused reactions" %len(rxn_id_to_remove))
    [model_to_compare.removeReaction(rxn_id) for rxn_id in rxn_id_to_remove]
    cpd_id_to_preserve = set()
    for rxn_id, study_ga in list(dict_rxn_ga.items()):
        rxn = model_to_compare.getElementBySId(rxn_id)
        #update notes
        notes_in_dict = sp.parseNotes(rxn)
        notes_in_dict["GENE_ASSOCIATION"] = [study_ga]
        notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
        for k,v_list in list(notes_in_dict.items()):
            for v in v_list:
                notes += "<p>"+k+": "+v+"</p>"
        notes += "</body>"
        rxn.setNotes(notes)
        cpd_in_rxn = set([p.getSpecies() for p in rxn.getListOfProducts()]).union(\
                         set([r.getSpecies() for r in rxn.getListOfReactants()]))
        cpd_id_to_preserve.update(cpd_in_rxn)
    all_species = [cpd.id for cpd in model_to_compare.getListOfSpecies()]
    [model_to_compare.removeSpecies(cpd_id ) for cpd_id in all_species if cpd_id not in cpd_id_to_preserve]
    new_id = os.path.basename(os.path.splitext(output)[0])    
    model_to_compare.setId(new_id)
    libsbml.writeSBMLToFile(document_to_compare, output)   

if __name__ == "__main__":
    main()

