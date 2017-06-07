# -*- coding: utf-8 -*-
"""
@author: Meziane AITE, meziane.aite@inria.fr
Description:

usage:
    post_pantograph.py    --ptg_run=DIR --output=FILE [-v]
    post_pantograph.py    --model_metabolic=FILE    --study_metabolic=FILE    --inp=FILE --omcl=FILE --output=FILE [-v]

option:
    -h --help    Show help.
    --model_metabolic=FILE    pathname to the metabolic network of the model (sbml).
    --model_faa=FILE    pathname to the proteom of the model (faa)
    --cutoff=FLOAT    cutoff [0:1] for comparing model_metabolic and model_faa. [default: 0.70]. 
    --dict_ids_file=FILE    pathname to the dict associating genes ids from the model_metabolic to the model_faa. line = 
    --output=FILE    output of get_valid_faa (a faa) or get_dict_ids (a dictionnary of gene ids in tsv)
    -v   print info
"""
import re
from padmet.sbmlPlugin import parseNotes, parseGeneAssoc
from padmet.sbmlGenerator import check
from libsbml import *
import docopt
import os
import subprocess

def main():
    args = docopt.docopt(__doc__)
    fold = args["--ptg_run"]
    output = args["--output"]
    verbose = args["-v"]
    if fold:
        if fold.endswith("/"):
            dir_name = os.path.split(fold[:-1])[1]
        else:
            dir_name = os.path.split(fold)[1]
            fold += "/"
        model_metabolic = fold+"metabolic_model.sbml"
        study_metabolic = fold+"original_output_pantograph_"+dir_name+".sbml"
        omcl_rez = fold+"all_orthomcl.out"
        inp_rez = fold+"table.FAA_model.faa-FAA_study.faa"
    else:
        model_metabolic = args["--model_metabolic"]
        study_metabolic = args["--study_metabolic"]
        omcl_rez = args["--omcl"]
        inp_rez = args["--inp"]

    dir_path_gbr = os.path.dirname(os.path.realpath(__file__))+"/grammar-boolean-rapsody.py"
    
    #create dict, k = OrtoA gene_id, v = set of OrtoB orthologus
    inp_dict = {}
    with open(inp_rez, 'r') as f:
        inp_raw_data = [line.split("\t") for line in f.read().splitlines()][1:]
    for line in inp_raw_data:
        #line[2] or [3] contains genes id with score, delete score and blank
        OrtoA = set([i for i in line[2].split(" ") if "." not in i and len(i) != 0])
        OrtoB = set([i for i in line[3].split(" ") if "." not in i and len(i) != 0])
        for gene in OrtoA:
            if gene in inp_dict.keys(): print("%s multiple /!\\" %gene)
            inp_dict[gene] = OrtoB

    #create a dict K = FAA_model gene id, V = set of FAA_study gene ortho
    omcl_dict = {}
    with open(omcl_rez, 'r') as f:
        omcl_raw_data = [line.split("\t")[1] for line in f.read().splitlines()]
    ortho = []
    for line in omcl_raw_data:
        all_genes = [i for i in line.split(" ") if len(i) != 0]
        if any([("FAA_model" in i) for i in all_genes]) and any([("FAA_study" in i) for i in all_genes]):
            ortho.append(all_genes)
    for genes in ortho:
        study_genes = set([g.replace("(FAA_study)","") for g in genes if "FAA_study" in g])
        model_genes = [g.replace("(FAA_model)","") for g in genes if "FAA_model" in g]
        for g in model_genes:
            omcl_dict[g] = study_genes
            
    ortho_in_omcl = set(omcl_dict.keys())
    ortho_in_inp = set(inp_dict.keys())
    ortho_in_omcl_and_inp = ortho_in_omcl.intersection(ortho_in_inp)


    reader = SBMLReader()
    document = reader.readSBML(model_metabolic)
    for i in range(document.getNumErrors()):
        print (document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    TU_reactions = [rxn for rxn in listOfReactions if "or" in parseNotes(rxn).get("GENE_ASSOCIATION","") 
    and "and" in parseNotes(rxn).get("GENE_ASSOCIATION","")]
    if verbose: print("nb TU reactions: %s" %len(TU_reactions))
    
    if TU_reactions:
        count = 0
        rxn_to_add = {}
        for rxn in TU_reactions:
            match_subsets = []
            count += 1
            if verbose: print("reaction %s/%s %s" %(count, len(TU_reactions), rxn.id))
            #rxn = TU_reactions[0]
            ga = parseNotes(rxn)["GENE_ASSOCIATION"][0]
            all_genes = parseGeneAssoc(ga)
            ga_for_gbr = re.sub(r" or " , "|", ga)
            ga_for_gbr = re.sub(r" and " , "&", ga_for_gbr)
            ga_for_gbr = re.sub(r"\s" , "", ga_for_gbr)
            ga_for_gbr = "\"" + ga_for_gbr + "\""
            #for edi test only
            #ga_subsets = eval(subprocess.check_output("python3 grammar-boolean-rapsody.py "+ga_for_gbr, shell=True))
            ga_subsets = eval(subprocess.check_output("python3 "+dir_path_gbr+" "+ga_for_gbr, shell=True))
            
            [match_subsets.append(subset) for subset in ga_subsets if set(subset).issubset(ortho_in_omcl_and_inp)]
            if match_subsets:
                print("\tTo add, OMCL & INPA valide")
                new_ga = " or ".join(["("+" and ".join(subset)+")" for subset in match_subsets])
                print("\t"+new_ga)
                rxn_to_add[rxn] = new_ga
        print("%s/%s reactions to add" %(len(rxn_to_add),len(TU_reactions)))
    
        reader_study = SBMLReader()
        document_study = reader_study.readSBML(study_metabolic)
        for i in range(document_study.getNumErrors()):
            print (document_study.getError(i).getMessage())
        model_study = document_study.getModel()
        
        for rxn, new_ga in rxn_to_add.items():
            notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
            notes += "<p>"+"GENE_ASSOCIATION:" + new_ga + "</p>"
            notes += "</body>"
            check(rxn.setNotes(notes), 'set notes %s' %notes)
            allSpecies = set([p.getSpecies() for p in rxn.getListOfProducts()]).union(set([r.getSpecies() for r in rxn.getListOfReactants()]))
            for species_id in allSpecies:
                if model_study.getSpecies(species_id) is None:
                    print("%s not in study model" %species_id)
                    species_sbml = model.getSpecies(species_id)
                    model_study.addSpecies(species_sbml)
            model_study.addReaction(rxn)

    reactions_to_remove = []
    listOfReactions = model_study.getListOfReactions()
    for reaction in listOfReactions:
        if "GENE_ASSOCIATION" not in parseNotes(reaction).keys():
            reactions_to_remove.append(reaction.getId())
    for rId in reactions_to_remove:
        print("Removing %s without gene association" %rId)
        listOfReactions.remove(rId)

    writeSBMLToFile(document_study, output)


if __name__ == "__main__":
    main()

