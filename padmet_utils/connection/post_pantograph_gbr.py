# -*- coding: utf-8 -*-
"""
Description:
    #TODO

::
    
    usage:
        post_pantograph.py    --ptg_run=DIR --output=FILE [-v]
        post_pantograph.py    --model_metabolic=FILE    --study_metabolic=FILE    --inp=FILE --omcl=FILE --output=FILE [-v]
    
    option:
        -h --help    Show help.
        --model_metabolic=FILE    pathname to the metabolic network of the model (sbml).
        --study_metabolic=FILE    ****.
        --inp=FILE    ****.
        --omcl=FILE    ****.
        --output=FILE    ****.
        -v   print info.
"""
import re
from padmet.utils.sbmlPlugin import parseNotes, parseGeneAssoc
import libsbml
import docopt
import os
import subprocess
import csv

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
        reader = csv.DictReader(f, delimiter = "\t")
        #line[2] or [3] contains genes id with score, delete score and blank
        for line in reader:
            OrtoA = set()
            for gene_id in line["OrtoA"].split(" "):
                if len(gene_id) != 0:
                    try:
                        float(gene_id)
                    except ValueError:
                        OrtoA.add(gene_id)
            OrtoB = set()
            for gene_id in line["OrtoB"].split(" "):
                if len(gene_id) != 0:
                    try:
                        float(gene_id)
                    except ValueError:
                        OrtoB.add(gene_id)
            for gene in OrtoA:
                if gene in list(inp_dict.keys()): print(("%s multiple /!\\" %gene))
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

    #read in model metabolic reactions with or & and in gene assoc'
    reader = libsbml.SBMLReader()
    document = reader.readSBML(model_metabolic)
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    TU_reactions = [rxn for rxn in listOfReactions if "or" in parseNotes(rxn).get("GENE_ASSOCIATION",[""])[0]
    and "and" in parseNotes(rxn).get("GENE_ASSOCIATION",[""])[0]]
    if verbose: print("nb TU reactions: %s" %len(TU_reactions))

    reader_study = libsbml.SBMLReader()
    document_study = reader_study.readSBML(study_metabolic)
    for i in range(document_study.getNumErrors()):
        print(document_study.getError(i).getMessage())
    model_study = document_study.getModel()

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
                new_ga = []
                print("\tTo add, OMCL & INPA valide")
                for subset in match_subsets:
                    correspondance = "("+" and ".join([" or ".join(inp_dict.get(gene,set()).union(inp_dict.get(gene,set()))) for gene in subset])+")"
                    new_ga.append(correspondance)
                new_ga = " or ".join(new_ga)
                print("\t"+new_ga)
                rxn_to_add[rxn.id] = new_ga
        print("%s/%s reactions to add" %(len(rxn_to_add),len(TU_reactions)))

        reader_study = libsbml.SBMLReader()
        document_study = reader_study.readSBML(study_metabolic)
        for i in range(document_study.getNumErrors()):
            print(document_study.getError(i).getMessage())
        model_study = document_study.getModel()

        for rxn_id, new_ga in list(rxn_to_add.items()):
            rxn = model.getReaction(rxn_id)
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
    listOfSpecies = model_study.getListOfSpecies()
    for reaction in listOfReactions:
        if "GENE_ASSOCIATION" not in list(parseNotes(reaction).keys()):
            reactions_to_remove.append(reaction.getId())
    for rId in reactions_to_remove:
        print("Removing %s without gene association" %rId)
        listOfReactions.remove(rId)

    species_in_rxn_temp = [[r.getSpecies() for r in rxn.getListOfReactants()]+[p.getSpecies() for p in rxn.getListOfProducts()] for rxn in listOfReactions]
    species_in_rxn = set()
    [species_in_rxn.update(set(x)) for x in species_in_rxn_temp]
    [listOfSpecies.remove(sId) for sId in set([x.id for x in listOfSpecies]).difference(species_in_rxn)]    

    libsbml.writeSBMLToFile(document_study, output)

def check(value, message):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    if value == None:
        raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
    elif type(value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + libsbml.OperationReturnValue_toString(value).strip() + '"'
            raise TypeError(err_msg)
    else:
        return

if __name__ == "__main__":
    main()

