#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Contains all necessary functions to generate wikiPages from a padmet file and update 
    a wiki online. Require WikiManager module (with wikiMate,Vendor)

::

    usage:
        wikiGenerator.py --padmetSpec=FILE --output=DIR --model_id=STR --model_name=STR [--padmetRef=FILE] [--log_file=FILE] [-v]
        wikiGenerator.py --aureme_run=DIR --padmetSpec=ID -v
    
    options:
        -h --help     Show help.
        --padmetSpec=FILE    path to padmet file.
        --output=DIR    path to folder to create with all wikipages in subdir.
        --model_id=STR    id of the model to use in mainpage of wiki.
        --model_name=FILE    name of the model to use in mainpage of wiki.
        --padmetRef=FILE    path to padmet of reference, ex: metacyc_xx.padmet, if given, able to calcul pathway rate completion.
        --log_file=FILE    log file from an aureme run, use this file to create a wikipage with all the command used during the aureme run.
        --aureme_run=DIR    can use an aureme run as input, will use from config file information for model_id and log_file and padmetRef.
        -v    print info.
"""
from padmet.classes import PadmetRef
from padmet.classes import PadmetSpec
import os
import shutil
from itertools import chain
from collections import Iterable
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import math
import docopt
import re

def main(): 
    global padmetSpec, padmetRef, wiki_folder, full_sources_dict, all_categories, all_tools, all_sources, total_pwy_id, total_cpd_id, all_rxns, all_genes, all_pwys, all_cpds, verbose, ext_link
    all_categories, all_tools, all_sources, total_pwy_id, total_cpd_id = set(), set(), set(), set(), set()
    full_sources_dict = dict()
    #files to upload: folder genomic_data, all sbml in output ortho, annot, external, seeds, targets
    args = docopt.docopt(__doc__)
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    try:
        db = padmetSpec.info["DB_info"]["DB"].lower()
        if db == "metacyc":
            ext_link = {"Reaction": "http://metacyc.org/META/NEW-IMAGE?object=",
                        "Pathway": "http://metacyc.org/META/NEW-IMAGE?object=",
                        "Metabolite": "http://metacyc.org/META/NEW-IMAGE?object="}
        elif db == "bigg":
            ext_link = {"Reaction":"http://bigg.ucsd.edu/universal/reactions/",
                        "Metabolite":"http://bigg.ucsd.edu/universal/metabolites/",
                        "Pathway":"http://www.genome.jp/dbget-bin/www_bget?"}
        else:
            raise KeyError
    except KeyError:
        ext_link = {}
    if args["--padmetRef"]:#TODO: finir case aureme_run vs all args given in param
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    verbose = args["-v"]
    model_id, model_name = args["--model_id"], args["--model_name"]
    wiki_folder = args["--output"]
    if not wiki_folder.endswith("/"): wiki_folder += "/"  
    log_file = args["--log_file"]
    createDirectory()
    all_rxns = [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "reaction"]
    all_genes = [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "gene"]
    for rxn_node in all_rxns:
        create_biological_page("Reaction", rxn_node, wiki_folder+"reactions/")
    for gene_node in all_genes:    
        create_biological_page("Gene", gene_node, wiki_folder+"genes/")
    all_pwys = [node for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_pwy_id]
    all_cpds = [node for (node_id, node) in padmetSpec.dicOfNode.items() if node_id in total_cpd_id]
    for pwy_node in all_pwys:
        create_biological_page("Pathway", pwy_node, wiki_folder+"pathways/")
    for cpd_node in all_cpds:
        create_biological_page("Metabolite", cpd_node, wiki_folder+"metabolites/")

    create_navigation_page(wiki_folder+"/navigation/")
    #create_venn()
    create_main(model_id, model_name)
    if log_file:
        create_log_page(log_file, wiki_folder+"/navigation/")
    

def createDirectory():
    """
    create the folders genes, reactions, metabolites, pathways in the folder dirPath/
    if already exist, it will replace old folders (and delete old files)
    """
    global wiki_folder
    #simple check that dirPath is a dir:
    dirNames = ["genes","reactions","metabolites","pathways","navigation","files"]
    #creatings the directory which will contains the wiki pages
    for d in dirNames:
        if not os.path.exists(wiki_folder+d):
            if verbose: print("Creating directory: "+wiki_folder+d)
            os.makedirs(wiki_folder+d)
        else:
            if verbose: print("The directory "+wiki_folder+d+" already exist. Old pages will be deleted")
            shutil.rmtree(wiki_folder+d)
            os.makedirs(wiki_folder+d)
    
def create_venn():
    """
    
    """
    if verbose: print("Venn Diagramm")
    #'ARGSUCCINSYN-RXN'
    categories_dict ={}
    """
    all_categories = ["orthology","annotation","gap-filling","manual","microbiont"]
    for category in all_categories:
        categories_dict[category] = set()
    """
    for rxn_id, rxn_src_dict in list(full_sources_dict.items()):
        for category in list(rxn_src_dict.keys()):
            try:
                categories_dict[category].add(rxn_id)
            except KeyError:
                categories_dict[category] = set(rxn_id)
                
    
    labels = get_labels(list(categories_dict.values()))
    fig, ax = venn4(labels, names=list(categories_dict.keys()))
    fig.savefig(wiki_folder+"files/venn.png")

def copy_io_files():
    """
    """
    #toDo in futur
        
def create_main(model_id, model_name):
    if verbose: print("Main page")
    ### create main page
    for line in main_template:
        main_template[main_template.index(line)] = line.replace("MODEL_ID",model_id).replace("MODEL_NAME",model_name)
    final_network_index = main_template.index([line for line in main_template if line.startswith("The automatic")][0])
    main_template[final_network_index] = main_template[final_network_index].replace("NB_RXN", str(len(all_rxns))).replace("NB_CPD", str(len(all_cpds))).replace("NB_PWY", str(len(all_pwys))).replace("NB_GENE", str(len(all_genes)))
    reconstruct_summary = {"ANNOTATION":0,"ORTHOLOGY":{},"MANUAL":0,"GAP-FILLING":0,"MICROBIONT":0,"HOST":0}
    for rec_node in [node for node in list(padmetSpec.dicOfNode.values()) if node.type == "reconstructionData"]:
        cat = rec_node.misc["CATEGORY"][0]
        if cat == "ORTHOLOGY":
            source = rec_node.misc["SOURCE"][0].replace("OUTPUT_PANTOGRAPH_","")
            try:
                reconstruct_summary["ORTHOLOGY"][source] += 1
            except KeyError:
                reconstruct_summary["ORTHOLOGY"][source] = 1
        else:
            reconstruct_summary[cat] += 1
    index = 1
    if reconstruct_summary["ANNOTATION"] != 0:
        main_template.insert(final_network_index+index, "* Based on annotation data:")
        index += 1
        main_template.insert(final_network_index+index, "** Tool: [http://bioinformatics.ai.sri.com/ptools/ Pathway tools]")
        index += 1
        main_template.insert(final_network_index+index, "*** Creation of a metabolic network containing "+str(reconstruct_summary["ANNOTATION"])+" reactions")
        index += 1
    if reconstruct_summary["ORTHOLOGY"]:
        main_template.insert(final_network_index+index, "* Based on orthology data:")
        index += 1
        main_template.insert(final_network_index+index, "** Tool: [http://pathtastic.gforge.inria.fr Pantograph]")
        index += 1        
        for k,v in list(reconstruct_summary["ORTHOLOGY"].items()):
            main_template.insert(final_network_index+index, "*** From template ''"+k+"'' creation of a metabolic network containing: "+str(v)+" reactions")
            index += 1
    if reconstruct_summary["MANUAL"] != 0:
        main_template.insert(final_network_index+index, "* Based on expertise:")
        index += 1
        main_template.insert(final_network_index+index, "*** "+str(reconstruct_summary["MANUAL"])+" reaction(s) added")
        index += 1
    if reconstruct_summary["GAP-FILLING"] != 0:
        main_template.insert(final_network_index+index, "* Based on gap-filling:")
        index += 1
        main_template.insert(final_network_index+index, "** Tool: [https://pypi.python.org/pypi/meneco meneco]")
        index += 1
        main_template.insert(final_network_index+index, "*** "+str(reconstruct_summary["GAP-FILLING"])+" reaction(s) added")

    fileName = wiki_folder+"/navigation/Main_Page"
    with open(fileName,'w') as f:
        for line in main_template:
            f.write(line+"\n")                

def create_navigation_page(output_folder):
    """
    
    """
    sideBarData = ["* navigation","** mainpage|mainpage-description","** workflow|workflow command history","** randompage-url|randompage","** Special:ListFiles|Files","* Metabolic network components"]

    category = "Reaction"
    sideBarData.append("** Category:"+category+"|"+category)
    fileName = output_folder+"Category:Reaction"
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:Reaction]]","| ?common name","| ?ec number",
                      "| ?reconstruction category","| ?reconstruction tool","| ?reconstruction source","| ?gene associated","| ?in pathway","}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "Gene"
    sideBarData.append("** Category:"+category+"|"+category)
    fileName = output_folder+"Category:Gene"
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:Gene]]", "| ?reaction associated", "| ?pathway associated","}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "Pathway"
    sideBarData.append("** Category:"+category+"|"+category)
    fileName = output_folder+"Category:Pathway"
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:Pathway]]","| ?common name","| ?reaction found","| ?total reaction","| ?completion rate","}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "Metabolite"
    sideBarData.append("** Category:"+category+"|"+category)
    fileName = output_folder+"Category:Metabolite"
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:Metabolite]]","| ?common name","| ?consumed by","| ?produced by","| ?reversible reaction associated","}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")




    sideBarData.append("* Reconstruction categories")
    [sideBarData.append("** "+rec_category+"|"+rec_category) for rec_category in sorted(all_categories)]
    for rec_category in all_categories:
        fileName = output_folder+rec_category
        if verbose: print("Reconstruction category: %s" %rec_category)
        dataInArray = ["{{#ask: [[Category:Reaction]] [[reconstruction category::"+rec_category+"]]","| ?common name","| ?ec number",
                          "| ?reconstruction category","| ?reconstruction tool","| ?reconstruction source","| ?reconstruction comment","| ?gene associated","| ?in pathway","}}"]
        with open(fileName,'w') as f:
            for line in dataInArray:
                f.write(line+"\n")

    sideBarData.append("* Reconstruction tools")
    [sideBarData.append("** "+rec_tool+"|"+rec_tool) for rec_tool in sorted(all_tools)]
    for rec_tool in all_tools:
        fileName = output_folder+rec_tool
        if verbose: print("Reconstruction tool: %s" %rec_tool)
        dataInArray = ["{{#ask: [[Category:Reaction]] [[reconstruction tool::"+rec_tool+"]]","| ?COMMON NAME","| ?ec number",
                          "| ?reconstruction category","| ?reconstruction tool","| ?reconstruction source","| ?reconstruction comment","| ?gene associated","| ?in pathway","}}"]
        with open(fileName,'w') as f:
            for line in dataInArray:
                f.write(line+"\n")

    sideBarData.append("* Reconstruction sources")
    [sideBarData.append("** "+rec_source+"|"+rec_source) for rec_source in sorted(all_sources)]
    for rec_source in all_sources:
        fileName = output_folder+rec_source
        if verbose: print("Reconstruction source: %s" %rec_source)
        dataInArray = ["{{#ask: [[Category:Reaction]] [[reconstruction source::"+rec_source+"]]","| ?COMMON NAME","| ?ec number",
                          "| ?reconstruction category","| ?reconstruction tool","| ?reconstruction source","| ?reconstruction comment","| ?gene associated","| ?in pathway","}}"]
        with open(fileName,'w') as f:
            for line in dataInArray:
                f.write(line+"\n")

    if verbose: print("SideBar page")
    with open(output_folder+"MediaWiki:Sidebar", 'w') as f:
        for line in sideBarData:
            f.write(line+"\n")

def create_biological_page(category, page_node, output_folder):
    """
    
    """
    global padmetSpec, all_pwy_id, all_cpd_id

    if "/" in page_node.id:
        print("%s contains not allowed caractere" %page_node.id)
        fileName = output_folder + page_node.id.replace("/","__47__")
    else:
        fileName = output_folder + page_node.id
    if verbose: print("%s: %s" %(category, page_node.id))
    #stock in properties: all properties associated to the current page.
    #properties = [{{#set PROPERTY_X:VALUE_1|...|VALUE_N}}, ...]
    properties = []
    #ext_link is used to create external link to the database of reference if known
    if ext_link.get(category):
        dataInArray = ['[[Category:'+category+']]',
        '== '+category+' ['+ext_link.get(category)+page_node.id+' '+page_node.id+'] ==']
    else:
        dataInArray = ['[[Category:'+category+']]',
        '== '+category+' '+page_node.id+' ==']
    #extracting all data in dict misc of node. ex: misc= {"A":["X","Y"]}
    #adding in page: ** A: *X *Y
    for k,v in page_node.misc.items():
        k = k.replace("-"," ").lower()
        line = '* '+k+':'
        dataInArray.append(line)
        for i in v:
            if k == "ec number":
                line = "** [http://enzyme.expasy.org/EC/"+i.replace("EC-","")+" "+i+"]"
            elif k == "taxonomic range" and ext_link.get(category):
                line = "** ["+ext_link.get(category)+i+" "+i+"]"
            else:
                line = '** '+i
            #add_property is used to stock all added properties and add them at the end of the page
            add_property(properties, k, [i])
            dataInArray.append(line)
    #if node is linked to a name node, extracts synonyms from name node misc dict, k = 'LABEL'
    #keyError: no relations in for current node [node => X]
    #names is empty: bool(names) == False
    try:
        names = [padmetSpec.dicOfNode[rlt.id_out].misc["LABEL"] for rlt in padmetSpec.dicOfRelationIn[page_node.id] if rlt.type == "has_name"]
        if names:
            names = names[0]
        else:
            raise KeyError
    except KeyError:
        names = None
    line = "* Synonym(s):"
    dataInArray.append(line)
    if names:
        add_property(properties, "common name", names)
        for name in names:
            line = "** "+name
            dataInArray.append(line)

    #For each category, extract in a specific way the information
    dataInArray.append("")
    if category == "Reaction":
        #extract direction information
        direction = page_node.misc["DIRECTION"][0]
        if direction == "UNKNOWN":
            direction = " '''=>/<=>''' "
        elif direction == "REVERSIBLE":
            direction = " '''<=>''' "
        elif direction == "LEFT-TO-RIGHT":
            direction = " '''=>''' "
        
        #global var total_cpd_id will contains all ids of compounds involved in a reaction
        total_cpd_id.update([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[page_node.id] if rlt.type in ["consumes","produces"]])
        # Recovering the formula
        dataInArray.append('== Reaction Formula ==', )
        #reactants and products: list of str: ['stoichio [[cpd_id]][cpd_compart]',...]
        reactants = [rlt.misc["STOICHIOMETRY"][0]+" [["+rlt.id_out+"]]["+rlt.misc["COMPARTMENT"][0]+"]"
        for rlt in padmetSpec.dicOfRelationIn.get(page_node.id, None) if rlt.type == "consumes"]
        products = [rlt.misc["STOICHIOMETRY"][0]+" [["+rlt.id_out+"]]["+rlt.misc["COMPARTMENT"][0]+"]"
        for rlt in padmetSpec.dicOfRelationIn.get(page_node.id, None) if rlt.type == "produces"]
        #join each list with '+' and direction in center
        formula_id = " '''+''' ".join(reactants)+direction+" '''+''' ".join(products)
        dataInArray.append("* With identifiers:")
        dataInArray.append("** "+formula_id)            
        dataInArray.append("* With common name(s):")
        #if all involved compounds have a common name, creating the same formula but with common name
        try:
            reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetSpec.dicOfNode[rlt.id_out].misc["COMMON-NAME"][0]+"["+rlt.misc["COMPARTMENT"][0]+"]" 
            for rlt in padmetSpec.dicOfRelationIn.get(page_node.id, None) if rlt.type == "consumes"]
            products = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetSpec.dicOfNode[rlt.id_out].misc["COMMON-NAME"][0]+"["+rlt.misc["COMPARTMENT"][0]+"]"
            for rlt in padmetSpec.dicOfRelationIn.get(page_node.id, None) if rlt.type == "produces"]
            formula_cname = " '''+''' ".join(reactants)+direction+" '''+''' ".join(products)+"\n"
            dataInArray.append("** "+formula_cname)
        except KeyError:
            dataInArray.append("**")

        dataInArray.append('== Genes associated with this reaction  ==')
        #get all relations, type == "is_linked_to"
        linked_rlt = [rlt for rlt in padmetSpec.dicOfRelationIn[page_node.id] if rlt.type == "is_linked_to"]
        if linked_rlt:
            dataInArray.append('Genes have been associated with this reaction based on different elements listed below.')
            add_property(properties, "gene associated", [rlt.id_out for rlt in linked_rlt])
            for rlt in linked_rlt:
                gene_id = rlt.id_out
                dataInArray.append("* Gene: [[%s]]" %gene_id)
                #a is_linked_to rlt have in misc a key "SOURCE:ASSIGNMENT"
                #the value can be only the source, ex: OUTPUT_PANTOGRAPH_X
                #or the source and the known assignment: SILI_ANNOTATION:EC-NUMBER
                sources = rlt.misc["SOURCE:ASSIGNMENT"]
                for src_data in sources:
                    #ValueError: can't split and get 2 value == no known assignment
                    try:
                        src, assignment = src_data.split(":")
                    except ValueError:
                        src = src_data
                        assignment = None
                    category = next(node.misc["CATEGORY"][0] for node in list(padmetSpec.dicOfNode.values()) if (node.type == "reconstructionData" and node.misc.get("SOURCE",[None])[0] == src))
                    src = src.lower()
                    category = category.lower()
                    if src.startswith("output_pantograph_"): src = src.replace("output_pantograph_","")
                    src = "-".join([category, src])
                    dataInArray.append("** Source: [[%s]]" %src)
                    if assignment:
                        dataInArray.append("*** Assignment: %s" %assignment)
        
        dataInArray.append('== Pathways  ==')
        #set of pathways id associated to the reaction
        pathways_ids = set([rlt.id_out for rlt in padmetSpec.dicOfRelationIn[page_node.id]
        if rlt.type == "is_in_pathway"])
        #update global var total_pwy_id containing all pathways involved in the metabolic network
        total_pwy_id.update(pathways_ids)
        for pwy_id in pathways_ids:
            #recovering the nb of reactions associated to the pathway
            if padmetRef is not None:
                try:
                    nbReactionsTotal = len([rlt for rlt in padmetRef.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway"])
                # If keyError: pathway not in padmetRef, pathway added manualy
                except KeyError: 
                    nbReactionsTotal = "NA"
            else:
                nbReactionsTotal = "NA"
            nbReactionsFound = len([rlt for rlt in padmetSpec.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway"])

            #extract pwy common name
            pwy_cname = padmetSpec.dicOfNode[pwy_id].misc.get("COMMON-NAME",[None])[0]
            if pwy_cname:
                line = "* [["+pwy_id+"]], "+pwy_cname+":"
            else:
                line = "* [["+pwy_id+"]]:"

            #if external link known, adding external link to the pathway
            #line: PWY_ID, common name, extern link to PWY
            if ext_link.get("Pathway"):
                line += ' ['+ext_link.get("Pathway")+pwy_id+' '+pwy_id+']'
            dataInArray.append(line)
            dataInArray.append("** '''"+str(nbReactionsFound)+"''' reactions found over '''"+str(nbReactionsTotal)+"''' reactions in the full pathway")
        add_property(properties, "in pathway", pathways_ids)

        dataInArray.append('== Reconstruction information  ==')
        reconstruction_data = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[page_node.id] if rlt.type == "has_reconstructionData"]
        #src_data = {category:{source:{comment:comment, tool:tool}}}
        src_data = {}
        for category in set([rec_data_node.misc["CATEGORY"][0].lower() for rec_data_node in reconstruction_data]):
            src_data[category] = dict()
        rxn_srcs, rxn_tools, rxn_comments, rxn_categories = set(), set(), set(), set()
        for rec_data_node in reconstruction_data:
            #if found, lower to standardize
            category = rec_data_node.misc.get("CATEGORY",[None])[0]
            if category: 
                category = category.lower()
                rxn_categories.add(category)
            tool = rec_data_node.misc.get("TOOL",[None])[0]
            if tool: 
                tool = tool.lower()
                rxn_tools.add(tool)
            comment = rec_data_node.misc.get("COMMENT",[None])[0]
            if comment: 
                comment = comment.lower()
                rxn_comments.add(comment)

            source = rec_data_node.misc.get("SOURCE",[None])[0]
            source = source.lower()
            if source.startswith("output_pantograph_"):
                source = source.replace("output_pantograph_","")
            source = category+"-"+source
            rxn_srcs.add(source)
            src_data[category][source] = {"comment":comment,"tool":tool}
        all_categories.update(rxn_categories)
        add_property(properties, "reconstruction category", rxn_categories)
        if rxn_srcs:
            all_sources.update(rxn_srcs)
            add_property(properties, "reconstruction source", rxn_srcs)
        if rxn_tools:
            all_tools.update(rxn_tools)
            add_property(properties, "reconstruction tool", rxn_tools)
        if rxn_comments:
            add_property(properties, "reconstruction comment", rxn_comments)

        #udpating global var full_sources_dict: k = rxn_id, v = previous src_data dict
        full_sources_dict[page_node.id] = src_data
        for category, category_data in list(src_data.items()):
            dataInArray.append("* Category: [["+category+"]]")
            for source, source_data in list(category_data.items()):
                dataInArray.append("** Source: [["+source+"]]")
                if source_data["tool"]:
                    dataInArray.append("*** Tool: [["+source_data["tool"]+"]]")
                    if source_data["comment"]:
                        dataInArray.append("**** Comment: [["+source_data["comment"]+"]]")
                elif source_data["comment"]:
                    dataInArray.append("*** Comment: [["+source_data["comment"]+"]]")

    elif category == "Gene":
        dataInArray.append("== Reactions associated ==")
        linked_rlt = [rlt for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type == "is_linked_to"]
        if linked_rlt:
            pwy_assoc = set()
            add_property(properties, "reaction associated", [rlt.id_in for rlt in linked_rlt])
            for rlt in linked_rlt:
                rxn_id = rlt.id_in
                [pwy_assoc.add(rlt_pwy.id_out) for rlt_pwy in padmetSpec.dicOfRelationIn[rxn_id] if rlt_pwy.type == "is_in_pathway"]
                dataInArray.append("* Reaction: [[%s]]" %rxn_id)
                sources = rlt.misc["SOURCE:ASSIGNMENT"]
                for src_data in sources:
                    try:
                        src, assignment = src_data.split(":")
                    except ValueError:
                        src = src_data
                        assignment = None
                    category = next(node.misc["CATEGORY"][0] for node in list(padmetSpec.dicOfNode.values()) if (node.type == "reconstructionData" and node.misc.get("SOURCE",[None])[0] == src))
                    src = src.lower()
                    category = category.lower()
                    if src.startswith("output_pantograph_"): src = src.replace("output_pantograph_","")
                    src = "-".join([category, src])
                    dataInArray.append("** Source: [[%s]]" %src)
                    if assignment:
                        assignment = assignment.lower()
                        dataInArray.append("*** Assignment: "+assignment)
            dataInArray.append("== Pathways associated ==")
            if pwy_assoc :
                add_property(properties, "pathway associated", pwy_assoc)
                for pwy_id in pwy_assoc:
                    dataInArray.append("* [["+pwy_id+"]]")

    elif category == "Pathway":
        dataInArray.append("== Reaction(s) found ==")
        #recovering the nb of reactions associated to the pathway
        if padmetRef is not None:
            try:
                reactionsTotal = [rlt.id_in for rlt in padmetRef.dicOfRelationOut[page_node.id] if rlt.type == "is_in_pathway"]
            # If keyError: pathway not in padmetRef, pathway added manualy
            except KeyError: 
                reactionsTotal = [rlt.id_in for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type == "is_in_pathway"]
        else:
                reactionsTotal = [rlt.id_in for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type == "is_in_pathway"]
        reactionsFound = [rlt.id_in for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type == "is_in_pathway"]
        reactionsMissing = [rxn_id for rxn_id in reactionsTotal if rxn_id not in reactionsFound]
        pwy_ratio = round(float(len(reactionsFound))/float(len(reactionsTotal)),2)*100
        dataInArray.append(" '''%s''' reactions found over '''%s''' reactions in the full pathway" %(len(reactionsFound), len(reactionsTotal)))

        add_property(properties, "reaction found", [len(reactionsFound)])
        add_property(properties, "total reaction", [len(reactionsTotal)])
        add_property(properties, "completion rate", [pwy_ratio])
        for rxn_id in reactionsFound:
            gene_assoc = [rlt.id_out for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"]
            dataInArray.append("* [["+rxn_id+"]]")
            if gene_assoc:
                dataInArray.append("** %s associated gene(s):" %(len(gene_assoc)))
                for gene_id in gene_assoc:
                    dataInArray.append("*** [["+gene_id+"]]")
            else:
                dataInArray.append("** 0 associated gene:")

            try:
                src_data = full_sources_dict[rxn_id]
                sources = set()
                [sources.update(list(category_data.keys())) for category_data in list(src_data.values())]
                if sources:
                    dataInArray.append("** %s reconstruction source(s) associated:" %(len(sources)))
                    for src in sources:
                        dataInArray.append("*** [["+src+"]]")
            except KeyError:
                #if keyError, not a reaction but a pathway in a pathway
                pass

        dataInArray.append("== Reaction(s) not found ==")
        if reactionsMissing:
            if ext_link.get("Reaction"):
                for rxn_id in reactionsMissing:
                    dataInArray.append("* ["+ext_link.get("Reaction")+rxn_id+" "+rxn_id+"]")
            else:
                for rxn_id in reactionsMissing:
                    dataInArray.append("* "+rxn_id)

    elif category == "Metabolite":
        rxn_cp = {"c":set(), "p": set(), "cp": set()}
        for rlt in [rlt for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type in ["consumes", "produces"]]:
            rxn_dir = padmetSpec.dicOfNode[rlt.id_in].misc["DIRECTION"][0]
            if rxn_dir == "REVERSIBLE":
                rxn_cp["cp"].add(rlt.id_in)
            else:
                if rlt.type == "consumes":
                    rxn_cp["c"].add(rlt.id_in)
                elif rlt.type == "produces":
                    rxn_cp["p"].add(rlt.id_in)
        
        dataInArray.append("== Reaction(s) known to consume the compound ==")
        if rxn_cp["c"]:
            add_property(properties, "consumed by",rxn_cp["c"])
            for rxn_id in rxn_cp["c"]:
                dataInArray.append("* [["+rxn_id+"]]")
        dataInArray.append("== Reaction(s) known to produce the compound ==")
        if rxn_cp["p"]:
            add_property(properties, "produced by",rxn_cp["p"])
            for rxn_id in rxn_cp["p"]:
                dataInArray.append("* [["+rxn_id+"]]")
        dataInArray.append("== Reaction(s) of unknown directionality ==")
        if rxn_cp["cp"]:
            add_property(properties, "reversible reaction associated",rxn_cp["cp"])
            for rxn_id in rxn_cp["cp"]:
                dataInArray.append("* [["+rxn_id+"]]")

            
    dataInArray.append('== External links  ==')
    try:
        xref_node = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[page_node.id] if rlt.type == "has_xref"][0]
        for db, ids in list(xref_node.misc.items()):
            xrefLink(dataInArray, db, ids)
    except (IndexError, KeyError) as e:
        pass
    
    dataInArray.extend(properties)
    with open(fileName,'w', encoding="utf8") as f:
        for line in dataInArray:
            f.write(line+"\n")
    """                
    for i in dataInArray:
        print(i)
    print("\n")
    """
def add_property(properties, prop_id, prop_values):
    prop_values = [str(i) for i in prop_values]
    start_line = "{{#set: "+prop_id+"="
    values_part = "|".join(prop_values)
    end_line = "}}"
    toInsert = start_line + values_part + end_line
    properties.append(toInsert)

def xrefLink(dataInArray, db, ids):
    if db == "METACYC":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid="+_id+" "+_id+"]"
            dataInArray.append(toInsert)
    elif db == "UNIPROT":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.uniprot.org/uniprot/"+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "KEGG" or db.startswith("LIGAND"):
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.genome.jp/dbget-bin/www_bget?"+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "RHEA":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.ebi.ac.uk/rhea/reaction.xhtml?id="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "WIKIPEDIA":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://en.wikipedia.org/wiki/"+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "CHEBI":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.ebi.ac.uk/chebi/searchId.do?chebiId="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "PUBCHEM":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "ECOCYC":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://metacyc.org/ECOLI/NEW-IMAGE?object="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "CHEMSPIDER":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://www.chemspider.com/Chemical-Structure."+_id+".html "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "umbbd-compounds":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://umbbd.ethz.ch/servlets/pageservlet?ptype=c&compID="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "ARACYC":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://metacyc.org/ARA/NEW-IMAGE?object="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "PIR":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://pir.georgetown.edu/cgi-bin/nbrfget?uid="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "NCI":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://cactus.nci.nih.gov/ncidb2.2/?nsc="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    elif db == "knapsacK":
        dataInArray.append("* "+db+":")
        for _id in ids:
            toInsert = "** [http://kanaya.naist.jp/knapsack_jsp/information.jsp?word="+_id+" "+_id+"]"
            dataInArray.append(toInsert)

    else:
        for _id in ids:
            toInsert = "* "+db+" : "+_id
            dataInArray.append(toInsert)

def get_labels(data, fill=["number"]):
    """    
    get a dict of labels for groups in data
    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill=["number"])
    Out[12]:
    {'001': '0',
    '010': '5',
    '011': '0',
    '100': '3',
    '101': '2',
    '110': '2',
    '111': '3'}
    
    Parameters
    ----------    
    data: list
        data to get label for
    fill: 
        ["number"|"logic"|"percent"]

    Returns
    -------
    dict:
        a dict of labels for different sets
    """
    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                             # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    labels = {k: "" for k in set_collections}
    if "logic" in fill:
        for k in set_collections:
            labels[k] = k + ": "
    if "number" in fill:
        for k in set_collections:
            labels[k] += str(len(set_collections[k]))
    if "percent" in fill:
        data_size = len(s_all)
        for k in set_collections:
            labels[k] += "(%.1f%%)" % (100.0 * len(set_collections[k]) / data_size)

    return labels

def venn4(labels, names=['A', 'B', 'C', 'D'], **options):
    """
    plots a 4-set Venn diagram

    Parameters
    ----------    
    labels: dict
        a label dict where keys are identified via binary codes ('0001', '0010', '0100', ...),
        hence a valid set could look like: {'0001': 'text 1', '0010': 'text 2', '0100': 'text 3', ...}.
        unmentioned codes are considered as ''.
    names: list
        group names

    Returns
    -------
    set
        (Figure, AxesSubplot), pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(4)])
    figsize = options.get('figsize', (12, 12))
    dpi = options.get('dpi', 96)
    
    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)
    
    # body   
    draw_ellipse(fig, ax, 0.350, 0.400, 0.72, 0.45, 140.0, colors[0])
    draw_ellipse(fig, ax, 0.450, 0.500, 0.72, 0.45, 140.0, colors[1])
    draw_ellipse(fig, ax, 0.544, 0.500, 0.72, 0.45, 40.0, colors[2])
    draw_ellipse(fig, ax, 0.644, 0.400, 0.72, 0.45, 40.0, colors[3])
    draw_text(fig, ax, 0.85, 0.42, labels.get('0001', ''))
    draw_text(fig, ax, 0.68, 0.72, labels.get('0010', ''))
    draw_text(fig, ax, 0.77, 0.59, labels.get('0011', ''))
    draw_text(fig, ax, 0.32, 0.72, labels.get('0100', ''))
    draw_text(fig, ax, 0.71, 0.30, labels.get('0101', ''))
    draw_text(fig, ax, 0.50, 0.66, labels.get('0110', ''))
    draw_text(fig, ax, 0.65, 0.50, labels.get('0111', ''))
    draw_text(fig, ax, 0.14, 0.42, labels.get('1000', ''))
    draw_text(fig, ax, 0.50, 0.17, labels.get('1001', ''))
    draw_text(fig, ax, 0.29, 0.30, labels.get('1010', ''))
    draw_text(fig, ax, 0.39, 0.24, labels.get('1011', ''))
    draw_text(fig, ax, 0.23, 0.59, labels.get('1100', ''))
    draw_text(fig, ax, 0.61, 0.24, labels.get('1101', ''))
    draw_text(fig, ax, 0.35, 0.50, labels.get('1110', ''))
    draw_text(fig, ax, 0.50, 0.38, labels.get('1111', ''))
    
    # legend
    draw_text(fig, ax, 0.13, 0.18, names[0], colors[0])
    draw_text(fig, ax, 0.18, 0.83, names[1], colors[1])
    draw_text(fig, ax, 0.82, 0.83, names[2], colors[2])
    draw_text(fig, ax, 0.87, 0.18, names[3], colors[3])
    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
    return fig, ax

def draw_ellipse(fig, ax, x, y, w, h, a, fillcolor):
    e = patches.Ellipse(
        xy=(x, y),
        width=w,
        height=h,
        angle=a,
        color=fillcolor)
    ax.add_patch(e)

def draw_text(fig, ax, x, y, text, color=[0, 0, 0, 1]):
    ax.text(
        x, y, text,
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=14,
        color=color)

default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [246, 236, 86, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]


main_template = ["== MODEL_IDGEM description ==",
     "== Automatic reconstruction with [http://aureme.genouest.org AuReMe] ==",
     "Model summary: [[MEDIA:summary.txt|summary]]",
     "",
     "Download '''AuReMe''' Input/Output [LINK OR MEDIA data]",
     "",
     "The automatic reconstruction of ''MODEL_NAME'' results to a Genome scale [[MEDIA:model.xml|Model]] containing NB_RXN reactions, NB_CPD metabolites, NB_GENE genes and NB_PWY pathways. This GeM was obtained based on the following sources:",
     "",
     "[[FILE:venn.png|frameless|border]]",
     "",
     "== Collaborative curation == ",
     "* Suggest reactions to add or remove:",
     "** Download this [[MEDIA:Add_delete_reaction.csv|form]]",
     "* Suggest new reactions to create and add:",
     "** Download this [[MEDIA:Reaction_creator.csv|form]]",
     "* '''Follow the examples given in the form(s) to correctly share your suggestions'''",
     "* Send the filled form(s) to: CONTACT_MAIL"]

def create_log_page(log_file, output_folder):
    """
    
    """
    cmd_regex = '--cmd=\"(.*)\"'
    fileName = output_folder+"workflow"
    log_page = ["=Workflow command history=","","==Command sequence=="]
    with open(log_file, 'r') as f:
        log_data = [line for line in f.read().splitlines() if not line.startswith("#")]
    for cmd_line in log_data:
        re_result = re.search(cmd_regex, cmd_line)
        if re_result:
            full_cmd = re_result.groups(1)[0]
            cmd = full_cmd.split(" ")[0]
            cmd_label,desc = get_cmd_label(cmd)
            if cmd_label and desc:
                if cmd == "curation":
                    curation_file = re.search('DATA=(.*\.\w+)\s?',full_cmd).groups(1)[0]
                    desc = desc.replace('FORM_FILE_NAME', curation_file)
                log_page.append("* '''%s''':" %cmd_label)
                log_page.append("''%s''" %desc)
    log_page.extend(["==Downloads==","You can download the [[MEDIA:log.txt|command log file here]]"])
    #todo downalosalsjka
    with open(fileName, 'w') as f:
        for line in log_page:
            f.write(line+"\n")
        
                 

def get_cmd_label(cmd):
    """
    """
    cmd_label_dict = {'getdb':{'CMD_LABEL':'Get database','DESC':"Display the available reference databases"},
                      'check_input':{'CMD_LABEL':'Check input','DESC':"Check the validity, consistency and presence of input files"},
                      'check_studied_organism_input':{'CMD_LABEL':'Check studied organism input','DESC':"Check if FAA or GBK was given for the studied organism"},
                      'check_model_organism_input':{'CMD_LABEL':'Check model organism input','DESC':"Check (if existing) each folder in orthology based reconstruction"},
                      'sbml_validity':{'CMD_LABEL':"Check SBML validity", 'DESC':"Check the SBML validity of $(METABOLIC_MODEL)."},
                      'faa_validity':{'CMD_LABEL':"Check Fasta validity", 'DESC':"Check if genes IDs in the model metabolic network are the same than in the Fasta file.<BR>If the rate of validated genes IDs is lower than the CUTOFF (see config.txt), raise an error and break the workflow."},
                      'check_gap_filling_input':{'CMD_LABEL':"Check inputs for gap-filling", 'DESC':"Check the seeds, targets and artefacts files."},
                      'curation':{'CMD_LABEL':"Manual curation",'DESC':"Apply the curation described in the form file FORM_FILE_NAME."},
                      'annotation_based':{'CMD_LABEL':"Annotation based reconstruction", 'DESC':"Extract network data from Pathway Tools annotation output."},
                      'pathwaytools':{'CMD_LABEL':"Pathway Tools", 'DESC':"For each folder in /annotation_based_reconstruction, check or generate from pgdb files (.dat) the SBML file."},
                      'orthology_based':{'CMD_LABEL':"Orthology based reconstruction", 'DESC':"Run the orthology based reconstruction."},
                      'inparanoid':{'CMD_LABEL':"Inparanoid", 'DESC':"Run Inparanoid : search for orthologs."},
                      'OMCL':{'CMD_LABEL':"OrthoMCL", 'DESC':"Run OrthoMCL : search for orthologs."},
                      'mp_pantograph' :{'CMD_LABEL':"Pantograph multiprocessed", 'DESC':"Run Pantograph : merges OrthoMCL and inparanoid results.<BR>Quickened by mulitprocessing."},
                      'pantograph':{'CMD_LABEL':"Pantograph", 'DESC':"Run Pantograph : merges OrthoMCL and inparanoid results."},
                      'draft':{'CMD_LABEL':"Create draft network", 'DESC':"Merges all available networks from the /networks directory into one metabolic network.<BR>Merge all data on the studied species."},
                      'gap_filling':{'CMD_LABEL':"Run gap-filling", 'DESC':"Calculate the gap-filling solution and generate the metabolic network, completed with the gap-filling solution."},
                      'gap_filling_solution':{'CMD_LABEL':"Calculate gap-filling solution", 'DESC':"Only calculate the gap-filling solution."},
                      'meneco':{'CMD_LABEL':"Meneco", 'DESC':"Run Meneco : a gap-filling reconstruction method."},
                      'final':{'CMD_LABEL':"Final network", 'DESC':"Generate the final metabolic network, once applyed all the reconstruction methods."},
                      'gbk_to_faa':{'CMD_LABEL':"GBK to Fasta", 'DESC':"Export a GeneBank (.gbk) file in Fasta (faa) format."},
                      'pgdb_to_padmet':{'CMD_LABEL':"PGDB to PADMet", 'DESC':"Export a PGDB (.dat Pathway Tools files) in PADMet (.padmet) format."},
                      'padmet_to_sbml':{'CMD_LABEL':"PADMet to SBML", 'DESC':"Export a PADMet (.padmet) file in the SBML format."},
                      'compounds_to_sbml':{'CMD_LABEL':"Compounds to SBML", 'DESC':"Export a list of compounds (.txt) in the SBML format."},
                      'sbml_mapping':{'CMD_LABEL':"SBML mapping", 'DESC':"Map an SBML file (all the entities IDs) to a reference database (specified in config.txt)."},
                      'get_medium':{'CMD_LABEL':"Get medium", 'DESC':"Display the studied species specified growth medium."},
                      'set_medium':{'CMD_LABEL':"Set medium", 'DESC':"Set the growth medium for the studied species."},
                      'del_medium':{'CMD_LABEL':"Delete medium", 'DESC':"Delete the growth medium for the studied species."},
                      'get_compart':{'CMD_LABEL':"Get compartments", 'DESC':"Display all the compartments of the metabolic network."},
                      'del_compart':{'CMD_LABEL':"Delete compartment", 'DESC':"Remove a compartment from the metabolic network."},
                      'change_compart':{'CMD_LABEL':"Change compartment", 'DESC':"Modify a compartment in the metabolic network."},
                      'report':{'CMD_LABEL':"Report", 'DESC':"Generate reports on the metabolic network reconstruction."},
                      'set_fba':{'CMD_LABEL':"Set FBA", 'DESC':"Set the biomass reaction to run flux balance analysis on the network."},
                      'summary':{'CMD_LABEL':"Test FBA", 'DESC':"Run flux balance analysis on the network."},
                      'menecheck':{'CMD_LABEL':"Menecheck", 'DESC':"Run topological analysis on the network."},
                      'shogen':{'CMD_LABEL':"Shogen", 'DESC':"Run Shogen : find shortest genome segments that regulate metabolic pathways."},
                      'wiki_pages':{'CMD_LABEL':"Create Wiki pages", 'DESC':"Create Wiki pages to display the metabolic network."},
                      'wiki_run':{'CMD_LABEL':"Run Wiki", 'DESC':"Create a Docker container for the Wiki."},
                      'wiki_init':{'CMD_LABEL':"Wiki initialization", 'DESC':"Send data on the metabolic network to the Docker container to fill in the Wiki."},
                      'send_all_page':{'CMD_LABEL':"Send all pages to Wiki", 'DESC':"Send all the generated pages on the metabolic network to the Wiki."},
                      'tsv':{'CMD_LABEL':"PADMet to tsv",'DESC':"Convert a PADMet (.padmet) file to a tsv files (for Askomics)."}
                      }
    current_cmd_dict = cmd_label_dict.get(cmd)
    if current_cmd_dict:
        cmd_label = current_cmd_dict["CMD_LABEL"]
        desc = current_cmd_dict["DESC"]
    else:
        cmd_label, desc = None, None
    return(cmd_label,desc)



if __name__ == "__main__":
    main()