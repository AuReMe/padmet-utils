#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of padmet.

padmet is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

padmet is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with padmet. If not, see <http://www.gnu.org/licenses/>.

@author: Meziane AITE, meziane.aite@inria.fr
Description:
Contains all necessary functions to generate wikiPages from a padmet file and update 
a wiki online. Require WikiManager module (with wikiMate,Vendor)

usage:
    wikiGenerator.py --padmetSpec=FILE --output=DIR --model_id=STR --model_name=STR [--padmetRef=FILE] -v

options:
    -h --help     Show help.
    --reaction_file=FILE    pathname of the file containing the reactions id, 1/line 
    --padmetRef=FILE    pathname of the padmet representing the database.
    --output=FILE    pathname of the file with line = pathway id, all reactions id, reactions ids from reaction file, ratio. sep = "\t"
"""
from padmet.padmetRef import PadmetRef
from padmet.padmetSpec import PadmetSpec
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
    if args["--padmetRef"]:
        padmetRef = PadmetRef(args["--padmetRef"])
    else:
        padmetRef = None
    verbose = args["-v"]
    model_id, model_name = args["--model_id"], args["--model_name"]
    wiki_folder = args["--output"]
    if not wiki_folder.endswith("/"): wiki_folder += "/"    

    createDirectory()
    all_rxns = [node for node in padmetSpec.dicOfNode.values() if node.type == "reaction"]
    all_genes = [node for node in padmetSpec.dicOfNode.values() if node.type == "gene"]
    for rxn_node in all_rxns:
        create_biological_page("Reaction", rxn_node, wiki_folder+"reactions/")
    for gene_node in all_genes:    
        create_biological_page("Gene", gene_node, wiki_folder+"genes/")
    all_pwys = [node for (node_id, node) in padmetSpec.dicOfNode.iteritems() if node_id in total_pwy_id]
    all_cpds = [node for (node_id, node) in padmetSpec.dicOfNode.iteritems() if node_id in total_cpd_id]
    for pwy_node in all_pwys:
        create_biological_page("Pathway", pwy_node, wiki_folder+"pathways/")
    for cpd_node in all_cpds:
        create_biological_page("Metabolite", cpd_node, wiki_folder+"metabolites/")

    create_navigation_page(wiki_folder+"/navigation/")
    create_venn()
    create_main(model_id, model_name)

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
    all_categories = ["orthology","annotation","gap-filling","manual"]
    for category in all_categories:
        categories_dict[category] = set()
    for rxn_id, rxn_src_dict in full_sources_dict.items():
        for category in rxn_src_dict.keys():
            categories_dict[category].add(rxn_id)
    
    labels = get_labels(categories_dict.values())
    fig, ax = venn4(labels, names=categories_dict.keys())
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
    reconstruct_summary = {"ANNOTATION":0,"ORTHOLOGY":{},"MANUAL":0,"GAP-FILLING":0}
    for rec_node in [node for node in padmetSpec.dicOfNode.values() if node.type == "reconstructionData"]:
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
        for k,v in reconstruct_summary["ORTHOLOGY"].items():
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
    sideBarData = ["* navigation","** mainpage|mainpage-description","** randompage-url|randompage","** Special:ListFiles|Files","* Metabolic network components"]

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
    dataInArray = ["{{#ask: [[Category:Pathway]]","| ?common name","| ?reaction found","| ?reaction not found","| ?completion rate","}}"]
    with open(fileName,'w') as f:
        for line in dataInArray:
            f.write(line+"\n")

    category = "Metabolite"
    sideBarData.append("** Category:"+category+"|"+category)
    fileName = output_folder+"Category:Metabolite"
    if verbose: print("Category: %s" %category)
    dataInArray = ["{{#ask: [[Category:Metabolite]]","| ?common name","| ?consumed by","| ?produced by","| ?consumed or produced by","}}"]
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

    fileName = output_folder + page_node.id.replace("/",".")
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
    for k,v in page_node.misc.iteritems():
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
                dataInArray.append("* [["+gene_id+"]]")
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
                    #TODO: not only for pantograph... 
                    #if reconstruction source start with output_pantograph_...
                    #reconverting it to a more readable format: [[PANTOGRAPH]]-[['TEMPLATE']]
                    if src.startswith("OUTPUT_PANTOGRAPH_"):
                        src = "[[pantograph]]-[["+src.replace("OUTPUT_PANTOGRAPH_","").lower()+"]]"
                    dataInArray.append("** "+src)
                    if assignment:
                        dataInArray.append("***"+assignment)
        
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
            if comment in ["added to manage seeds from boundary to extracellular compartment", "added to manage seeds from extracellular to cytosol compartment"]:
                source = "MANUAL-IMPORT_FROM_MEDIUM"
            if source:
                source = source.lower()
                if source.startswith("output_pantograph_"):
                    source = source.replace("output_pantograph_","")
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
        for category, category_data in src_data.items():
            dataInArray.append("* Category: [["+category+"]]:")
            for source, source_data in category_data.items():
                dataInArray.append("** Source: [["+source+"]]:")
                if source_data["tool"]:
                    dataInArray.append("*** Tool: [["+source_data["tool"]+"]]:")
                    if source_data["comment"]:
                        dataInArray.append("**** Comment: [["+source_data["comment"]+"]]:")
                elif source_data["comment"]:
                    dataInArray.append("*** Comment: [["+source_data["comment"]+"]]:")

    elif category == "Gene":
        dataInArray.append("== Reactions associated ==")
        linked_rlt = [rlt for rlt in padmetSpec.dicOfRelationOut[page_node.id] if rlt.type == "is_linked_to"]
        if linked_rlt:
            pwy_assoc = set()
            add_property(properties, "reaction associated", [rlt.id_in for rlt in linked_rlt])
            for rlt in linked_rlt:
                rxn_id = rlt.id_in
                [pwy_assoc.add(rlt_pwy.id_out) for rlt_pwy in padmetSpec.dicOfRelationIn[rxn_id] if rlt_pwy.type == "is_in_pathway"]
                dataInArray.append("* [["+rxn_id+"]]")
                sources = rlt.misc["SOURCE:ASSIGNMENT"]
                for src_data in sources:
                    try:
                        src, assignment = src_data.split(":")
                    except ValueError:
                        src = src_data
                        assignment = None
                    src = src.lower()
                    if src.startswith("output_pantograph_"):
                        src = "[[pantograph]]-[["+src.replace("output_pantograph_","")+"]]"
                    dataInArray.append("** "+src)
                    if assignment:
                        assignment = assignment.lower()
                        dataInArray.append("***"+assignment)
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
        add_property(properties, "reaction not found", [len(reactionsMissing)])
        add_property(properties, "completion rate", [pwy_ratio])
        for rxn_id in reactionsFound:
            dataInArray.append("* [["+rxn_id+"]]")
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
            add_property(properties, "consumed or produced by",rxn_cp["cp"])
            for rxn_id in rxn_cp["cp"]:
                dataInArray.append("* [["+rxn_id+"]]")

            
    dataInArray.append('== External links  ==')
    try:
        xref_node = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[page_node.id] if rlt.type == "has_xref"][0]
        for db, ids in xref_node.misc.items():
            xrefLink(dataInArray, db, ids)
    except (IndexError, KeyError) as e:
        pass
    
    dataInArray.extend(properties)
    with open(fileName,'w') as f:
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
    
    @type data: list[Iterable]    
    @rtype: dict[str, str]
    input
      data: data to get label for
      fill: ["number"|"logic"|"percent"]
    return
      labels: a dict of labels for different sets
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
        
    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)
    
    input
      labels: a label dict where keys are identified via binary codes ('0001', '0010', '0100', ...),
              hence a valid set could look like: {'0001': 'text 1', '0010': 'text 2', '0100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi
    return
      pyplot Figure and AxesSubplot object
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

if __name__ == "__main__":
    main()