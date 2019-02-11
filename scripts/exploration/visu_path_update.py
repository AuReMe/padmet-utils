# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 11:34:54 2018

@author: maite
"""

import csv
import matplotlib.pylab as plt
import networkx as nx
import os
from matplotlib import colors as mcolors

from padmet.classes import PadmetRef, PadmetSpec
from networkx.drawing.nx_agraph import graphviz_layout

def main():
    global src_tags, padmetRef, padmetSpec
    args = {"--padmetSpec": "/home/maite/Forge/docker/aureme_workspace/chondrus/networks/chondrus_1.3.padmet",
            "--padmetRef" : "/home/maite/Aureme/aureme/database/BIOCYC/METACYC/22.0/metacyc_22.0_enhanced.padmet",
            "--output_folder" : "/home/maite/Forge/docker/aureme_workspace/chondrus/analysis/report/1.3_path_visu/",
            "--select_from": 'CATEGORY',
            "--pwys":None,
            "--rxns":None,
            "--split":False}
    #"--pwys":"['late_SSR','early_SSR','PWY-81', 'IDNCAT-PWY', 'PWY-6364', 'PWY-5703', 'PWY-5915']

    #src_dict= {orginal_src_id:tag}
    padmetSpec = PadmetSpec(args["--padmetSpec"])
    padmetRef = PadmetRef(args["--padmetRef"])
    output_folder = args["--output_folder"]
    #pwy_csv = output_folder+"pathways.tsv"
    png_folder = output_folder+"pwys/"   
    select_from=args["--select_from"].upper()
    verbose = True
    split = args["--split"]

    if select_from not in ["SOURCE","TOOL","CATEGORY"]:
        print("please select from [SOURCE,TOOL,CATEGORY] only")
        exit()

    #k=pwy_id, v = {'total':[rxn_ids], 'ulva':[rxn_ids], 'bact_1':[rxn_ids], 'bact_2':[rxn_ids]}
    pwy_dict = {}
    src_tags = set()
    if args["--pwys"]:
        all_pwy = eval(args["--pwys"])
        if verbose:
            print("%s pathways selected: %s"%(len(all_pwy), all_pwy))
    else:
        if verbose:
            print("all pathways mode")
        all_pwy = set([rlt.id_out for rlt in padmetSpec.getAllRelation() if rlt.type == "is_in_pathway" and padmetSpec.dicOfNode[rlt.id_in].type == "reaction"])
    
    count = 0
    print("Extracting data")
    for pwy_id in list(all_pwy):
        count += 1
        if verbose:
            print("%s/%s %s" %(count, len(all_pwy), pwy_id))
        #pwy_id = 'late_SSR'
        pwy_cname = padmetSpec.dicOfNode[pwy_id].misc.get("COMMON-NAME",[pwy_id])[0]
        pwy_dict[pwy_id] = {'name':pwy_cname}
        try:
            total_reactions = set([rlt.id_in for rlt in padmetRef.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway" and padmetRef.dicOfNode[rlt.id_in].type == "reaction"])
        except KeyError:
            total_reactions = set([rlt.id_in for rlt in padmetSpec.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway" and padmetSpec.dicOfNode[rlt.id_in].type == "reaction"])            
        pwy_dict[pwy_id]['total'] = total_reactions
        for rxn_id in [rxn_id for rxn_id in total_reactions if rxn_id in padmetSpec.dicOfNode.keys()]:
            rxn_srcs = set([padmetSpec.dicOfNode[rlt.id_out].misc[select_from][0] for rlt in padmetSpec.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"])
            if rxn_srcs:
                for src in rxn_srcs:
                    src_tags.add(src)
                    try:
                        pwy_dict[pwy_id][src].add(rxn_id)
                    except KeyError:
                        pwy_dict[pwy_id][src] = set([rxn_id])
            else:
                try:
                    pwy_dict[pwy_id]['unknown_source'].add(rxn_id)
                except KeyError:
                    pwy_dict[pwy_id]['unknown_source'] = set([rxn_id])
                    
    if verbose:
        print("Creating pngs")
    count = 0
    init_attribs = set_init_attribs()
    if split:
        for pwy_id, pwy_data in pwy_dict.items():
            count += 1
            if verbose:
                print("%s/%s %s" %(count, len(all_pwy), pwy_id))
            DG=nx.DiGraph()
            nodes_attrib = dict()
            pwy_name = pwy_data['name']
            pwy_data.pop('name')
            for src, reaction_ids in pwy_data.items():
                nodes_attrib[src] = {'reaction_ids':set()}
                if src == 'total':
                    nodes_attrib[src]['metabolite_ids'] = set()
                    for reaction_id in reaction_ids:
                        # Reactants & products for each reaction
                        if padmetRef and reaction_id in padmetRef.dicOfNode.keys():
                            padmet = padmetRef
                        else:
                            padmet = padmetSpec
                        nodes_attrib['total']['reaction_ids'].add(reaction_id)
                        reactants = ["%s:%s[%s]"%(rlt.misc["STOICHIOMETRY"][0].replace(".0",""), rlt.id_out, rlt.misc["COMPARTMENT"][0])
                            for rlt in padmet.dicOfRelationIn.get(reaction_id, None)
                                if rlt.type == "consumes"]
                        products = ["%s:%s[%s]"%(rlt.misc["STOICHIOMETRY"][0].replace(".0",""), rlt.id_out, rlt.misc["COMPARTMENT"][0])
                            for rlt in padmet.dicOfRelationIn.get(reaction_id, None)
                                if rlt.type == "produces"]
                        reac_rev = padmet.dicOfNode[reaction_id].misc['DIRECTION'][0]
                        for reac in reactants:
                            nodes_attrib['total']['metabolite_ids'].add(reac)
                            DG.add_edge(reac, reaction_id)
                            if reac_rev == "REVERSIBLE":
                                DG.add_edge(reaction_id, reac)
    
                        for prod in products:
                            nodes_attrib['total']['metabolite_ids'].add(prod)
                            DG.add_edge(reaction_id, prod)
                            if reac_rev == "REVERSIBLE":
                                DG.add_edge(prod, reaction_id)
                else:
                    for reaction_id in reaction_ids:
                        # Reactants & products for each reaction
                        if padmetRef and reaction_id in padmetRef.dicOfNode.keys():
                            padmet = padmetRef
                        else:
                            padmet = padmetSpec
                        nodes_attrib[src]['reaction_ids'].add(reaction_id)
    
            plt.figure(figsize=(60,60))
            nx.draw_networkx_nodes(DG,
                                   graphviz_layout(DG, prog='neato'),
                                   nodelist=nodes_attrib['total']['metabolite_ids'],
                                   node_color=init_attribs['total']['metabolite_attribs']['color'],
                                   node_size=init_attribs['total']['metabolite_attribs']['size'],
                                   node_shape='o',
                                   label="Metabolite",
                               alpha=0.8)
            for src_tag, attrib_data in sorted(init_attribs.items(), key=lambda x: x[1]['reaction_attribs']['size'], reverse=True):
                if src_tag == 'total':
                    nx.draw_networkx_nodes(DG,
                                           graphviz_layout(DG, prog='neato'),
                                           nodelist=nodes_attrib[src_tag]['reaction_ids'],
                                           node_color=attrib_data['reaction_attribs']['color'],
                                           node_size=attrib_data['reaction_attribs']['size'],
                                           node_shape='s',
                                           label="Total",
                                       alpha=0.2)
                elif src_tag in nodes_attrib.keys():
                    nx.draw_networkx_nodes(DG,
                                           graphviz_layout(DG, prog='neato'),
                                           nodelist=nodes_attrib[src_tag]['reaction_ids'],
                                           node_color=attrib_data['reaction_attribs']['color'],
                                           node_size=attrib_data['reaction_attribs']['size'],
                                           node_shape='s',
                                           label=src_tag,
                                       alpha=0.8)
            nx.draw_networkx_edges(DG,
                                   graphviz_layout(DG, prog='neato'),
                                   edge_color="black",
                                   alpha=0.5,
                                   width=2.0,
                                   arrow=True)
    
            
            nx.draw_networkx_labels(DG,
                                   graphviz_layout(DG, prog='neato'),
                                   font_size=15)
            plt.axis('off')
            plt.legend(title=pwy_name, scatterpoints=1, markerscale=0.1)
            png_path = "%s%s"%(png_folder, pwy_id)
            save_plot(plt, png_path)
    else:
        DG=nx.DiGraph()
        for pwy_id, pwy_data in pwy_dict.items():
            count += 1
            if verbose:
                print("%s/%s %s" %(count, len(all_pwy), pwy_id))
            nodes_attrib = dict()
            pwy_name = pwy_data['name']
            pwy_data.pop('name')
            for src, reaction_ids in pwy_data.items():
                nodes_attrib[src] = {'reaction_ids':set()}
                if src == 'total':
                    nodes_attrib[src]['metabolite_ids'] = set()
                    for reaction_id in reaction_ids:
                        # Reactants & products for each reaction
                        if padmetRef and reaction_id in padmetRef.dicOfNode.keys():
                            padmet = padmetRef
                        else:
                            padmet = padmetSpec
                        nodes_attrib['total']['reaction_ids'].add(reaction_id)
                        reactants = ["%s:%s[%s]"%(rlt.misc["STOICHIOMETRY"][0].replace(".0",""), rlt.id_out, rlt.misc["COMPARTMENT"][0])
                            for rlt in padmet.dicOfRelationIn.get(reaction_id, None)
                                if rlt.type == "consumes"]
                        products = ["%s:%s[%s]"%(rlt.misc["STOICHIOMETRY"][0].replace(".0",""), rlt.id_out, rlt.misc["COMPARTMENT"][0])
                            for rlt in padmet.dicOfRelationIn.get(reaction_id, None)
                                if rlt.type == "produces"]
                        reac_rev = padmet.dicOfNode[reaction_id].misc['DIRECTION'][0]
                        for reac in reactants:
                            nodes_attrib['total']['metabolite_ids'].add(reac)
                            DG.add_edge(reac, reaction_id)
                            if reac_rev == "REVERSIBLE":
                                DG.add_edge(reaction_id, reac)
    
                        for prod in products:
                            nodes_attrib['total']['metabolite_ids'].add(prod)
                            DG.add_edge(reaction_id, prod)
                            if reac_rev == "REVERSIBLE":
                                DG.add_edge(prod, reaction_id)
                else:
                    for reaction_id in reaction_ids:
                        # Reactants & products for each reaction
                        if padmetRef and reaction_id in padmetRef.dicOfNode.keys():
                            padmet = padmetRef
                        else:
                            padmet = padmetSpec
                        nodes_attrib[src]['reaction_ids'].add(reaction_id)
    
        plt.figure(figsize=(100,100))
        nx.draw_networkx_nodes(DG,
                               graphviz_layout(DG, prog='neato'),
                               nodelist=nodes_attrib['total']['metabolite_ids'],
                               node_color=init_attribs['total']['metabolite_attribs']['color'],
                               node_size=init_attribs['total']['metabolite_attribs']['size'],
                               node_shape='o',
                               label="Metabolite",
                           alpha=0.8)
        for src_tag, attrib_data in sorted(init_attribs.items(), key=lambda x: x[1]['reaction_attribs']['size'], reverse=True):
            if src_tag == 'total':
                nx.draw_networkx_nodes(DG,
                                       graphviz_layout(DG, prog='neato'),
                                       nodelist=nodes_attrib[src_tag]['reaction_ids'],
                                       node_color=attrib_data['reaction_attribs']['color'],
                                       node_size=attrib_data['reaction_attribs']['size'],
                                       node_shape='s',
                                       label="Total",
                                   alpha=0.2)
            elif src_tag in nodes_attrib.keys():
                nx.draw_networkx_nodes(DG,
                                       graphviz_layout(DG, prog='neato'),
                                       nodelist=nodes_attrib[src_tag]['reaction_ids'],
                                       node_color=attrib_data['reaction_attribs']['color'],
                                       node_size=attrib_data['reaction_attribs']['size'],
                                       node_shape='s',
                                       label=src_tag,
                                   alpha=0.8)
        nx.draw_networkx_edges(DG,
                               graphviz_layout(DG, prog='neato'),
                               edge_color="black",
                               alpha=0.5,
                               width=2.0,
                               arrow=True)

        
        nx.draw_networkx_labels(DG,
                               graphviz_layout(DG, prog='neato'),
                               font_size=15)
        plt.axis('off')
        plt.legend(title="FULL NETWORK", scatterpoints=1, markerscale=0.1)
        png_path = "%s%s"%(png_folder, "all_pwys")
        save_plot(plt, png_path)

def set_init_attribs():
    node_metabolite_color = mcolors.cnames['yellow']
    node_metabolite_size = 500
    node_reaction_base_color = mcolors.cnames['lightgreen']
    node_reaction_minal_size = 2000
    node_reaction_base_size = node_reaction_minal_size + 1000*len(src_tags)

    colors = [col for col in mcolors.cnames.values() if col not in [node_metabolite_color, node_reaction_base_color]]
    init_attribs = {}
    init_attribs['total'] = {'metabolite_attribs': {'color':node_metabolite_color, 'size':node_metabolite_size},
                        'reaction_attribs': {'color':node_reaction_base_color,'size':node_reaction_base_size}}
    node_reaction_size = node_reaction_base_size - 1000
    
    for index, src_tag in enumerate(src_tags):
        node_reaction_color = colors[index]
        init_attribs[src_tag] = {'reaction_attribs': {'color':node_reaction_color,'size':node_reaction_size}}
        node_reaction_size -= 1000
    return(init_attribs)


def save_plot(plot, filepath):
    """Saves plot in multiple formats

    :param arg1: plot object
    :param arg2: filename with its path
    :type arg1: <plot object>
    :type arg2: <str>

    """
    plot.savefig(filepath + ".svg",
                 dpi=144, format='svg', frameon=True)
    plot.close()
    #plot.savefig(filepath + ".pdf",
    #              dpi=144, format='pdf')
                 
if __name__ == "__main__":
    main()
