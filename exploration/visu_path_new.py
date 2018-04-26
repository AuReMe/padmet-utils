# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 11:34:54 2018

@author: maite
"""

from padmet.classes import PadmetRef, PadmetSpec
import networkx as nx
import matplotlib.pylab as plt
import csv



def main():
    ulva_bact = PadmetSpec("/home/maite/Documents/data/Ulva/ulva_bacteria.padmet")
    padmetRef = PadmetRef("/home/maite/Forge/docker/aureme_workspace/ulva/database/metacyc_21.5_enhanced.padmet")
    output_folder = "/home/maite/Documents/data/Ulva/pwy_analysis/"
    pwy_csv = output_folder+"pathways.tsv"
    png_folder = output_folder+"pwys/"    
    #k=pwy_id, v = {'total':[rxn_ids], 'ulva':[rxn_ids], 'bact_1':[rxn_ids], 'bact_2':[rxn_ids]}
    pwy_dict = {}
    all_pwy = set([rlt.id_out for rlt in ulva_bact.getAllRelation() if rlt.type == "is_in_pathway" and ulva_bact.dicOfNode[rlt.id_in].type == "reaction"])
    
    all_src = set()
    count = 0
    print("Extracting data")
    for pwy_id in list(all_pwy):
        count += 1
        print("%s/%s %s" %(count, len(all_pwy), pwy_id))
        #pwy_id = 'PWY-5664'
        pwy_cname = padmetRef.dicOfNode[pwy_id].misc.get("COMMON-NAME",['UNKNOWN'])[0]
        pwy_dict[pwy_id] = {'name':pwy_cname}
        total_reactions = set([rlt.id_in for rlt in padmetRef.dicOfRelationOut[pwy_id] if rlt.type == "is_in_pathway" and padmetRef.dicOfNode[rlt.id_in].type == "reaction"])
        pwy_dict[pwy_id]['total'] = total_reactions
        pwy_dict[pwy_id].update({'ULVA':set(),'MS6_BACTERIA':set(),'MS2_BACTERIA':set()})
        for rxn_id in [rxn_id for rxn_id in total_reactions if rxn_id in ulva_bact.dicOfNode.keys()]:
            src = set([ulva_bact.dicOfNode[rlt.id_out].misc["SOURCE"][0] for rlt in ulva_bact.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"])
            all_src.update(src)
            if 'MS6_BACTERIA' in src:
                pwy_dict[pwy_id]['MS6_BACTERIA'].add(rxn_id)
            if 'MS2_BACTERIA' in src:
                pwy_dict[pwy_id]['MS2_BACTERIA'].add(rxn_id)
            if src.difference(set(['MS6_BACTERIA', 'MS2_BACTERIA'])):
                pwy_dict[pwy_id]['ULVA'].add(rxn_id)
    
    print("Creating pngs")
    count = 0
    for pwy_id, pwy_data in pwy_dict.items():
        count += 1
        print("%s/%s %s" %(count, len(all_pwy), pwy_id))
        #pwy_id = 'PWY-5664'
        #pwy_data =  pwy_dict[pwy_id]
        plt = visu_path(pwy_data, padmetRef)
        png_path = png_folder+"pathway_"+pwy_id
        save_plot(plt, png_path)

    print("Creating CSV")
    count = 0
    with open(pwy_csv, 'w') as csvfile:
        fieldnames = ['pathway_id','pathway_name',
                      'Nb_all_reactions','Nb_reactions_in_ULVA','Ratio_ULVA',
                      'Nb_reactions_in_ULVA_MS2_BACTERIA', 'Ratio_ULVA_MS2_BACTERIA',
                      'Nb_reactions_in_ULVA_MS6_BACTERIA', 'Ratio_ULVA_MS6_BACTERIA',
                      'Nb_reactions_in_ULVA_MS2_MS6_BACTERIA','Ratio_ULVA_MS2_MS6_BACTERIA',
                      'All_reactions', 'Reactions_in_ULVA',
                      'Reactions_in_ULVA_MS2_BACTERIA', 'Reactions_in_ULVA_MS6_BACTERIA',
                      'Reactions_in_ULVA_MS2_MS6_BACTERIA']

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
                
        for pwy_id, pwy_data in pwy_dict.items():
            count += 1
            print("%s/%s %s" %(count, len(all_pwy), pwy_id))
            #pwy_id = 'PWY-5664'
            #pwy_data =  pwy_dict[pwy_id]
            ratio_ulva = round(float(len(pwy_data['ULVA']))/float(len(pwy_data['total'])),2)
            rxn_ulva_ms2 = pwy_data['ULVA'].union(pwy_data['MS2_BACTERIA'])
            ratio_ulva_ms2 = round(float(len(rxn_ulva_ms2))/float(len(pwy_data['total'])),2)
            rxn_ulva_ms6 = pwy_data['ULVA'].union(pwy_data['MS6_BACTERIA'])
            ratio_ulva_ms6 = round(float(len(rxn_ulva_ms6))/float(len(pwy_data['total'])),2)
            rxn_ulva_ms2_ms6 = rxn_ulva_ms2.union(rxn_ulva_ms6)
            ratio_ulva_ms2_ms6 = round(float(len(rxn_ulva_ms2_ms6))/float(len(pwy_data['total'])),2)
            

            dict_write = {'pathway_id':pwy_id, 'pathway_name':pwy_data['name']}
            dict_write.update({'All_reactions':';'.join(pwy_data['total']),'Nb_all_reactions':len(pwy_data['total'])})
            dict_write.update({'Reactions_in_ULVA':';'.join(pwy_data['ULVA']),'Nb_reactions_in_ULVA':len(pwy_data['ULVA']),'Ratio_ULVA':ratio_ulva})
            dict_write.update({'Reactions_in_ULVA_MS2_BACTERIA':';'.join(rxn_ulva_ms2),'Nb_reactions_in_ULVA_MS2_BACTERIA':len(rxn_ulva_ms2),'Ratio_ULVA_MS2_BACTERIA':ratio_ulva_ms2})
            dict_write.update({'Reactions_in_ULVA_MS6_BACTERIA':';'.join(rxn_ulva_ms6),'Nb_reactions_in_ULVA_MS6_BACTERIA':len(rxn_ulva_ms6),'Ratio_ULVA_MS6_BACTERIA': ratio_ulva_ms6})
            dict_write.update({'Reactions_in_ULVA_MS2_MS6_BACTERIA':';'.join(rxn_ulva_ms2_ms6),'Nb_reactions_in_ULVA_MS2_MS6_BACTERIA':len(rxn_ulva_ms2_ms6),'Ratio_ULVA_MS2_MS6_BACTERIA': ratio_ulva_ms2_ms6})
            writer.writerow(dict_write)            

def visu_path(pwy_data, padmetRef):
    DG=nx.DiGraph()
    nodes_attrib = {}
    for reaction_id in pwy_data['total']:
        # Reaction colors
        nodes_attrib[reaction_id] = {}
        nodes_attrib[reaction_id]['color'] = "lightgreen"
        nodes_attrib[reaction_id]['type'] = "reaction"
        # Reactants & products for each reaction
        reactants = [rlt.id_out
            for rlt in padmetRef.dicOfRelationIn.get(reaction_id, None)
                if rlt.type == "consumes"]
        products = [rlt.id_out
            for rlt in padmetRef.dicOfRelationIn.get(reaction_id, None)
                if rlt.type == "produces"]
        reac_rev = padmetRef.dicOfNode[reaction_id].misc['DIRECTION'][0]
        for reac in reactants:
            if reac not in nodes_attrib.keys():
                nodes_attrib[reac] = {}
                nodes_attrib[reac]['type'] = "metabolite"
            DG.add_edge(reac, reaction_id)
            if reac_rev == "REVERSIBLE":
                DG.add_edge(reaction_id, reac)

        for prod in products:
            if prod not in nodes_attrib.keys():
                nodes_attrib[prod] = {}
                nodes_attrib[prod]['type'] = "metabolite"
            DG.add_edge(reaction_id, prod)
            if reac_rev == "REVERSIBLE":
                DG.add_edge(prod, reaction_id)

    # https://networkx.github.io/documentation/latest/reference/generated/networkx.drawing.nx_pylab.draw_networkx.html
    # apt-get install graphviz graphviz-dev (python-pygraphviz)
    # pip install pygraphviz
    plt.figure(figsize=(20,20)) 
    nx.draw_networkx_nodes(DG,
                           nx.graphviz_layout(DG, prog='neato'),
                           nodelist=[k for k,v in nodes_attrib.iteritems() if v['type'] == "reaction"],
                           node_color="lightgreen",
                           node_size=2000,
                           node_shape='s',
                       alpha=0.2)
    nx.draw_networkx_nodes(DG,
                           nx.graphviz_layout(DG, prog='neato'),
                           nodelist=[k for k,v in nodes_attrib.iteritems() if v['type'] == "reaction" and k in pwy_data['ULVA']], 
                           node_color="green",
                           node_size=1500,
                           node_shape='s',
                           label="In Ulva",
                       alpha=0.8)
    nx.draw_networkx_nodes(DG,
                           nx.graphviz_layout(DG, prog='neato'),
                           nodelist=[k for k,v in nodes_attrib.iteritems() if v['type'] == "reaction" and k in pwy_data['MS6_BACTERIA']], 
                           node_color="red",
                           node_size=1000,
                           node_shape='s',
                           label="In MS6_BACTERIA",
                       alpha=0.8)
    nx.draw_networkx_nodes(DG,
                           nx.graphviz_layout(DG, prog='neato'),
                           nodelist=[k for k,v in nodes_attrib.iteritems() if v['type'] == "reaction" and k in pwy_data['MS2_BACTERIA']], 
                           node_color="blue",
                           node_size=500,
                           node_shape='s',
                           label="In MS2_BACTERIA",
                       alpha=0.8)
    nx.draw_networkx_nodes(DG,
                           nx.graphviz_layout(DG, prog='neato'),
                           nodelist=[k for k,v in nodes_attrib.iteritems() if v['type'] == "metabolite"],
                           node_color="yellow",
                           node_size=500,
                           node_shape='o',
                           label="Metabolite",
                       alpha=0.8)
    nx.draw_networkx_edges(DG,
                           nx.graphviz_layout(DG, prog='neato'),
                           edge_color="black",
                           alpha=0.5,
                           width=2.0,
                           arrow=True)
    nx.draw_networkx_labels(DG,
                           nx.graphviz_layout(DG, prog='neato'),
                           font_size=15)
    plt.axis('off')
    plt.legend(scatterpoints=1, markerscale=0.6)
    return plt

def save_plot(plot, filepath):
    """Saves plot in multiple formats

    :param arg1: plot object
    :param arg2: filename with its path
    :type arg1: <plot object>
    :type arg2: <str>

    """
    plot.savefig(filepath + ".png",
                 dpi=144, format='png', frameon=True)
    plot.close()
    #plot.savefig(filepath + ".pdf",
    #              dpi=144, format='pdf')
                 
if __name__ == "__main__":
    main()
