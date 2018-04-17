# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 10:47:07 2018

@author: maite
Try to find required reactions to produce biomass.
Need a functionnal metabolic network, (flux in biomass rxn)
1./ From model, get list of reactions with flux
2./ Add all reactions
3./ Remove this reactions 1 by 1. If flux not working then select this reaction in final set.
"""
from cobra import *
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from padmet.classes import PadmetSpec
from padmet.utils.sbmlPlugin import convert_from_coded_id

def main():
    repair_file = "/home/maite/Forge/docker/aureme_workspace/chondrus/networks/ESiliculosus_final_alanine.sbml"
    repair_model = create_cobra_model_from_sbml_file(repair_file)
    to_repair_file = "/home/maite/Forge/docker/aureme_workspace/chondrus/networks/draft_gapfilled.sbml"
    to_repair_model = create_cobra_model_from_sbml_file(to_repair_file)
    rxn_flux = get_flux(repair_model)
    #check rxn with flux not in to_repair
    not_in = set(rxn_flux).difference([rxn.id for rxn in to_repair_model.reactions])
    print("nb rxn with flux in model: %s" %len(rxn_flux))
    print("nb rxn not in to_repai: %s" %len(not_in))
    padmet_repair = PadmetSpec("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/ESiliculosus_final_alanine.padmet")
    not_in_with_ga = list()
    for rxn_id in not_in:
        dec, typ, comp =convert_from_coded_id(rxn_id)
        if comp:
            dec += "_"+comp
        ga = [rlt.id_out for rlt in padmet_repair.dicOfRelationIn[dec] if rlt.type == "is_linked_to"]
        if ga:
            not_in_with_ga.append(dec)
    for i in not_in_with_ga:
        print i
    

    with open("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/not_in.txt", 'w') as f:
        for rxn_id in not_in:
            dec, typ, comp =convert_from_coded_id(rxn_id)
            if comp:
                dec += "_"+comp
            line = "\t".join([rxn_id,'//',dec])
            f.write(line+"\n")
    sjap_model = create_cobra_model_from_sbml_file("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/sjap_alanine.sbml")
    not_in_sjap = set(rxn_flux).difference([rxn.id for rxn in sjap_model.reactions])
    rxn_flux_sjap = get_flux(sjap_model)
    inter_sjap_ecto = set(rxn_flux_sjap).intersection(set(rxn_flux))
    
    rxn_to_keep = list()
    for rxn_id in rxn_flux:
        dec, typ, comp =convert_from_coded_id(rxn_id)
        if comp:
            dec += "_"+comp
        print dec,typ,comp
        rxn_to_keep.append(dec)
    for rxn_id in [node.id for node in padmet.dicOfNode.values() if node.type == "reaction"]:
        if rxn_id not in rxn_to_keep:
            padmet.delNode(rxn_id)
    padmet.generateFile("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/repair.padmet")


def get_flux(model):
    flux = flux_analysis.flux_variability_analysis(model)
    rxn_flux = list()
    for index, row in flux.iterrows():
        max_flux = row['maximum']
        min_flux = row['minimum']
        if max_flux != 0 or min_flux != 0:
            rxn_flux.append(index)
    return rxn_flux


ecto = create_cobra_model_from_sbml_file("/home/maite/Forge/docker/aureme_workspace/chondrus/orthology_based_reconstruction/Ectocarpus_siliculosus/metabolic_model.sbml")
chondrus = create_cobra_model_from_sbml_file("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/draft_gapfilled.sbml")
not_in = set(rxn_postifif).difference(set([i.id for i in chondrus.reactions]))
ecto_padmet = PadmetSpec("/home/maite/Forge/docker/aureme_workspace/e_siliculosus/networks/ESiliculosus_final.padmet")
chondrus_padmet = PadmetSpec("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/draft_gapfilled.padmet")
for rxn in not_in:
    dec, typ, comp =convert_from_coded_id(rxn)
    if comp:
        dec += "_"+comp
    print dec,typ,comp
    chondrus_padmet.copyNode(ecto_padmet, dec)
chondrus_padmet.generateFile("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/add_rxn_flux_ecto.padmet")
for rxn in not_in:
    dec, typ, comp =convert_from_coded_id(rxn)
    if comp:
        dec += "_"+comp
    print dec,typ,comp
    chondrus_padmet.delNode(dec)
    

chondrus_2 = create_cobra_model_from_sbml_file("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/add_rxn_flux_ecto.sbml")

chondrus_2_positif = set()
for index, row in flux_analysis.flux_variability_analysis(chondrus_2).iterrows():
    max_flux = row['maximum']
    min_flux = row['minimum']
    if max_flux != 0 or min_flux != 0:
        chondrus_2_positif.add(index)
new_flux = set(chondrus_2_positif).difference(set([i.id for i in chondrus.reactions]))
to_conserv = set()
count = 0
temp = chondrus_2.copy()
for rxn in new_flux:
    count += 1
    print("%s/%s %s" %(count, len(new_flux), rxn))
    temp.reactions.remove(rxn)
    if temp.optimize().objective_value < 0:
        print rxn
        break

ecto_alanine = create_cobra_model_from_sbml_file("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/ecto_test.sbml")
alanine_pos = set()
for index, row in flux_analysis.flux_variability_analysis(ecto_alanine).iterrows():
    max_flux = row['maximum']
    min_flux = row['minimum']
    if max_flux != 0 or min_flux != 0:
        alanine_pos.add(index)
pos_alanine_not_in_chondrus = set(alanine_pos).difference(set([i.id for i in chondrus.reactions]))
ess_alanine = flux_analysis.find_essential_reactions(ecto_alanine)
ess_ala_not_in_chondrus = set([i.id for i in ess_alanine]).difference(set([i.id for i in chondrus.reactions]))

a = create_cobra_model_from_sbml_file("/home/maite/Forge/docker/aureme_workspace/chondrus/networks/test_chondrus.sbml")
a.optimize().objective_value
set([i.id for i in ess_alanine]).difference(set([i.id for i in a.reactions]))