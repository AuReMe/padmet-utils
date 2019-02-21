#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 16:40:04 2018

@author: maite
"""

from padmet.classes import PadmetSpec
import padmet.utils.sbmlPlugin as sp
import os
import docopt
import libsbml


reader = libsbml.SBMLReader()
document = reader.readSBML("/home/maite/Documents/data/recon/773_networks/773_networks/Campylobacter_rectus_RM3267.xml")
for i in range(document.getNumErrors()):
    print(document.getError(i).getMessage())
model = document.getModel()

listOfGenesProduct = [i for i in model.getListOfAllElements() if type(i) == libsbml.ListOfGeneProducts][0]
dict_gene_id_label = dict([(gene.id, gene.label) for gene in listOfGenesProduct])

dict_rxn_genes =  {}
for rxn in model.model.getListOfReactions():
    dict_rxn_genes[rxn.id] = set()
    try:
        geneAssoc = [i for i in rxn.getListOfAllElements() if type(i) == libsbml.GeneProductAssociation][0]
        for geneProduct in [i for i in geneAssoc.getListOfAllElements() if type(i) == libsbml.GeneProductRef]:
            dict_rxn_genes[rxn.id].add(geneProduct.getGeneProduct())
    except IndexError:
        pass

#check if the data are all the same
#get all
        