# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 16:34:38 2018

@author: maite
toy for ptg tests
Benchmark:
T_lutea
Athalia, Chlamy, e_sill, synecho
Pour chaque model:
    1./ Lire output pantograph, récupérer N reactions et faire un sous sbml
    2./ Récuéprer la lsit des genes assoc de ces N reactions, lire le fasta Study et faire un sous fasta avec X% de bruit
    3./ Lire le sbml du Model, recup la liste des genes assoc et faire
avec du bruit X: % de genes totaux or cela à gader dans le fasta


Tiso: Tiso_gene_9894 / R_ARYLDIALKYL__45__PHOSPHATASE__45__RXN
Ecto: Ec-27_006880
EC-27_006


"""

import libsbml
from Bio import SeqIO


def main():
    output_pantograph = "/home/maite/Forge/docker/aureme_workspace/Sjap/networks/output_orthology_based_reconstruction/pantograph/output_pantograph_Nannochloropsis_salina.sbml"
    study_fasta = "/home/maite/Forge/docker/aureme_workspace/Sjap/orthology_based_reconstruction/Nannochloropsis_salina/FAA_model.faa"
    model_fasta = "/home/maite/Forge/docker/aureme_workspace/Sjap/orthology_based_reconstruction/Nannochloropsis_salina/FAA_model.faa"


reader = libsbml.SBMLReader()
document = reader.readSBML(output_pantograph)
for i in range(document.getNumErrors()):
    print (document.getError(i).getMessage())
model = document.getModel()
listOfReactions = model.getListOfReactions()


def sub_sbml(model, list_rxn):
    return new_model

def sub_fasta(fasta, list_genes, noise):
    return new_fasta