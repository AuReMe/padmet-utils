#!/usr/bin/python
#-*- coding: utf-8 -*-
"""
This file is part of padmet-utils.

padmet-utils is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

padmet-utils is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with padmet-utils. If not, see <http://www.gnu.org/licenses/>.

@author: Meziane AITE, meziane.aite@inria.fr
Description:
For the a given metabolic network: create for each gene un network without the reactions
associated to this gene. Then process to fba et menescope on each network.

usage:
    fba_all_species.py --sbml=FILE

option:
    -h --help    Show help.
    --sbml=FILE    metabolic network in SBML format.
"""
import os
import shutil
from libsbml import *
from padmet.utils.sbmlPlugin import parseGeneAssoc, parseNotes
from multiprocessing import Pool, cpu_count
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from cobra import *

def main():
    #metabolic network in folder X
    #create folders all_networks and scope
    global original_model, wdir, fba_path, menescope_path
    fba_path = os.path.dirname(os.path.realpath(__file__))+"/fba_all_species.py"
    wdir = "/shared/faecalis_v583/analysis/mutagene/"
    sbml_name = "faecalis_v583.sbml"
    sbml_path = wdir+sbml_name
    """
    if os.path.exists(wdir+"all_networks"):
        shutil.rmtree(wdir+"all_networks")
    os.makedirs(wdir+"all_networks")
    if os.path.exists(wdir+"scopes"):
        shutil.rmtree(wdir+"scopes")
    os.makedirs(wdir+"scopes")
    """
    #for g in genes, get the reactions where g is needed and remove them
    original_model=create_cobra_model_from_sbml_file(sbml_path)
    gene_rxn_dict = dict([(gene.id,[r.id for r in gene.reactions]) for gene in original_model.genes])

    print("creating gene_rxn.txt")
    with open(wdir+"gene_rxn.txt", 'w') as f:
        line = "\t".join(["gene","reactions"])+"\n"
        f.write(line)
        for gene, rxns in list(gene_rxn_dict.items()):
            line = "\t".join([gene, ";".join(rxns)])+"\n"
    
    print("multiprocess all models")
    """
    p = Pool(processes=cpu_count())
    resultats = p.map(konetwork, gene_rxn_dict.items())
    [None for _ in resultats]
    """
    
    print("scope")
    all_sbml_path = [wdir+"all_networks/"+i for i in next(os.walk(wdir+"all_networks"))[2]]
    scope(all_sbml_path[0])


def konetwork(gene_rxn_tuple):
    gene, rxns = gene_rxn_tuple
    temp_model = original_model.copy()
    temp_model_path = wdir + "all_networks/"+"faecalis_v583-"+gene+".sbml"
    temp_model.remove_reactions(rxns, remove_orphans=True)
    write_cobra_model_to_sbml_file(temp_model, temp_model_path)


def scope(sbml_path):
    #sbml_path = all_sbml_path[0]
    fba_cmd = "python "+fba_path+" --sbml="+sbml_path+"> "+wdir+"scopes/"+os.path.split(sbml_path)[-1].repalce(".sbml",".FBAscope") 
    os.system(fba_cmd)
    men_cmd = "menescope.py -d "+sbml_path+" -s "+wdir+"seeds_artefacts.sbml"
    os.system(men_cmd)



if __name__ == '__main__':
    main()
