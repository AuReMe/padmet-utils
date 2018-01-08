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
The module sbmlGenerator contains functions to generate sbml files from padmet and txt
usign the libsbml package

usage:
    sbmlGenerator.py --padmet=FILE --output=FILE [--model_id=STR] [--model_name=STR] [--obj_fct=STR] [--sbml_lvl=x] [--sbml_version=y] [--mnx_chem_prop=FILE] [--mnx_chem_xref=FILE] [-v]
    sbmlGenerator.py --compound=FILE --output=FILE [--padmetRef=FILE] [-v]
    sbmlGenerator.py --reaction=FILE --output=FILE --padmetRef=FILE [-v]

option:
    -h --help    Show help.
    --padmet=FILE    pathname of the padmet file to convert into sbml
    --output=FILE    pathanme of the sbml file to generate.
    --obj_fct=STR    id of the reaction objective.
    --sbml_lvl=x    sbml level 2 is sufficient for FBA [default: 2].
    --sbml_version=y    sbml version 1 is sufficient for FBA [default: 1].
    -v   print info.
"""
from padmet.padmetSpec import PadmetSpec
import padmet.sbmlPlugin as sp
import os
import docopt
import cobra
from cobra.io.sbml import write_cobra_model_to_sbml_file
import libsbml
  
#default variables
global def_max_lower_bound, def_max_upper_bound
def_max_upper_bound = 1000
def_max_lower_bound = -1000

def main():
    args = docopt.docopt(__doc__)
    output = args["--output"]
    verbose = args["-v"]
    if args["--padmet"]:
        padmet = PadmetSpec(args["--padmet"])
        if args["--model_id"]:
            model_id = args["--model_id"]
        else:
            model_id = "GEM"
        model_name = args["--model_name"]
        sbml_lvl = int(args["--sbml_lvl"])
        sbml_version = int(args["--sbml_version"])
        obj_fct = args["--obj_fct"]
        mnx_chem_xref = args["--mnx_chem_xref"]
        mnx_chem_prop = args["--mnx_chem_prop"]
        padmet_to_sbml(padmet, output, model_id, model_name, obj_fct, sbml_lvl, sbml_version, mnx_chem_xref, mnx_chem_prop, verbose)
    elif args["--reaction"]:
        with open(args["--reaction"], 'r') as f:
            reactions = set(f.read().splitlines())
        padmetRef = PadmetSpec(args["--padmetRef"])
        reaction_to_sbml(reactions, output, padmetRef, verbose)
    elif args["--compound"]:
        with open(args["--compound"], 'r') as f:
            species_compart = [line.split("\t") for line in f.read().splitlines()]
        compound_to_sbml(species_compart, output, verbose)



def padmet_to_sbml(padmet, output, model_id = None, model_name = None, obj_fct = None, sbml_lvl = 3, sbml_version = 1, mnx_chem_xref = None, mnx_chem_prop = None, verbose = False):
    """
    Convert padmet file to sbml file.
    Specificity: 
    - ids are encoded for sbml using functions sbmlPlugin.convert_to_coded_id
    @param padmet_file: the pathname to the padmet file to convert
    @param output: the pathname to the sbml file to create
    @param obj_fct: the identifier of the objection function, the reaction to test in FBA
    @param sbml_lvl: the sbml level
    @param sbml_version: the sbml version
    @param verbose: print informations
    @type padmet_file, output, verbose: str
    @type sbml_lvl, sbml_version: int
    @return: check return of writeSBMLToFile
    @rtype: int
    """
    #k = xref id, v = mnx_id
    if mnx_chem_xref and mnx_chem_prop:
        dict_mnx_chem_xref = {}
        with open(mnx_chem_xref, 'r') as f:
            for k,v in [line.split("\t")[:2] for line in f.read().splitlines() if not line.startswith("#") and not line.startswith("MNX")]:
                if ":" in k: k = k.split(":")[1]
                try:
                    dict_mnx_chem_xref[v].add(k)
                except KeyError:
                    dict_mnx_chem_xref[v] = set([k])
                    
        #k=mnx_id, v = dict of data
        dict_mnx_chem_prop = {}
        with open(mnx_chem_prop, 'r') as f:
            for line in (line.split("\t") for line in f.read().splitlines() if not line.startswith("#")):
                dict_mnx_chem_prop[line[0]] = {"description":line[1],"formula":line[2],"charge":line[3],"mass":line[4],"inchi":line[5],"smiles":line[6],"source":line[7],"inchikey":line[8]}
    else: dict_mnx_chem_xref, dict_mnx_chem_prop = dict(), dict()
    if not model_id:
        model_id = os.path.splitext(os.path.basename(output))[0]
    model = cobra.Model(id_or_model = model_id, name = model_name)
    rxn_to_test = obj_fct
    if rxn_to_test and rxn_to_test not in padmet.dicOfNode.keys():
        raise KeyError("%s reaction not found in padmet" %rxn_to_test)
    
    species_compart = [(rlt.id_out, rlt.misc["COMPARTMENT"][0]) for rlt in padmet.getAllRelation() 
    if rlt.type in ["consumes","produces"]]
    
    if verbose: print("%s species" %len(species_compart))
    for species_id, compart in species_compart:
        if compart != "C-BOUNDARY":
            sbml_species_id = sp.convert_to_coded_id(species_id, None,compart)
            #from mnx, get charge and formula.
            try:
                mnx_id = [k for k,v in dict_mnx_chem_xref.items() if species_id in v][0]
                species_prop = dict(dict_mnx_chem_prop[mnx_id])
                species_formula = str(species_prop["formula"])
                for k,v in species_prop.items():
                    if k in ["formula", "source", "description"] or v in ["NA",""]:
                        species_prop.pop(k)
            except (IndexError, KeyError) as e:
                print(species_id)
                species_prop = None
                species_formula = None
            species_name = padmet.dicOfNode[species_id].misc.get("COMMON-NAME",[species_id])[0]
            sbml_species = cobra.Metabolite(sbml_species_id, species_formula, species_name, None, compart)
            sbml_species.notes = species_prop
            model.add_metabolites(sbml_species)
    
    all_rxn = [node for node in padmet.dicOfNode.values() if node.type == "reaction"]
    if verbose: print("%s reactions" %len(all_rxn))
    for rxn_node in all_rxn:
        #set rxn id
        #set rxn name
        rxn_id = rxn_node.id
        sbml_rxn_id = sp.convert_to_coded_id(rxn_id)
        rxn_name = rxn_node.misc.get("COMMON-NAME",[sbml_rxn_id])[0]
        sbml_rxn = cobra.Reaction(sbml_rxn_id, rxn_name)
        #set rxn reversibility
        rxn_rev = rxn_node.misc["DIRECTION"][0]
        if rxn_rev == "REVERSIBLE":
            sbml_rxn.lower_bound = def_max_lower_bound
            sbml_rxn.upper_bound = def_max_upper_bound
        elif rxn_rev == "LEFT-TO-RIGHT":
            sbml_rxn.lower_bound = 0
            sbml_rxn.upper_bound = def_max_upper_bound
    
        #get consumed and produced compounds
        consumed = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc["COMPARTMENT"][0]) 
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "consumes" and rlt.misc["COMPARTMENT"][0] != "C-BOUNDARY")
                
        produced = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc["COMPARTMENT"][0]) 
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "produces" and rlt.misc["COMPARTMENT"][0] != "C-BOUNDARY")
    
        metabolitedict = {}
        for species, stoichio, compart in consumed:
            try:
                stoichio = float(stoichio)*-1
            except ValueError:
                stoichio = float(-1.0)
            metabolitedict[sp.convert_to_coded_id(species,None,compart)] = stoichio
        for species, stoichio, compart in produced:
            try:
                stoichio = float(stoichio)
            except ValueError:
                stoichio = float(1.0)
            metabolitedict[sp.convert_to_coded_id(species,None,compart)] = stoichio
        model.add_reaction(sbml_rxn)
        sbml_rxn.add_metabolites(metabolitedict)
    
        #set gene reaction rule
        linked_genes = set([rlt.id_out for rlt in padmet.dicOfRelationIn.get(rxn_id, [])
        if rlt.type == "is_linked_to"])
        if len(linked_genes) != 0:
            gene_reaction_rule = "("+" or ".join(linked_genes)+")"
            sbml_rxn.gene_reaction_rule = gene_reaction_rule
        
        #add subsystem from pathways
        pathways = set([rlt.id_out for rlt in padmet.dicOfRelationIn.get(rxn_id, [])
        if rlt.type == "is_in_pathway"])
        if len(pathways) != 0:
            subsystem = " , ".join(pathways)
            sbml_rxn.subsystem = subsystem        
        if rxn_id == rxn_to_test:
            sbml_rxn.objective_coefficient = 1.0
    
    write_cobra_model_to_sbml_file(model, output, sbml_lvl, sbml_version)





#################################

def reaction_to_sbml(reactions, output, padmetRef, verbose = False):
    """
    convert a list of reactions to sbml format based on a given padmet of reference.
    - ids are encoded for sbml using functions sbmlPlugin.convert_to_coded_id
    @param reactions: list of reactions ids
    @param padmetRef: padmet of reference
    @param output: the pathname to the sbml file to create
    @param sbml_lvl: the sbml level
    @param sbml_version: the sbml version
    @param verbose: print informations
    @type reactions: set
    @type output, verbose: str
    @type padmetRef: <Padmet>
    @type sbml_lvl, sbml_version: int
    @return: check return of writeSBMLToFile
    @rtype: int
    """
    #check if all rxn id are in padmetRef.
    all_rxn = set([k for k,v in padmetRef.dicOfNode.iteritems() if v.type == "reaction"])
    rxn_not_in_ref = reactions.difference(all_rxn)
    if len(rxn_not_in_ref) == len(reactions):
        raise KeyError("None of the reactions is in padmetRef")
    else:
        for rxn_id in rxn_not_in_ref:
            if verbose: print("%s not in padmetRef" % rxn_id)
            reactions.remove(rxn_id)
    padmet = PadmetSpec()
    [padmet.copyNode(padmetRef, rxn_id) for rxn_id in reactions]
    padmet_to_sbml(padmet, output, verbose = verbose)


def compound_to_sbml(species_compart, output, verbose = False):
    """
    convert a list of compounds to sbml format
    if compart_name is not None, then the compounds id will by: M_originalID_compart_name
    if verbose and specified padmetRef and/or padmetSpec: will check if compounds are in one of the padmet files
    Ids are encoded for sbml using functions sbmlPlugin.convert_to_coded_id
    @param compounds_file: the pathname to the file containing the compounds ids and the compart, line = cpd-id\tcompart.
    @param output: the pathname to the sbml file to create
    @param padmetRef_file: the pathname to the file padmet of reference
    @param padmetRef_file: the pathname to the file padmet of a species
    @param compart_name: the default compart to concatenate
    @param sbml_version: the sbml version
    @param verbose: print informations
    @type compounds_file, output, padmetRef_file, padmetSpec_file, verbose: str
    @type sbml_lvl, sbml_version: int
    @return: check return of writeSBMLToFile
    @rtype: int
    """
    document = libsbml.SBMLDocument(2, 1)
    model = document.createModel()

    if verbose: print("%s species" %len(species_compart))
    for data in species_compart:
        species_id = data[0]
        if len(data) == 1:
            compart = "c"
            sId_encoded = sp.convert_to_coded_id(species_id,"M",compart)
        else:
            compart = data[1]
            if compart == "C-BOUNDARY":
                compart = "e"
                sId_encoded = sp.convert_to_coded_id(species_id,"M","e_boundary")
            else:
                sId_encoded = sp.convert_to_coded_id(species_id,"M",compart)
        s = model.createSpecies()
        check(s, 'create species')
        check(s.setId(sId_encoded), 'set species id')
        check(s.setName(species_id), 'set species name')
        check(s.setCompartment(sp.convert_to_coded_id(compart)), 'set species compartment')
    
    libsbml.writeSBMLToFile(document, output)

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
            raise SystemExit(err_msg)
    else:
        return


if __name__ == "__main__":
    main()