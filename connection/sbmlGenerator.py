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
    sbmlGenerator.py --padmet=FILE --output=FILE [--model_id=STR] [--model_name=STR] [--obj_fct=STR] [--sbml_lvl=x] [--mnx_chem_prop=FILE] [--mnx_chem_xref=FILE] [-v]
    sbmlGenerator.py --compound=FILE --output=FILE [--padmetRef=FILE] [-v]
    sbmlGenerator.py --reaction=FILE --output=FILE --padmetRef=FILE [-v]

option:
    -h --help    Show help.
    --padmet=FILE    pathname of the padmet file to convert into sbml
    --output=FILE    pathanme of the sbml file to generate.
    --obj_fct=STR    id of the reaction objective.
    --sbml_lvl=x    sbml level 2 is sufficient for FBA.
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
        obj_fct = args["--obj_fct"]
        mnx_chem_xref = args["--mnx_chem_xref"]
        mnx_chem_prop = args["--mnx_chem_prop"]
        if sbml_lvl == 2:
            padmet_to_sbml2(padmet, output, obj_fct, verbose)
        elif sbml_lvl == 3:
            padmet_to_sbml3(padmet, output, model_id, model_name, obj_fct, mnx_chem_xref, mnx_chem_prop, verbose)
    elif args["--reaction"]:
        with open(args["--reaction"], 'r') as f:
            reactions = set(f.read().splitlines())
        padmetRef = PadmetSpec(args["--padmetRef"])
        reaction_to_sbml(reactions, output, padmetRef, verbose)
    elif args["--compound"]:
        with open(args["--compound"], 'r') as f:
            species_compart = [line.split("\t") for line in f.read().splitlines()]
        compound_to_sbml(species_compart, output, verbose)



def padmet_to_sbml3(padmet, output, model_id = None, model_name = None, obj_fct = None, mnx_chem_xref = None, mnx_chem_prop = None, verbose = False):
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
                dict_mnx_chem_xref[k] = v
                    
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
                mnx_id = dict_mnx_chem_xref[species_id]
                species_prop = dict(dict_mnx_chem_prop[mnx_id])
                species_formula = str(species_prop["formula"])
                for k,v in species_prop.items():
                    if k in ["formula", "source", "description"] or v in ["NA",""]:
                        species_prop.pop(k)
            except (IndexError, KeyError) as e:
                #print(species_id)
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
        #get all categories
        categories = set([padmet.dicOfNode[rlt.id_out].misc["CATEGORY"][0] for rlt in padmet.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"])
        #get consumed and produced compounds
        consumed = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc["COMPARTMENT"][0]) 
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "consumes" and rlt.misc["COMPARTMENT"][0] != "C-BOUNDARY")
                
        produced = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc["COMPARTMENT"][0]) 
        for rlt in padmet.dicOfRelationIn.get(rxn_id, None) if rlt.type == "produces" and rlt.misc["COMPARTMENT"][0] != "C-BOUNDARY")
    
        model.add_reaction(sbml_rxn)
        for species, stoichio, compart in consumed:
            metabolitedict = {}
            try:
                stoichio = float(stoichio)*-1
            except ValueError:
                stoichio = float(-1.0)
            metabolitedict[sp.convert_to_coded_id(species,None,compart)] = stoichio
            sbml_rxn.add_metabolites(metabolitedict)
        for species, stoichio, compart in produced:
            metabolitedict = {}
            try:
                stoichio = float(stoichio)
            except ValueError:
                stoichio = float(1.0)
            metabolitedict[sp.convert_to_coded_id(species,None,compart)] = stoichio
            sbml_rxn.add_metabolites(metabolitedict)
        model.repair()
    
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
        #if lvl 2: 
        #sbml_rxn.notes  = {"CATEGORY": " and ".join(categories)}
        #sbml_rxn.notes.update(dict({"CATEGORY": "test"}))
    #use_fbc_package = not (sbml_lvl == 2)
    write_cobra_model_to_sbml_file(model, output, sbml_level=2, use_fbc_package=True)


def padmet_to_sbml2(padmet, output, obj_fct = None, verbose = False):
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
    #create an empty sbml model
    document = libsbml.SBMLDocument(2, 1)
    model = document.createModel()
    model_id = os.path.splitext(os.path.basename(output))[0]
    model.id = model_id

    BOUNDARY_ID = 'C-BOUNDARY'
    default_lower_bound = -1000
    default_upper_bound = 1000
    math_ast = libsbml.parseL3Formula('FLUX_VALUE')
    check(math_ast, 'create AST for rate expression')

    # Create a unit definition
    mmol_per_gDW_per_hr = model.createUnitDefinition()
    check(mmol_per_gDW_per_hr, 'create unit definition')
    check(mmol_per_gDW_per_hr.setId('mmol_per_gDW_per_hr'), 'set unit definition id')
    
    unit = mmol_per_gDW_per_hr.createUnit()
    check(unit, 'create mole unit')
    check(unit.setKind(libsbml.UNIT_KIND_MOLE), 'set unit kind')
    check(unit.setScale(-3), 'set unit scale')
    check(unit.setMultiplier(1), 'set unit multiplier')
    check(unit.setOffset(0), 'set unit offset')
    
    unit = mmol_per_gDW_per_hr.createUnit()
    check(unit, 'create gram unit')
    check(unit.setKind(libsbml.UNIT_KIND_GRAM), 'set unit kind')
    check(unit.setExponent(-1), 'set unit exponent')
    check(unit.setMultiplier(1), 'set unit multiplier')
    check(unit.setOffset(0), 'set unit offset')
    
    unit = mmol_per_gDW_per_hr.createUnit()
    check(unit, 'create second unit')
    check(unit.setKind(libsbml.UNIT_KIND_SECOND), 'set unit kind')
    check(unit.setExponent(-1), 'set unit exponent')
    check(unit.setMultiplier(0.00027777), 'set unit multiplier')
    check(unit.setOffset(0), 'set unit offset')


    #generator of tuple: (x,y) x=species id,y=value of compart, if not defined=""
    species = [(rlt.id_out, rlt.misc.get("COMPARTMENT",[None])[0]) for rlt in padmet.getAllRelation() 
    if rlt.type in ["consumes","produces"]]
    if verbose: print("%s species" %len(species))
    #compart_dict: k = id_encoded, v = original id
    compart_dict = {}    
    #species_dict: k = species_id_encoded, v = dict: k' = {species_id, compart, name}, v' = value or None 
    species_dict = {}
    for species_id, compart in species:
        #encode id for sbml
        species_id_encoded = sp.convert_to_coded_id(species_id, "M", compart)
        
        #encode compart id for sbml
        #try to get the common_name, if non value return None
        name = padmet.dicOfNode[species_id].misc.get("COMMON-NAME",[species_id])[0]
        #update dicts
        species_dict[species_id_encoded] = {"species_id":species_id, "compart":compart, "name":name}
        
    
    for k, v in species_dict.iteritems():
        compart = v["compart"]
        name = v["name"]
        
        s = model.createSpecies()
        check(s, 'create species')
        check(s.setId(k), 'set species id %s' %k)
        check(s.setBoundaryCondition(False), 'set boundaryCondition to False')
        #check(s.setMetaId(metaId), 'set species MetaId %s' %metaId)
        if name is not None:
            check(s.setName(name), 'set species Name %s' %name)
        else:
            check(s.setName(name), 'set species Name %s' %species_id)

        if compart is not None:
            compart_encoded = sp.convert_to_coded_id(compart)
            compart_dict[compart_encoded] = compart
            check(s.setCompartment(compart_encoded), 'set species compartment %s' %compart_encoded)
            if compart == BOUNDARY_ID:
                check(s.setBoundaryCondition(True), 'set boundaryCondition to True')

    for k, v in compart_dict.iteritems():
        compart = model.createCompartment()
        check(compart,'create compartment')
        check(compart.setId(k),'set compartment id %s' %k)
        if v == "c":
            check(compart.setName("cytosol"),'set compartment name cytosol')
        elif v == "e":
            check(compart.setName("extracellular"),'set compartment name extracellular')
        elif v == "p":
            check(compart.setName("periplasm"),'set compartment name periplasm')
        elif v != k:
            check(compart.setName(v),'set compartment id %s' %v)


    if obj_fct is not None:
        obj_fct_encoded = sp.convert_to_coded_id(obj_fct)
        if verbose: print("the objectif reaction is: %s initialy: %s" %(obj_fct_encoded, obj_fct))
        
    reactions = [node for node in padmet.dicOfNode.itervalues() if node.type == "reaction"]
    nb_reactions = str(len(reactions))    
    # Create reactions
    if verbose: print("%s reactions" %nb_reactions)

    for rNode in reactions:
        rId = rNode.id
        rId_encoded = sp.convert_to_coded_id(rId,"R")
        rName = rNode.misc.get("COMMON-NAME",[rId])[0]

        #generator of tuple (reactant_id,stoichiometry,compart)
        try:
            consumed = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc.get("COMPARTMENT",[None])[0]) 
            for rlt in padmet.dicOfRelationIn.get(rId, None) if rlt.type == "consumes")
        except TypeError:
            print rId
            exit()
            
        #generator of tuple (product_id,stoichiometry,compart)        
        produced = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc.get("COMPARTMENT",[None])[0]) 
        for rlt in padmet.dicOfRelationIn.get(rId, None) if rlt.type == "produces")
        direction = rNode.misc["DIRECTION"][0]
        if direction == "LEFT-TO-RIGHT":
            reversible = False
        #include if direction = unknown
        else:
            reversible = True
        
        reaction = model.createReaction()
        check(reaction, 'create reaction')
        check(reaction.setId(rId_encoded), 'set reaction id %s' %rId_encoded)
        if rName is not None:
            check(reaction.setName(rName), 'set reaction name %s' %rName)
        check(reaction.setReversible(reversible), 'set reaction reversibility flag %s' %reversible)
        
        kinetic_law = reaction.createKineticLaw()
        check(kinetic_law, 'create kinetic law')
        check(kinetic_law.setMath(math_ast), 'set math on kinetic law')
        #add parameter flux_value
        flux_value_k = kinetic_law.createParameter()
        check(flux_value_k, 'create parameter flux_value_k')
        check(flux_value_k.setId('FLUX_VALUE'), 'set parameter flux_value_k id')
        check(flux_value_k.setValue(0), 'set parameter flux_value_k value')
        check(flux_value_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter flux_value_k units')
        #add parameter upper/lower_bound, lower value depend on reversibility
        upper_bound_k = kinetic_law.createParameter()
        check(upper_bound_k, 'create parameter upper_bound_k')
        check(upper_bound_k.setId('UPPER_BOUND'), 'set parameter upper_bound_k')
        check(upper_bound_k.setValue(default_upper_bound),'set parameter upper_bounp_k value')
        check(upper_bound_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter uppper_bound_k units')

        if reversible:
            lower_bound_k = kinetic_law.createParameter()
            check(lower_bound_k, 'create parameter lower_bound_k')
            check(lower_bound_k.setId('LOWER_BOUND'), 'set parameter lower_bound_k id')
            check(lower_bound_k.setValue(default_lower_bound), 'set parameter lower_bound_k value')
            check(lower_bound_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter lower_bound_k units')
        else:
            lower_bound_k = kinetic_law.createParameter()
            check(lower_bound_k, 'create parameter lower_bound_k')
            check(lower_bound_k.setId('LOWER_BOUND'), 'set parameter lower_bound_k id')
            check(lower_bound_k.setValue(0), 'set parameter lower_bound_k value')
            check(lower_bound_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter lower_bound_k units')
        #objective_coeeficient
        if rId == obj_fct:
            obj_fct_k = kinetic_law.createParameter()
            check(obj_fct_k, 'create parameter obj_fct_k')
            check(obj_fct_k.setId('OBJECTIVE_COEFFICIENT'), 'set parameter obj_fct_k id')
            check(obj_fct_k.setValue(1), 'set parameter obj_fct_k value')
        else:
            obj_fct_k = kinetic_law.createParameter()
            check(obj_fct_k, 'create parameter obj_fct_k')
            check(obj_fct_k.setId('OBJECTIVE_COEFFICIENT'), 'set parameter obj_fct_k id')
            check(obj_fct_k.setValue(0), 'set parameter obj_fct_k value')

        for cId, stoich, compart in consumed:
            cId_encoded = sp.convert_to_coded_id(cId,"M",compart)
            try:
                stoich = float(stoich)
            #for case stoich = n
            except ValueError:
                stoich = float(1)
            species_ref = reaction.createReactant()
            check(species_ref, 'create reactant')
            check(species_ref.setSpecies(cId_encoded), 'assign reactant species %s' %cId_encoded)
            check(species_ref.setStoichiometry(stoich), 'set stoichiometry %s' %stoich)

        for pId, stoich, compart in produced:
            pId_encoded = sp.convert_to_coded_id(pId,"M",compart)
            try:
                stoich = float(stoich)
            except ValueError:
                stoich = float(1)
            species_ref = reaction.createProduct()
            check(species_ref, 'create product')
            check(species_ref.setSpecies(pId_encoded), 'assign product species %s' %pId_encoded)
            check(species_ref.setStoichiometry(stoich), 'set stoichiometry %s' %stoich)

        try:
            linked_genes = set([rlt.id_out for rlt in padmet.dicOfRelationIn.get(rId, [])
            if rlt.type == "is_linked_to"])
            if len(linked_genes) != 0:
                if verbose: print rId, "is linked to", len(linked_genes), "genes."
                notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
                linked_genes = " or ".join(linked_genes)
                notes += "<p>"+"GENE_ASSOCIATION:" + linked_genes + "</p>"
                notes += "</body>"
                check(reaction.setNotes(notes), 'set notes %s' %notes)
        except IndexError:
            pass

    if verbose: print("Done, creating sbml file: %s" %output)
    libsbml.writeSBMLToFile(document, output)


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
        else:
            compart = data[1]
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