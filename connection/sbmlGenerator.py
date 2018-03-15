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
    sbmlGenerator.py --padmet=FILE --output=FILE --sbml_lvl=STR [--model_id=STR] [--obj_fct=STR] [--mnx_chem_prop=FILE] [--mnx_chem_xref=FILE] [-v]
    sbmlGenerator.py --padmet=FILE --output=FILE [--init_source=STR] [-v]
    sbmlGenerator.py --compound=FILE --output=FILE [--padmetRef=FILE] [-v]
    sbmlGenerator.py --reaction=FILE --output=FILE --padmetRef=FILE [-v]

option:
    -h --help    Show help.
    --padmet=FILE    pathname of the padmet file to convert into sbml
    --output=FILE    pathanme of the sbml file to generate.
    --obj_fct=STR    id of the reaction objective.
    --sbml_lvl=STR    sbml level 2 is sufficient for FBA.
    -v   print info.
"""
from padmet.padmetSpec import PadmetSpec
import padmet.sbmlPlugin as sp
import os
import docopt
import libsbml
import re
#default variables
global def_max_lower_bound, def_max_upper_bound, BOUNDARY_ID
def_max_upper_bound = 1000
def_max_lower_bound = -1000
BOUNDARY_ID = 'C-BOUNDARY'

def main():
    args = docopt.docopt(__doc__)
    output = args["--output"]
    verbose = args["-v"]
    if args["--padmet"]:
        padmet = PadmetSpec(args["--padmet"])
        if args["--init_source"]:
            rxn_to_del = set()
            init_source = args["--init_source"].lower()
            rec_rlts = [rlt for rlt in padmet.getAllRelation() if rlt.type == "has_reconstructionData"]
            [rxn_to_del.add(rlt.id_in) for rlt in rec_rlts if padmet.dicOfNode[rlt.id_out].misc.get("SOURCE",[""])[0].lower() == init_source]            
            reactions_to_add = set([node.id for node in padmet.dicOfNode.values() if node.type == "reaction"]).difference(rxn_to_del)
            reaction_to_sbml(reactions_to_add, output, padmet, verbose)
            return
        if args["--model_id"]:
            model_id = args["--model_id"]
        else:
            model_id = os.path.splitext(os.path.basename(output))[0]
        if args["--sbml_lvl"]:
            sbml_lvl = int(args["--sbml_lvl"])
        else:
            sbml_lvl = 3
        obj_fct = args["--obj_fct"]
        mnx_chem_xref = args["--mnx_chem_xref"]
        mnx_chem_prop = args["--mnx_chem_prop"]
        padmet_to_sbml(padmet, output, model_id, obj_fct, sbml_lvl, mnx_chem_prop, mnx_chem_xref, verbose)
    elif args["--reaction"]:
        with open(args["--reaction"], 'r') as f:
            reactions = set(f.read().splitlines())
        padmetRef = PadmetSpec(args["--padmetRef"])
        reaction_to_sbml(reactions, output, padmetRef, verbose)
    elif args["--compound"]:
        with open(args["--compound"], 'r') as f:
            species_compart = [line.split("\t") for line in f.read().splitlines()]
        compound_to_sbml(species_compart, output, verbose)

def padmet_to_sbml(padmet, output, model_id = None, obj_fct = None, sbml_lvl = 3, mnx_chem_prop = None, mnx_chem_xref = None, verbose = False):
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
    global all_ga
    all_ga = []
    #create an empty sbml model
    with_mnx = False
    if mnx_chem_prop and mnx_chem_xref:
        with_mnx = True
        dict_mnx_chem_xref = parse_mnx_chem_xref(mnx_chem_xref)
        dict_mnx_chem_prop = parse_mnx_chem_prop(mnx_chem_prop)
    if sbml_lvl == 2:
        sbmlns = libsbml.SBMLNamespaces(2,1)
        document = libsbml.SBMLDocument(sbmlns)
        model = document.createModel()
        association = None
    elif sbml_lvl == 3:
        sbmlns = libsbml.SBMLNamespaces(3,1,"fbc",1)
        document = libsbml.SBMLDocument(sbmlns)
        document.setPackageRequired("fbc", False)
        model = document.createModel()
        mplugin = model.getPlugin("fbc")
        association = ['<annotation>', 
        '<listOfGeneAssociations xmlns="http://www.sbml.org/sbml/level3/version1/fbc/version1">']
    if not model_id: model_id = os.path.splitext(os.path.basename(output))[0]
    model.setId(model_id)
    
    """
    #set units
    #model.setSubstanceUnits("substance")
    model.setTimeUnits("second")
    #model.setVolumeUnits("volume")
    model.setAreaUnits("area")
    model.setLengthUnits("metre")
    #model.setExtentUnits("substance")
    """
    check(model,                              'create model')
    check(model.setTimeUnits("second"),       'set model-wide time units')
    check(model.setExtentUnits("mole"),       'set model units of extent')
    check(model.setSubstanceUnits('mole'),    'set model substance units')
    math_ast = libsbml.parseL3Formula('FLUX_VALUE')
    check(math_ast, 'create AST for rate expression')

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
        original_id = v["species_id"]
        s = model.createSpecies()
        check(s, 'create species')
        check(s.setId(k), 'set species id %s' %k)
        check(s.setBoundaryCondition(False), 'set boundaryCondition to False')
        check(s.setHasOnlySubstanceUnits(False), 'set setHasOnlySubstanceUnits to False')
        check(s.setConstant(False), 'set setConstant to False')
        check(s.setInitialAmount(0.0), 'set initAmount')
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
        if with_mnx:
            try:
                mnx_id = dict_mnx_chem_xref[original_id]
                species_prop = dict(dict_mnx_chem_prop[mnx_id])
            except (IndexError, KeyError) as e:
                #print(species_id)
                species_prop = None
            if species_prop:
                try:
                    charge = int(species_prop["charge"])
                except ValueError:
                    charge = 0
                formula = species_prop["formula"]
                #remove () from forumla
                #formula = re.sub("\(|\)|\.","",formula).upper()
                if re.findall("\(|\)|\.",formula): formula = None
                if sbml_lvl == 3:
                    splugin = s.getPlugin("fbc")
                    check(splugin.setCharge(charge), 'set charge')
                    if formula:
                        check(splugin.setChemicalFormula(formula), 'set Formula')
                    for prop, prop_v in species_prop.items():
                        if prop in ["charge", "formula", "source", "description"] or prop_v in ["NA",""]:
                            species_prop.pop(prop)
                notes = create_note(species_prop)
                check(s.setNotes(notes), 'set Notes')

    for k, v in compart_dict.iteritems():
        compart = model.createCompartment()
        check(compart,'create compartment')
        check(compart.setId(k),'set compartment id %s' %k)
        check(compart.setSize(1),'set size for compartment id %s' %k)
        check(compart.setConstant(True),'set constant for compartment id %s' %k)

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
        if verbose: print("the objectif reaction is: %s" %(obj_fct_encoded))
        
    reactions = [node for node in padmet.dicOfNode.itervalues() if node.type == "reaction"]
    nb_reactions = str(len(reactions))    
    # Create reactions
    if verbose: print("%s reactions" %nb_reactions)
    for rNode in reactions:
        rId = rNode.id
        rId_encoded = sp.convert_to_coded_id(rId,"R")
        rName = rNode.misc.get("COMMON-NAME",[rId])[0]
        reaction = model.createReaction()
        check(reaction, 'create reaction')
        check(reaction.setId(rId_encoded), 'set reaction id %s' %rId_encoded)
        if rName is not None:
            check(reaction.setName(rName), 'set reaction name %s' %rName)
        check(reaction.setFast(False), 'set fast')

        #generator of tuple (reactant_id,stoichiometry,compart)
        consumed = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc.get("COMPARTMENT",[None])[0]) 
        for rlt in padmet.dicOfRelationIn.get(rId, None) if rlt.type == "consumes")
        #generator of tuple (product_id,stoichiometry,compart)        
        produced = ((rlt.id_out, rlt.misc["STOICHIOMETRY"][0], rlt.misc.get("COMPARTMENT",[None])[0]) 
        for rlt in padmet.dicOfRelationIn.get(rId, None) if rlt.type == "produces")
        #set reversibility
        direction = rNode.misc["DIRECTION"][0]
        if direction == "LEFT-TO-RIGHT":
            reversible = False
        else:
            reversible = True
        check(reaction.setReversible(reversible), 'set reaction reversibility flag %s' %reversible)
        if sbml_lvl == 3:
            bound= mplugin.createFluxBound()
            bound.setReaction(rId_encoded)
            bound.setOperation("lessEqual")
            bound.setValue(def_max_upper_bound)

            bound= mplugin.createFluxBound()
            bound.setReaction(rId_encoded)
            bound.setOperation("greaterEqual")
            if reversible:
                bound.setValue(def_max_lower_bound)
            else:
                bound.setValue(0)
            if rId == obj_fct:
                 objective = mplugin.createObjective()
                 objective.setId("obj1")
                 objective.setType("maximize")
                 mplugin.setActiveObjectiveId("obj1")
                 fluxObjective = objective.createFluxObjective()
                 fluxObjective.setReaction(rId_encoded)
                 fluxObjective.setCoefficient(1)
        elif sbml_lvl == 2:
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
            check(upper_bound_k.setValue(def_max_upper_bound),'set parameter upper_bounp_k value')
            check(upper_bound_k.setUnits('mmol_per_gDW_per_hr'), 'set parameter uppper_bound_k units')
    
            if reversible:
                lower_bound_k = kinetic_law.createParameter()
                check(lower_bound_k, 'create parameter lower_bound_k')
                check(lower_bound_k.setId('LOWER_BOUND'), 'set parameter lower_bound_k id')
                check(lower_bound_k.setValue(def_max_lower_bound), 'set parameter lower_bound_k value')
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
            check(species_ref.setStoichiometry(stoich), 'set stoichiometry %s' %stoich)
            check(species_ref.setConstant(False), 'set constant %s' %False)

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
            check(species_ref.setConstant(False), 'set constant %s' %False)

        #set notes
        notes_dict = {}
        linked_genes = set([rlt.id_out for rlt in padmet.dicOfRelationIn.get(rId, [])
        if rlt.type == "is_linked_to"])
        if linked_genes:
            notes_dict["GENE ASSOCIATION"] = " or ".join(linked_genes)
            if association:
                add_ga(list(linked_genes), rId_encoded)
        categories = set([padmet.dicOfNode[rlt.id_out].misc["CATEGORY"][0] for rlt in padmet.dicOfRelationIn.get(rId,[]) if rlt.type == "has_reconstructionData"])
        if categories:
            notes_dict["CATEGORIES"] = " and ".join(categories)
        pathways = set([rlt.id_out for rlt in padmet.dicOfRelationIn.get(rId, [])
        if rlt.type == "is_in_pathway"])
        if len(pathways) != 0:
            notes_dict["SUBSYSTEM"] = " , ".join(pathways)
        if notes_dict.keys():
            notes = create_note(notes_dict)            
            check(reaction.setNotes(notes), 'set notes %s' %notes)
    if all_ga:
        for ga in all_ga:
            association.extend(ga)
        association.extend(['</listOfGeneAssociations>', '</annotation>'])
        association = " ".join(association)
        model.setAnnotation(association)
    if verbose: print("Done, creating sbml file: %s" %output)
    libsbml.writeSBMLToFile(document, output)

def parse_mnx_chem_xref(mnx_chem_xref):
    #k = xref id, v = mnx_id
    dict_mnx_chem_xref = {}
    with open(mnx_chem_xref, 'r') as f:
        for k,v in [line.split("\t")[:2] for line in f.read().splitlines() if not line.startswith("#") and not line.startswith("MNX")]:
            if ":" in k: k = k.split(":")[1]
            dict_mnx_chem_xref[k] = v
    return dict_mnx_chem_xref

def parse_mnx_chem_prop(mnx_chem_prop):
    #k=mnx_id, v = dict of data
    dict_mnx_chem_prop = {}
    with open(mnx_chem_prop, 'r') as f:
        for line in (line.split("\t") for line in f.read().splitlines() if not line.startswith("#")):
            dict_mnx_chem_prop[line[0]] = {"description":line[1],"formula":line[2],"charge":line[3],"mass":line[4],"inchi":line[5],"smiles":line[6],"source":line[7],"inchikey":line[8]}
    return dict_mnx_chem_prop

def create_note(dict_data):
    notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">"
    for k,v in dict_data.items():
        notes += "<p>"+k+": "+v+"</p>"
    notes += "</body>"
    return notes
    
def add_ga(linked_genes, rId_encoded):
    global all_ga
    ga_count = len(all_ga) + 1
    ga = [' <geneAssociation id="ga_'+str(ga_count)+'" reaction="'+rId_encoded+'">']
    if len(linked_genes) == 1:
        ga.append('<gene reference="'+linked_genes[0]+'"/>')
    else:
        ga.append('<or>')
        for gene_id in linked_genes:
            ga.append('<gene reference="'+gene_id+'"/>')
        ga.append('</or>')
    ga.append('</geneAssociation>')
    all_ga.append(ga)        
    

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
    padmet_to_sbml(padmet, output, sbml_lvl=2, verbose = verbose)


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
            raise TypeError(err_msg)
    else:
        return


if __name__ == "__main__":
    main()