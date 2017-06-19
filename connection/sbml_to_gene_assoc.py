# -*- coding: utf-8 -*-
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
usage:
    sbml_to_gene_assoc.py --sbml=FILE [--output=FILE]

option:
    -h --help    Show help.
    --sbml=FILE    pathname to the sbml file
    --output=FILE    pathname to the output containing the gene association, line = reaction id,genes associated, sep = \t
"""
import re
from libsbml import *
import docopt

def main():
    args = docopt.docopt(__doc__)
    sbml_file = args["--sbml"]
    #using libSbml to read sbml_file
    print ("loading sbml file: %s" %sbml_file)
    reader = SBMLReader()
    document = reader.readSBML(sbml_file)
    for i in range(document.getNumErrors()):
        print (document.getError(i).getMessage())

    model = document.getModel()
    listOfReactions = model.getListOfReactions()
    nbReactions = len(listOfReactions)
    print("nb reactions: %s" %nbReactions)
    print("Reaction id\tGene association")
    for reaction in listOfReactions:
        notes = parseNotes(reaction)
        try:
            gene_assoc = parseGeneAssoc(notes["GENE_ASSOCIATION"][0])
            print (reaction.getId()+"\t"+";".join(gene_assoc))
        except KeyError:
            print ("No gene association for %s" %reaction.getId())
            
            






def parseNotes(element):
    """
    From an SBML element (ex: species or reaction) will return all the section
    note in a dictionnary.
    ex:
    <notes>
        <html:body>
            <html:p>BIOCYC: |Alkylphosphonates|</html:p>
            <html:p>CHEBI: 60983</html:p>
        </html:body>
     </notes>
    output: {'BIOCYC': ['Alkylphosphonates'],'CHEBI':['60983']}
    value is a list in case diff lines for the same type of info

    @param element: an element from libsbml
    @type element: libsbml.element
    @return: the dictionnary of note
    @rtype: dict
    """
    notes = element.getNotesString()
    notesList = notes.splitlines()
    notesDict = {}
    for line in notesList:
        try:
            #line = <html:p>BIOCYC: |Alkylphosphonates|</html:p>
            start = line.index(">")+1
            end = line.index("<",start)
            line = line[start:end]
            #line = BIOCYC: |Alkylphosphonates|
            line = line.split(":")
            #line = [BIOCYC,|Alkylphosphonates|]
            line[0] = re.sub(" ","_",line[0])
            line[1] = re.sub("\s|\|","",line[1])
            if len(line[1]) != 0:
                line[1] = line[1].split(",")
                notesDict[line[0]] = line[1]
        except ValueError:
            continue
    return notesDict

def parseGeneAssoc(GeneAssocStr):
    """
    Given a grammar of 'and', 'or' and '(' ')'. Extracts genes ids to a list.
    (geneX and geneY) or geneW' => [geneX,geneY,geneW]
    @param GeneAssocStr: the string containing genes ids
    @type GeneAssocStr: str
    @return: the list of unique ids
    @rtype: list
    """
    #sub '(',')',' ' by ''   sub "and" by "or"
    resultat = re.sub("\(|\)|\s","",GeneAssocStr).replace("and","or")
    #create a set by spliting 'or' then convert to list, set for unique genes
    if len(resultat) != 0:
        resultat = list(set(resultat.split("or")))
    else:
        resultat = []
    return resultat

if __name__ == "__main__":
    main()