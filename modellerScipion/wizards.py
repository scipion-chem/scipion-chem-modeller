# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from .protocols import modellerMutateResidue
from pwem.wizards import GetStructureChainsWizard, pwobj
import pwem.convert as emconv
import requests, os

from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog

class selectChainWizard(GetStructureChainsWizard):
  _targets = [(modellerMutateResidue, ['mutChain'])]

  @classmethod
  def getModelsChainsStep(cls, protocol):
    """ Returns (1) list with the information
       {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
       (2) list with residues, position and chain (modelsFirstResidue)"""
    structureHandler = emconv.AtomicStructHandler()
    fileName = ""
    if hasattr(protocol, 'pdbId'):
      if protocol.pdbId.get() is not None:
        pdbID = protocol.pdbId.get()
        url = "https://www.rcsb.org/structure/"
        URL = url + ("%s" % pdbID)
        try:
          response = requests.get(URL)
        except:
          raise Exception("Cannot connect to PDB server")
        if (response.status_code >= 400) and (response.status_code < 500):
          raise Exception("%s is a wrong PDB ID" % pdbID)
        fileName = structureHandler.readFromPDBDatabase(
          os.path.basename(pdbID), dir="/tmp/")
      else:
        fileName = protocol.pdbFile.get()
    else:
      if protocol.inputAtomStruct.get() is not None:
        fileName = os.path.abspath(protocol.inputAtomStruct.get(
        ).getFileName())

    structureHandler.read(fileName)
    structureHandler.getStructure()
    # listOfChains, listOfResidues = structureHandler.getModelsChains()
    return structureHandler.getModelsChains()

  def show(self, form, *params):
    protocol = form.protocol
    try:
      listOfChains, listOfResidues = self.getModelsChainsStep(protocol)
    except Exception as e:
      print("ERROR: ", e)
      return

    self.editionListOfChains(listOfChains)
    finalChainList = []
    for i in self.chainList:
      finalChainList.append(pwobj.String(i))
    provider = ListTreeProviderString(finalChainList)
    dlg = dialog.ListDialog(form.root, "Model chains", provider,
                            "Select one of the chains (model, chain, "
                            "number of chain residues)")
    form.setVar('mutChain', dlg.values[0].get())

class selectResidueWizard(selectChainWizard):
  _targets = [(modellerMutateResidue, ['mutPosition'])]

  def editionListOfResidues(self, modelsFirstResidue, model, chain):
    self.residueList = []
    for modelID, chainDic in modelsFirstResidue.items():
      if int(model) == modelID:
        for chainID, seq_number in chainDic.items():
          if chain == chainID:
            for i in seq_number:
              self.residueList.append(
                '{"residue": %d, "%s"}' % (i[0], str(i[1])))

  def getResidues(self, form):
    protocol = form.protocol
    try:
      modelsLength, modelsFirstResidue = self.getModelsChainsStep(protocol)
    except Exception as e:
      print("ERROR: ", e)
      return
    selection = protocol.mutChain.get()
    print(selection)

    model = selection.split(',')[0].split(':')[1].strip()
    chain = selection.split(',')[1].split(':')[1].split('"')[1]
    self.editionListOfResidues(modelsFirstResidue, model, chain)
    finalResiduesList = []
    for i in self.residueList:
      finalResiduesList.append(pwobj.String(i))
    return finalResiduesList

  def show(self, form, *params):
    finalResiduesList = self.getResidues(form)
    provider = ListTreeProviderString(finalResiduesList)
    dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                            "Select one residue (residue number, "
                            "residue name)")
    form.setVar('mutPosition', dlg.values[0].get())