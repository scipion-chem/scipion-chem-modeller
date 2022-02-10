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

from .protocols import ModellerMutateResidue
from pwem.wizards import GetStructureChainsWizard, pwobj, EmWizard
import pwem.convert as emconv
import requests, os

from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog
from .constants import AA_LIST

from pwchem.wizards import SelectChainWizard

class SelectChainWizardModeller(SelectChainWizard):
  _targets = [(ModellerMutateResidue, ['mutChain'])]

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

class SelectResidueWizardModeller(SelectChainWizard):
  _targets = [(ModellerMutateResidue, ['mutPosition'])]

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

class AddMutationWizard(EmWizard):
  _targets = [(ModellerMutateResidue, ['addMutation'])]

  def show(self, form, *params):
    protocol = form.protocol
    chain, pos, resId = protocol.mutChain.get(), protocol.mutPosition.get(), protocol.mutResidue.get()
    res = AA_LIST[resId]
    form.setVar('toMutateList', protocol.toMutateList.get() +
                '{} | {} | {}\n'.format(chain, pos, res))
    
class ClearMutationList(EmWizard):
  _targets = [(ModellerMutateResidue, ['clearLabel'])]

  def show(self, form, *params):
    form.setVar('toMutateList', '')
