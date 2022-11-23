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

import json, os, random, requests

from .protocols import ModellerMutateResidue, ProtModellerComparativeModelling
from pwem.wizards import SelectChainWizard, SelectResidueWizard, EmWizard, VariableWizard
import pwem.objects as emobj
import pwem.convert as emconv

from pwchem.utils import downloadPDB
from pwchem.wizards import SelectChainWizardQT, SelectResidueWizardQT, SelectMultiChainWizard

from .constants import AA_LIST

SelectChainWizardQT().addTarget(protocol=ModellerMutateResidue,
                              targets=['mutChain'],
                              inputs=['inputAtomStruct'],
                              outputs=['mutChain'])

SelectResidueWizardQT().addTarget(protocol=ModellerMutateResidue,
                                targets=['mutPosition'],
                                inputs=['inputAtomStruct', 'mutChain'],
                                outputs=['mutPosition'])

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


class AddStructSequenceWizard(SelectResidueWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getModelsChainsStep(self, protocol, inputObj):
      """ Returns (1) list with the information
         {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
         (2) list with residues, position and chain (modelsFirstResidue)"""
      structureHandler = emconv.AtomicStructHandler()
      if type(inputObj) == str:
        if os.path.exists(inputObj):
          fileName = inputObj
        else:
          fileName = downloadPDB(inputObj)

      elif str(type(inputObj).__name__) == 'SchrodingerAtomStruct':
        fileName = os.path.abspath(inputObj.convert2PDB())
      else:
        fileName = os.path.abspath(inputObj.getFileName())

      structureHandler.read(fileName)
      structureHandler.getStructure()
      return structureHandler.getModelsChains()

    def getResidues(self, form, inputObj, modelChain):
      protocol = form.protocol

      if type(inputObj) == str:
          # Select the residues if the input structure parameter is a str (PDB id or PDB file)
          structureHandler = emconv.AtomicStructHandler()
          if os.path.exists(inputObj):
              fileName = inputObj
          else:
              fileName = downloadPDB(inputObj)

          structureHandler.read(fileName)
          structureHandler.getStructure()
          modelsLength, modelsFirstResidue = structureHandler.getModelsChains()

          model, chain = modelChain.split('-')
          residueList = self.editionListOfResidues(modelsFirstResidue, model, chain)
          finalResiduesList = []
          for i in residueList:
              finalResiduesList.append(emobj.String(i))

      elif issubclass(type(inputObj), emobj.AtomStruct):
          try:
            modelsLength, modelsFirstResidue = self.getModelsChainsStep(protocol, inputObj)
          except Exception as e:
            print("ERROR: ", e)
            return
          model, chain = modelChain.split('-')

          residueList = self.editionListOfResidues(modelsFirstResidue, model, chain)
          finalResiduesList = []
          for i in residueList:
            finalResiduesList.append(emobj.String(i))
      return finalResiduesList

    def show(self, form, *params):
        protocol = form.protocol
        pseudoProtId = random.randint(0, 1000)  #we cannot obtain prot Id because it might not be launched yet
        inputParams, outputParam = self.getInputOutput(form)

        # StructName
        inputTemplate = getattr(protocol, inputParams[0]).get()
        pdbFile = ''
        if issubclass(type(inputTemplate), str):
            outName = inputTemplate.lower()
        else:
            pdbFile = inputTemplate.getFileName()
            outName = os.path.splitext(os.path.basename(pdbFile))[0].lower()

        if not protocol.multiChain:
            # Chain
            chainJson = getattr(protocol, inputParams[1]).get()
            modelId, chainId = json.loads(chainJson)['model'], json.loads(chainJson)['chain']

            # Positions
            posJson = getattr(protocol, inputParams[2]).get()
            if posJson:
                posIdxs = json.loads(posJson)['index']
                seq = json.loads(posJson)['residues']
            else:
                posIdxs = 'FIRST-LAST'
                finalResiduesList = self.getResidues(form, inputTemplate, '{}-{}'.format(modelId, chainId))
                idxs = [json.loads(finalResiduesList[0].get())['index'],
                        json.loads(finalResiduesList[-1].get())['index']]
                seq = self.getSequence(finalResiduesList, idxs)

            seqFile = protocol.getProject().getTmpPath('{}_{}_{}_{}.fa'.format(outName, pseudoProtId, chainId, posIdxs))
            with open(seqFile, 'w') as f:
                f.write('>{}\n{}\n'.format(outName, seq))

            prevStr = getattr(protocol, outputParam[0]).get()
            lenPrev = len(prevStr.strip().split('\n')) + 1
            if prevStr.strip() == '':
              lenPrev -= 1
            elif not prevStr.endswith('\n'):
              prevStr += '\n'

            pdbStr = ''
            if pdbFile:
                pdbStr = ', "pdbFile": "{}"'.format(pdbFile)
            jsonStr = '%s) {"pdbName": "%s", "chain": "%s", "index": "%s", "seqFile": "%s"%s}\n' % \
                      (lenPrev, outName, chainId, posIdxs, seqFile, pdbStr)
            form.setVar(outputParam[0], prevStr + jsonStr)
        else:
            # Chain
            chainsList = json.loads(getattr(protocol, inputParams[1]).get())['model-chain']
            for i, chainStr in enumerate(chainsList.split(',')):
                chain = chainStr.split('-')[1]

                # Positions
                finalResiduesList = self.getResidues(form, inputTemplate, chainStr)
                posStr = getattr(protocol, inputParams[2]).get()
                if posStr:
                    posIdxs = posStr.split(',').strip()
                else:
                    posIdxs = [json.loads(finalResiduesList[0].get())['index'],
                               json.loads(finalResiduesList[-1].get())['index']]

                seq = self.getSequence(finalResiduesList, posIdxs)

                posIdxStr = '-'.join(list(map(str, posIdxs)))
                seqFile = protocol.getProject().getTmpPath('{}_{}_{}_{}.fa'.format(outName, pseudoProtId, chain,
                                                                                   posIdxStr))
                with open(seqFile, 'w') as f:
                    f.write('>{}_{}\n{}\n'.format(outName, chain, seq))

            seqFileTemplate = seqFile = protocol.getProject().getTmpPath('{}_{}_*_*.fa'.format(outName, pseudoProtId))

            prevStr = getattr(protocol, outputParam[0]).get()
            lenPrev = len(prevStr.strip().split('\n')) + 1
            if prevStr.strip() == '':
              lenPrev -= 1
            elif not prevStr.endswith('\n'):
              prevStr += '\n'

            pdbStr = ''
            if pdbFile:
              pdbStr = ', "pdbFile": "{}"'.format(pdbFile)
            jsonStr = '%s) {"pdbName": "%s", "chains": "%s", "seqFiles": "%s"%s}\n' % \
                      (lenPrev, outName, chainsList, seqFile, pdbStr)
            form.setVar(outputParam[0], prevStr + jsonStr)

AddStructSequenceWizard().addTarget(protocol=ProtModellerComparativeModelling,
                                    targets=['addTemplate'],
                                    inputs=[{'templateOrigin': ['inputAtomStruct', 'pdbTemplate']},
                                            {'multiChain': ['tempChain', 'tempMChain']},
                                            {'multiChain': ['tempPositions', 'tempMPositions']}],
                                    outputs=['templateList'])

SelectChainWizardQT().addTarget(protocol=ProtModellerComparativeModelling,
                                targets=['tempChain'],
                                inputs=[{'templateOrigin': ['inputAtomStruct', 'pdbTemplate']}],
                                outputs=['tempChain'])

SelectResidueWizardQT().addTarget(protocol=ProtModellerComparativeModelling,
                                  targets=['tempPositions'],
                                  inputs=[{'templateOrigin': ['inputAtomStruct', 'pdbTemplate']},
                                          'tempChain'],
                                  outputs=['tempPositions'])

SelectMultiChainWizard().addTarget(protocol=ProtModellerComparativeModelling,
                                   targets=['tempMChain'],
                                   inputs=[{'templateOrigin': ['inputAtomStruct', 'pdbTemplate']}],
                                   outputs=['tempMChain'])

