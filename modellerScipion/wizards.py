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

import json, os

from .protocols import ModellerMutateResidue, ModellerComparativeModelling
from pwem.wizards import SelectChainWizard, SelectResidueWizard, EmWizard, VariableWizard

from pwchem.wizards import SelectChainWizardQT, SelectResidueWizardQT

from .constants import AA_LIST

SelectChainWizard().addTarget(protocol=ModellerMutateResidue,
                              targets=['mutChain'],
                              inputs=['inputAtomStruct'],
                              outputs=['mutChain'])

SelectResidueWizard().addTarget(protocol=ModellerMutateResidue,
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
  
    def show(self, form, *params):
        protocol = form.protocol
        inputParams, outputParam = self.getInputOutput(form)

        # StructName
        inputTemplate = getattr(protocol, inputParams[0]).get()
        pdbFile = ''
        if issubclass(type(inputTemplate), str):
            outStr = [inputTemplate]
        else:
            pdbFile = inputTemplate.getFileName()
            outStr = [os.path.splitext(os.path.basename(pdbFile))[0]]

        # Chain
        chainJson = getattr(protocol, inputParams[1]).get()
        chainId = json.loads(chainJson)['chain']
        outStr += [chainId]

        # Positions
        posJson = getattr(protocol, inputParams[2]).get()
        if posJson:
            posIdxs = json.loads(posJson)['index']
            seq = json.loads(posJson)['residues']
            outStr += [posIdxs]
        else:
            outStr += ['FIRST-LAST']
            finalResiduesList = self.getResidues(form, inputParams)
            idxs = [json.loads(finalResiduesList[0].get())['index'], json.loads(finalResiduesList[-1].get())['index']]
            seq = self.getSequence(finalResiduesList, idxs)

        prevStr = getattr(protocol, outputParam[0]).get()
        lenPrev = len(prevStr.strip().split('\n')) + 1
        if prevStr.strip() == '':
            lenPrev -= 1
        elif not prevStr.endswith('\n'):
            prevStr += '\n'

        seqFile = protocol.getProject().getTmpPath('{}_{}_{}.fa'.format(outStr[0], lenPrev, outStr[2]))
        with open(seqFile, 'w') as f:
            f.write('>{}\n{}\n'.format(outStr[0], seq))

        pdbStr = ''
        if pdbFile:
            pdbStr = ', "pdbFile": "{}"'.format(pdbFile)
        jsonStr = '%s) {"pdbName": "%s", "chain": "%s", "index": "%s", "seqFile": "%s"%s}\n' % \
                  (lenPrev, outStr[0], outStr[1], outStr[2], seqFile, pdbStr)
        form.setVar(outputParam[0], prevStr + jsonStr)

AddStructSequenceWizard().addTarget(protocol=ModellerComparativeModelling,
                                    targets=['addTemplate'],
                                    inputs=[{'templateOrigin': ['inputAtomStruct', 'pdbTemplate']},
                                            'tempChain', 'tempPositions'],
                                    outputs=['templateList'])

SelectChainWizardQT().addTarget(protocol=ModellerComparativeModelling,
                                targets=['tempChain'],
                                inputs=[{'templateOrigin': ['inputAtomStruct', 'pdbTemplate']}],
                                outputs=['tempChain'])

SelectResidueWizardQT().addTarget(protocol=ModellerComparativeModelling,
                                  targets=['tempPositions'],
                                  inputs=[{'templateOrigin': ['inputAtomStruct', 'pdbTemplate']},
                                          'tempChain'],
                                  outputs=['tempPositions'])