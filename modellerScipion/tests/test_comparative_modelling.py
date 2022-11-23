# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import json

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols.protocol_import import ProtImportSequence
import pwem.convert as emconv
import pwem.objects as emobj
from pwem.constants import RESIDUES3TO1, RESIDUES1TO3

from pwchem.utils import downloadPDB

from ..protocols import ProtModellerComparativeModelling

pdbDic = {'1a00': ['56', 'B', 'FIRST-LAST'], '1a01': ['496', 'D', 'FIRST-LAST'],
          '1a0u': ['356', 'B', 'FIRST-LAST'], '1aj9': ['613', 'B', 'FIRST-LAST'],}

templatesStr = ''
for i, pdbId in enumerate(pdbDic):
    templatesStr += '%s) {"pdbName": "%s", "chain": "%s", "index": "%s", "seqFile": "Tmp/%s_%s_%s_%s.fa"}\n' \
                    % (i+1, pdbId, pdbDic[pdbId][1], pdbDic[pdbId][2], pdbId, *pdbDic[pdbId])


def getResidueList(modelsFirstResidue, model, chain):
    residueList = []
    for modelID, chainDic in modelsFirstResidue.items():
        if int(model) == modelID:
            for chainID, seq_number in chainDic.items():
                if chain == chainID:
                    for i in seq_number:
                        residueList.append('{"index": %d, "residue": "%s"}' % (i[0], str(i[1])))
    return residueList

def getSequence(finalResiduesList, idxs):
    roiStr, inSeq = '', False
    for residue in finalResiduesList:
        resDic = json.loads(residue.get())
        if resDic['residue'] in RESIDUES3TO1:
            if resDic['index'] == idxs[0]:
              inSeq = True
            if resDic['index'] == idxs[-1]:
              inSeq = False
              roiStr += RESIDUES3TO1[resDic['residue']]
              break

            if inSeq:
                roiStr += RESIDUES3TO1[resDic['residue']]
    return roiStr

class TestModellerComparativeModelling(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportSeq()
        cls._writeSequenceTemplates()

    @classmethod
    def _writeSequenceTemplates(cls):
        structureHandler = emconv.AtomicStructHandler()

        for pdbId in pdbDic:
            fileName = downloadPDB(pdbId, structureHandler)

            structureHandler.read(fileName)
            modelsLength, modelsFirstResidue = structureHandler.getModelsChains()

            residueList = getResidueList(modelsFirstResidue, 0, 'B')
            finalResiduesList = []
            for i in residueList:
                finalResiduesList.append(emobj.String(i))

            idxs = [json.loads(finalResiduesList[0].get())['index'],
                    json.loads(finalResiduesList[-1].get())['index']]
            seq = getSequence(finalResiduesList, idxs)

            seqFile = cls.proj.getTmpPath('{}_{}_{}_{}.fa'.format(pdbId, *pdbDic[pdbId]))
            with open(seqFile, 'w') as f:
                f.write('>{}\n{}\n'.format(pdbId, seq))


    @classmethod
    def _runImportSeq(cls):
        protImportSeq = cls.newProtocol(
            ProtImportSequence,
            inputSequence=0, inputProteinSequence=3,
            uniProtSequence='P68871', inputSequenceName='hemo')
        cls.launchProtocol(protImportSeq)
        cls.protImportSeq = protImportSeq

    def _runModellerComparative(self):
        protModeller = self.newProtocol(
            ProtModellerComparativeModelling,
            inputSequence=self.protImportSeq.outputSequence,
            alignMethod=3, templateList=templatesStr)

        self.launchProtocol(protModeller)
        pdbOut = getattr(protModeller, 'outputAtomStruct_1', None)
        self.assertIsNotNone(pdbOut)

    def test_mutateResidue(self):
        self._runModellerComparative()




