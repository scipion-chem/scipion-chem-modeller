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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from ..protocols import ModellerMutateResidue
from ..constants import AA_LIST

textMutationListExample = '{"model": 0, "chain": "A", "residues": 141} | {"index": "1-1", "residues": "V"} | TRP\n\
{"model": 0, "chain": "B", "residues": 146} | {"index": "2-2", "residues": "H"} | PRO\n'

class TestModellerMutateRes(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    def _runModellerMutate(self):
        PRO = AA_LIST.index('PRO')
        protModeller = self.newProtocol(
            ModellerMutateResidue,
            inputAtomStruct=self.protImportPDB.outputPdb,
            mutChain='{"model": 0, "chain": "B", "residues": 146}',
            mutPosition='{"index": "2-2", "residues": "H"}', mutResidue=PRO,
            toMutateList=textMutationListExample)

        self.launchProtocol(protModeller)
        pdbOut = getattr(protModeller, 'mutatedAtomStruct', None)
        self.assertIsNotNone(pdbOut)

    def test_mutateResidue(self):
        self._runModellerMutate()




