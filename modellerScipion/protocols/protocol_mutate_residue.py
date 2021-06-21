# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to perform a residue mutation in a protein structure.
A energy optimization is performed over the mutated residue and its surroundings.

"""
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb
from pyworkflow.utils import Message
import pyworkflow.utils as pwutils
import os
from pwem.objects.data import AtomStruct
from modellerScipion import Plugin
from ..constants import AA_LIST

class modellerMutateResidue(EMProtocol):
    """
    Performs a residue substitution in a protein structure.
    https://salilab.org/modeller/wiki/Mutate%20model
    """
    _label = 'Mutate structure residue'

    # -------------------------- DEFINE param functions ----------------------
    def _addMutationForm(self, form):
      form.addParam('mutChain', params.StringParam,
                    allowsNull=False, label='Chain to mutate',
                    help='Specify the protein chain to mutate')
      form.addParam('mutPosition', params.StringParam,
                    allowsNull=False, label='Position to mutate',
                    help='Specify the residue position to mutate (int)')
      form.addParam('mutResidue', params.EnumParam,
                    choices=AA_LIST, label="Residue to introduce",
                    help='Select the substitute residue which will be introduced')
      form.addParam('addMutation', params.LabelParam,
                    label='Add defined mutation',
                    help='Add defined mutation to list')

    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputAtomStruct', params.PointerParam,
                       pointerClass='AtomStruct', allowsNull=False,
                       label="Input atom structure",
                       help='Select the atom structure to be fitted in the volume')
        group = form.addGroup('Add mutation')
        self._addMutationForm(group)
        group = form.addGroup('List of mutations')

        group.addParam('toMutateList', params.TextParam, width=70,
                      default='', label='List of mutations',
                      help='List of chain | position | residue to mutate. '
                           'They will be muted serially.')
        group.addParam('clearLabel', params.LabelParam,
                      label='Clear mutation list',
                      help='Clear mutations list')

        form.addParam('seed', params.IntParam, label='Random seed', expertLevel=params.LEVEL_ADVANCED,
                      default=-49837, help='Random seed for modeller')

    def _parseChain(self, chainLine):
      return chainLine.split(',')[1].split(':')[1].strip().split('"')[1]

    def _parsePosition(self, posLine):
      return posLine.split(":")[1].split(",")[0].strip()

    def _parseType(self, typeLine):
      return typeLine.strip()

    def parseMutations(self):
      chains, respos, restypes = [], [], []
      text = self.toMutateList.get()
      for line in text.split('\n'):
        if len(line.split('|'))==3:
          chain, resp, restyp = line.split('|')
          chains.append(self._parseChain(chain))
          respos.append(self._parsePosition(resp))
          restypes.append(self._parseType(restyp))
      return chains, respos, restypes

    def _getModellerArgs(self):
      self.outputFolder, self.pdbFile = os.path.abspath(self._getPath()), self._getPdbInputStruct()
      modelbase, ext = os.path.splitext(self.pdbFile.split('/')[-1])
      self.outputFile = '{}/{}_mutant.pdb'.format(self.outputFolder, modelbase)

      chains, respos, restypes = self.parseMutations()

      args = ['-i', self.pdbFile, '-p', respos, '-r', restypes, '-c', chains,
              '-o', self.outputFile]
      return args

    def _getScriptsFolder(self):
      thisFile = os.path.realpath(__file__)
      path = []
      for step in thisFile.split('/'):
        path.append(step)
        if step == 'scipion-chem-modeller':
          break
      return '/'+os.path.join(*path)+'/modellerScipion/scripts-10_1/'

    def _getPdbInputStruct(self):
      inpStruct = self.inputAtomStruct.get()
      name, ext = os.path.splitext(inpStruct.getFileName())
      if ext != '.pdb':
          cifFile = inpStruct.getFileName()
          pdbFile = self._getExtraPath(pwutils.replaceBaseExt(cifFile, 'pdb'))
          toPdb(cifFile, pdbFile)

      else:
        pdbFile = inpStruct.getFileName()
      return os.path.abspath(pdbFile)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('modellerStep')
        self._insertFunctionStep('createOutputStep')

    def modellerStep(self):
        Plugin.runModeller(self, self._getScriptsFolder()+'mutate_residue.py',
                          args=self._getModellerArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        mutatedPDB = AtomStruct(self.outputFile)
        self._defineOutputs(mutatedPDB=mutatedPDB)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods
