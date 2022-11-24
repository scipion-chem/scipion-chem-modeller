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

import os, json
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct

from modellerScipion import Plugin
from ..constants import AA_LIST

class ModellerMutateResidue(EMProtocol):
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
                    help='Here you can define a mutation which will be added to the list of mutations below.'
                         'Modeller will be used to sequentially perform the mutations you define in the list.')

    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', params.PointerParam,
                        pointerClass='AtomStruct', allowsNull=False, label="Input atom structure",
                        help='Select the atom structure to be fitted in the volume')

        group = form.addGroup('Define mutation')
        self._addMutationForm(group)
        group.addParam('toMutateList', params.TextParam, width=70,
                      default='', label='List of mutations',
                      help='List of chain | position | residue to mutate. '
                           'The mutations here will be performed sequentially using the modeller software.')
        group.addParam('clearLabel', params.LabelParam,
                      label='Clear mutation list',
                      help='Clear mutations list')

        form.addParam('seed', params.IntParam, label='Random seed', expertLevel=params.LEVEL_ADVANCED,
                      default=-49837, help='Random seed for modeller')

        form.addSection(label='Energy objective functions')
        group = form.addGroup('Distance parameters')
        group.addParam('contactShell', params.FloatParam, default=4.0,
                       label='Nonbond distance cutoff',
                       help='This defines the maximal distance between atoms that flags a non-bonded atom pair. '
                            'Such pairs are stored in the list of non-bonded atom pairs. Only those non-bonded pairs '
                            'that are sufficiently close to each other will result in an actual non-bonded restraint. '
                            'https://salilab.org/modeller/9.9/manual/node122.html')
        group.addParam('updateDynamic', params.FloatParam, default=0.39,
                       label='Nonbond recalculation threshold',
                       help='This sets the cumulative maximal atomic shift during optimization that triggers '
                            'recalculation of the list of atom-atom non-bonded pairs. It should be set in combination '
                            'with contact_shell. https://salilab.org/modeller/9.9/manual/node123.html')

        group = form.addGroup('Energy restrains')
        group.addParam('dynamicSphere', params.BooleanParam, default=True,
                       label='Calculate soft-sphere overlap restraints',
                       help='The dynamic soft-sphere overlap restraints are calculated. Note that they are only '
                            'calculated if the scaled standard deviation of the soft-sphere overlap restraints is '
                            'greater than zero. The soft-sphere potential is simply a lower bound harmonic restraint '
                            'https://salilab.org/modeller/9.9/manual/node125.html')
        group.addParam('sphereStdv', params.FloatParam, default=0.05, condition='dynamicSphere==True',
                       label='Soft-sphere standard deviation', expertLevel=params.LEVEL_ADVANCED,
                       help="The dynamic soft-sphere overlap restraints are calculated. "
                            "The soft-sphere potential is simply a lower bound harmonic restraint, "
                            "with standard deviation sphereStdv, dropping to zero at the sum of "
                            "the two atoms' van der Waals radii: https://salilab.org/modeller/9.9/manual/node124.html")
        group.addParam('lennardJones', params.BooleanParam, default=False,
                      label='Calculate Lennard-Jones restraints',
                      help='Dynamic Lennard-Jones restraints are calculated, using equation: '
                           'https://salilab.org/modeller/9.9/manual/node459.html#eq:lennard')
        group.addParam('LJSwitch1', params.FloatParam, default=6.5, condition='lennardJones==True',
                       label='Lennard-Jones switching parameter 1', expertLevel=params.LEVEL_ADVANCED,
                       help="This is the parameter f_1 to the Lennard-Jones switching function, "
                            "which smoothes the potential down to zero")
        group.addParam('LJSwitch2', params.FloatParam, default=7.5, condition='lennardJones==True',
                       label='Lennard-Jones switching parameter 2', expertLevel=params.LEVEL_ADVANCED,
                       help="This is the parameter f_2 to the Lennard-Jones switching function, "
                            "which smoothes the potential down to zero")
        group.addParam('coulomb', params.BooleanParam, default=False,
                       label='Calculate Coulomb restraints',
                       help='Dynamic Coulomb (electrostatic) restraints are calculated '
                            'https://salilab.org/modeller/9.9/manual/node459.html#eq:coulomb')
        group.addParam('CouSwitch1', params.FloatParam, default=6.5, condition='coulomb==True',
                       label='Coulomb switching parameter 1', expertLevel=params.LEVEL_ADVANCED,
                       help="This is the parameter f_1 to the Coulomb switching function, "
                            "which smoothes the potential down to zero")
        group.addParam('CouSwitch2', params.FloatParam, default=7.5, condition='coulomb==True',
                       label='Coulomb switching parameter 2', expertLevel=params.LEVEL_ADVANCED,
                       help="This is the parameter f_2 to the Coulomb switching function, "
                            "which smoothes the potential down to zero")
        group.addParam('relativeDielectric', params.FloatParam, default=1.0, condition='coulomb==True',
                       label='Relative dielectric', expertLevel=params.LEVEL_ADVANCED,
                       help="This sets the relative dielectric epsilon_r , used in the calculation"
                            "of the Coulomb energy")
        group.addParam('dynamicModeller', params.BooleanParam, default=False,
                       label='Calculate non-bonded spline restraints',
                       help='Dynamic MODELLER non-bonded spline restraints are calculated. These include the loop '
                            'modeling potential and DOPE: https://salilab.org/modeller/9.9/manual/node128.html')



    def _parseChain(self, chainLine):
      return json.loads(chainLine)['chain']

    def _parsePosition(self, posLine):
      return json.loads(posLine)['index'].split('-')[0]

    def _parseType(self, typeLine):
      return typeLine.strip()

    def parseMutations(self):
      chains, respos, restypes = [], [], []
      text = self.toMutateList.get()
      for line in text.split('\n'):
        if len(line.split('|'))==3:
          chain, resp, restyp = line.split('|')
          chains.append(self._parseChain(chain.strip()))
          respos.append(self._parsePosition(resp.strip()))
          restypes.append(self._parseType(restyp.strip()))
      return chains, respos, restypes

    def _getModellerArgs(self):
      self.outputFolder, self.ASFile = os.path.abspath(self._getPath()), self._getFileInputStruct()
      modelbase, ext = os.path.splitext(self.ASFile.split('/')[-1])
      self.outputFile = '{}/{}_mutant.pdb'.format(self.outputFolder, modelbase)

      chains, respos, restypes = self.parseMutations()

      args = ['-i', self.ASFile, '-p', respos, '-r', restypes, '-c', chains, '-s', self.seed.get(),
              '-o', self.outputFile]

      args += ['-contactShell', self.contactShell.get(), '-updateDynamic', self.updateDynamic.get()]
      if self.dynamicSphere.get():
        args += ['--dynamicSphere', '-sphereStdv', self.sphereStdv.get()]
      if self.lennardJones.get():
        args += ['--lennardJones', '-LJSwitch1', self.LJSwitch1.get(), '-LJSwitch2', self.LJSwitch2.get()]
      if self.coulomb.get():
        args += ['--coulomb', '-CouSwitch1', self.CouSwitch1.get(), '-CouSwitch2', self.CouSwitch2.get(),
                 '-relativeDielectric', self.relativeDielectric.get()]
      if self.dynamicModeller.get():
        args += ['--dynamicModeller']

      return args

    def _getScriptsFolder(self):
      thisFile = os.path.realpath(__file__)
      path = []
      for step in thisFile.split('/'):
        path.append(step)
        if step == 'scipion-chem-modeller':
          break
      return '/'+os.path.join(*path)+'/modellerScipion/scripts-10_1/'

    def _getFileInputStruct(self):
      structFile = self.inputAtomStruct.get().getFileName()
      return os.path.abspath(structFile)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('modellerStep')
        self._insertFunctionStep('createOutputStep')

    def modellerStep(self):
        Plugin.runModeller(self, self._getScriptsFolder()+'mutate_residue.py',
                          args=self._getModellerArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        mutatedAS = AtomStruct(self.outputFile)
        self._defineOutputs(mutatedAtomStruct=mutatedAS)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        if not self.toMutateList.get().strip():
            errors.append('You have not added any mutation to the list. Do so using the "add mutation" '
                  'wizard once you have defined it')
        else:
            for i, line in enumerate(self.toMutateList.get().strip().split('\n')):
                residue = line.split('|')[1].strip()
                idxs = json.loads(residue)['index']
                if idxs.split('-')[0] != idxs.split('-')[1]:
                    errors.append('Error in mutation nº {}: '
                                  'Modeller protocol designed to produce one point substitutions.'.format(i+1))
        return errors
