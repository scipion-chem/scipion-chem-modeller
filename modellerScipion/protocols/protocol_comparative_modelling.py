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
This protocol is used to perform comparative modelling to predict a atomic structure from a sequence and set of
alignable known structures

"""

import os, json, shutil
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
import pwem.convert as emconv
from pwchem.utils.utilsFasta import parseAlnFile, parseFasta

from modellerScipion import Plugin

AUTOMODELLER, CLUSTALO, MUSCLE, MAFFT, CUSTOM = 'AutoModeller', 'Clustal_Omega', 'Muscle', 'Mafft', 'Custom'

class ModellerComparativeModelling(EMProtocol):
    """
    Performs a comparative modelling prediction using modeller and a set of similar structures
    https://salilab.org/modeller/manual/node16.html
    """
    _label = 'Comparative modelling'

    # -------------------------- DEFINE param functions ----------------------
    def _addTemplateForm(self, form):
        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Template structure: ", condition='templateOrigin==0',
                      help='Structure to use as template with a sequence alignment')

        form.addParam('pdbTemplate', params.StringParam, condition='templateOrigin==1',
                      label='Template PDB: ', help='Specify the PDB code to use as template')

        form.addParam('tempChain', params.StringParam,
                      label='Template chain: ',
                      help='Specify the protein chain to use as template. If None, modeller will try to guess it')
        form.addParam('tempPositions', params.StringParam,
                      label='Template positions: ',
                      help='Specify the positions of the sequence to use in the alignment. '
                           'If None, modeller will use the whole sequence')

        form.addParam('addTemplate', params.LabelParam,
                      label='Add template',
                      help='Add structure as template for the sequence alignment')

    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSequence', params.PointerParam,
                       pointerClass='Sequence', allowsNull=False,
                       label="Input sequence to predict: ",
                       help='Select the sequence whose atomic structure will be predicted')
        form.addParam('nModels', params.IntParam, default=1,
                      label="Number of models: ",
                      help='Number of models to generate')

        group = form.addGroup('Structure templates')
        group.addParam('templateOrigin', params.EnumParam, default=0,
                       label='Templates origin', choices=['AtomStruct', 'PDB code'],
                       help='Structure template origin')
        self._addTemplateForm(group)
        group.addParam('templateList', params.TextParam, width=100,
                      default='', label='List of templates: ',
                      help='The list of templates to use for the comparative modelling.')

        group = form.addGroup('Alignment')
        group.addParam('alignMethod', params.EnumParam,
                       label='Alignment method: ', default=0, 
                       choices=[AUTOMODELLER, CLUSTALO, MUSCLE, MAFFT, CUSTOM],
                       help="How to generate the sequences alignment:\n 1) Modeller automatic alignment\n"
                            "2,3,4) Scipion-chem included programs with default params\n"
                            "5) Custom alignment as SetOfSequences aligned")
        
        group.addParam('inputAlign', params.PointerParam, pointerClass='SetOfSequencesChem',
                       label='Input sequence alignment: ', condition='alignMethod == 4', allowsNull=True,
                       help="Input sequence alignment containing all the sequences specified "
                            "in this protocol formulary (both templates and target).\n")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('modellerStep')
        self._insertFunctionStep('createOutputStep')

    def modellerStep(self):
        pdbsFile = self.buildPDBsFile()
        alignFile = self.buildAlignFile()

        Plugin.runModeller(self, Plugin.getScriptsDir('comparative_modelling.py'),
                           args=self._getModellerArgs(), cwd=self._getPath())

    def createOutputStep(self):
        for file in os.listdir(self._getPath()):
            if file.endswith('.pdb'):
                outFile = self._getPath(file)
                outId = int(os.path.splitext(outFile)[0][-2:])

                modellerAS = AtomStruct(outFile)
                self._defineOutputs(**{'outputAtomStruct_{}'.format(outId): modellerAS})

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        scoresFile = self._getPath('scores.txt')
        if os.path.exists(scoresFile):
            summary.append('DOPE scores for the generated models (energy-like, the lower the better)\n')
            with open(scoresFile) as fSc:
              summary.append(fSc.read())
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    def _warnings(self):
        warns = []
        if self.getEnumText('alignMethod') == CUSTOM:
            warns.append('The custom input alignment must be conformed by the same sequences and positions defined in '
                         'this formulary, including the target and template sequences. You can build such a set of '
                         'sequences with a similar form with the "Define set of sequences" protocol followed by the '
                         '"multiple sequence alignment" protocol, where you can define and check the custom alignment')
        return warns
    
    # --------------------------- UTILS functions ------------------------
    def getTargetID(self):
        return self.inputSequence.get().getId()

    def getTargetSequence(self):
        return self.inputSequence.get().getSequence()

    def _getModellerArgs(self):
        args = ''
        args += '-i {} '.format(self.getTargetID())
        args += '-af {} '.format(os.path.abspath(self.getAlignmentFile()))
        args += '-pf {} '.format(os.path.abspath(self.getPDBsFile()))
        args += '-pd {} '.format(os.path.abspath(self._getExtraPath()))
        args += '-n {} '.format(self.nModels.get())
        if self.getEnumText('alignMethod') == AUTOMODELLER:
            args += '--align '
        return args.split()
        
    def getPDBsFile(self):
        return self._getExtraPath('templatePDBs.txt')

    def getAlignmentFile(self):
        return os.path.abspath(self._getPath('alignment.pir'))

    def buildPDBsFile(self):
        pdbsFile = self.getPDBsFile()
        with open(pdbsFile, 'w') as f:
            for tempLine in self.templateList.get().split('\n'):
              if tempLine.strip():
                tempJson = json.loads(tempLine.split(')')[1].strip())
                pdbCode = tempJson['pdbName']
                if 'pdbFile' in tempJson:
                    shutil.copy(tempJson['pdbFile'], self._getExtraPath(os.path.basename(tempJson['pdbFile'])))
                else:
                    aSH = emconv.AtomicStructHandler()
                    print("Retrieving atomic structure with ID = %s" % pdbCode)
                    aSH.readFromPDBDatabase(pdbCode, type='mmCif', dir=self._getExtraPath())

                f.write(pdbCode + '\n')
        return pdbsFile

    def buildAlignFile(self):
        alignFile = self.getAlignmentFile()
        seqDic = {}
        if self.getEnumText('alignMethod') == CUSTOM:
            seqDic = self.parseInputAlignment()
        elif self.getEnumText('alignMethod') in [CLUSTALO, MUSCLE, MAFFT]:
            seqDic = self.makeScipionAlignment()

        with open(alignFile, 'w') as f:
            targetSeq = self.getTargetSequence()
            if self.getTargetID() in seqDic:
                targetSeq = seqDic[self.getTargetID()]

            f.write('>P1;{}\nsequence::::::{}:::\n{}*\n'.
                    format(self.getTargetID(), self.inputSequence.get().getSeqName(), targetSeq))

            for tempLine in self.templateList.get().split('\n'):
              if tempLine.strip():
                  tempJson = json.loads(tempLine.split(')')[1].strip())
                  pdbCode, chain = tempJson['pdbName'], tempJson['chain']
                  idxs = tempJson['index'].split('-')

                  tempSeq = ''
                  if pdbCode in seqDic:
                      tempSeq = seqDic[pdbCode]


                  f.write('>P1;{}\nstructureX:{}:{}:{}:{}:{}:{}:::\n{}*\n'.
                          format(pdbCode, pdbCode, idxs[0], chain, idxs[1], chain, pdbCode, tempSeq))

        return alignFile

    def getScipionAlignFile(self):
        return self._getExtraPath('alignment.aln')

    def parseInputAlignment(self):
        seqDic = {}
        for seq in self.inputAlign.get():
            seqDic[seq.getId()] = seq.getSequence()
        return seqDic

    def makeScipionAlignment(self):
        inpSeqsFile = self._getTmpPath('inputSeqs.fa')
        with open(inpSeqsFile, 'w') as f:
            f.write('>{}\n{}\n'.
                    format(self.getTargetID(), self.getTargetSequence()))

            for tempLine in self.templateList.get().split('\n'):
              if tempLine.strip():
                  tempJson = json.loads(tempLine.split(')')[1].strip())
                  seqFile = tempJson['seqFile']
                  with open(seqFile) as fIn:
                    f.write(fIn.read().strip() + '\n')

        alignFile = self.getScipionAlignFile()
        programName = self.getEnumText('alignMethod')
        if programName == CLUSTALO:
            cline = 'clustalo -i {} --auto -o {} --outfmt=clu'.format(inpSeqsFile, alignFile)
        elif programName == MUSCLE:
            alignFile = alignFile.replace('.aln', '.fa')
            cline = 'muscle -align {} -output {}'.format(inpSeqsFile, alignFile)
        elif programName == MAFFT:
            cline = 'mafft --auto --clustalout {} > {}'.format(inpSeqsFile, alignFile)
        self.runJob(cline, '')

        if programName == MUSCLE:
            seqDic = parseFasta(alignFile)
        else:
            seqDic = parseAlnFile(alignFile)

        return seqDic



