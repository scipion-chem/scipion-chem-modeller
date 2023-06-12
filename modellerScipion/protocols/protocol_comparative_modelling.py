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

import os, json, shutil, glob, string
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
import pwem.convert as emconv

from pwchem import Plugin as pwchemPlugin
from pwchem.utils.utilsFasta import parseAlnFile, parseFasta
from pwchem.constants import BIOCONDA_DIC

from modellerScipion import Plugin
from modellerScipion.constants import MODELLER_DIC

AUTOMODELLER, CLUSTALO, MUSCLE, MAFFT, CUSTOM = 'AutoModeller', 'Clustal_Omega', 'Muscle', 'Mafft', 'Custom'
chainAlph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
scoreChoices = ['DOPE', 'DOPE-HR', 'Normalized_DOPE', 'GA341']

class ProtModellerComparativeModelling(EMProtocol):
    """
    Performs a comparative modelling prediction using modeller and a set of similar structures
    https://salilab.org/modeller/manual/node15.html
    """
    _label = 'Comparative modelling'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _addTemplateForm(self, form):
        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Template structure: ", condition='templateOrigin==0',
                      help='Structure to use as template with a sequence alignment')

        form.addParam('pdbTemplate', params.StringParam, condition='templateOrigin==1',
                      label='Template PDB: ', help='Specify the PDB code to use as template')

        form.addParam('tempChain', params.StringParam,
                      label='Template chain: ', condition='not multiChain',
                      help='Specify the protein chain to use as template. If None, modeller will try to guess it')

        form.addParam('tempPositions', params.StringParam,
                      label='Template positions: ', condition='not multiChain',
                      help='Specify the positions of the sequence to use in the alignment. '
                           'If None, modeller will use the whole sequence')

        form.addParam('tempMChain', params.StringParam,
                      label='Template multi chain: ', condition='multiChain',
                      help='Specify the protein chains to use as template (Use CTRL from multiple selection). '
                           'If None, modeller will try to guess it')
        form.addParam('tempMPositions', params.StringParam,
                      label='Template multi positions: ', condition='multiChain',
                      help='Specify the positions of each of the chains to use in the alignment.\n'
                           'In the same order as in "Template multi chain: ", write the first and last index '
                           'comma-separated. e.g: 1-30, 4-35, 20-55\n '
                           'If None, FIRST-LAST will be used')

        form.addParam('addTemplate', params.LabelParam,
                      label='Add template',
                      help='Add structure as template for the sequence alignment')

    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('multiChain', params.BooleanParam, default=False,
                      label="Multiple chains: ",
                      help='Build a comparative modelling using several chains. \nThis way, the input must be the '
                           'target sequences of each of the chains selected from the templates '
                           '(and in the same order!).\nYou can build this set of sequences using the '
                           '"define set of sequences" protocol')
        group.addParam('inputSequence', params.PointerParam,
                       pointerClass='Sequence', allowsNull=True,
                       label="Input sequence to predict: ", condition='not multiChain',
                       help='Select the sequence whose atomic structure will be predicted')
        group.addParam('inputSequences', params.PointerParam,
                      pointerClass='SetOfSequences', allowsNull=True,
                      label="Input sequences to predict: ", condition='multiChain',
                      help='Select the sequences whose atomic structure will be predicted')

        group.addParam('adIni', params.BooleanParam, default=False,
                      label="Use initial model: ",
                      help='Use structure as initial model. It must have the same sequence as the target')
        group.addParam('iniModel', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Initial model: ", condition='adIni',
                      help='Initial model to use for the target. It must have the same sequence as the target')

        group = form.addGroup('Output')
        group.addParam('nModels', params.IntParam, default=1,
                      label="Number of models: ",
                      help='Number of models to generate. Their generation can be parallelized.')
        group.addParam('modelH', params.BooleanParam, default=False,
                       label="Build model hydrogens: ",
                       help='Build also model hydrogens')

        group = form.addGroup('Structure templates')
        group.addParam('templateOrigin', params.EnumParam, default=0,
                       label='Templates origin: ', choices=['AtomStruct', 'PDB code'],
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

        group.addParam('fromFile', params.BooleanParam, default=True,
                      label="Alignment from file: ", condition='alignMethod==4 and not multiChain',
                      help='Whether to choose custom alignment from file or scipion alignment object')
        group.addParam('inputAlign', params.PointerParam, pointerClass='SetOfSequencesChem', allowsNull=True,
                       label='Input sequence alignment: ',
                       condition='alignMethod==4 and not fromFile and not multiChain',
                       help="Input sequence alignment containing all the sequences specified "
                            "in this protocol formulary (both templates and target).\n")

        group.addParam('inputAlignFile', params.PathParam,
                       label='Input alignment file: ', condition='alignMethod==4 and fromFile', allowsNull=True,
                       help="Input sequence alignment containing all the sequences and chains in PIR format as"
                            " expected from modeller. https://salilab.org/modeller/manual/node29.html")


        form.addSection(label='Other parameters')
        group = form.addGroup('Renaming', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('renumberResidues', params.StringParam, default='', label="Renumber residues: ",
                       help='Renumber residues so each chain first residue index is the one specified.'
                            'You can either enter one value and every chain first index will be that or specify '
                            'comma-separated values for each of the chains. Blank will not renumber.')
        group.addParam('renameChains', params.StringParam, default='', label="Rename chains: ",
                       help='Rename chains with desired values. Must be comma separated values for the chain names, '
                            'in the order of the sequences. If you use symmetry, remember to keep using these new '
                            'names to refer to the chains. Blank will not rename.')

        group = form.addGroup('Symmetry', condition='multiChain')
        group.addParam('symChains', params.StringParam, default='',
                       label='Symmetry chains: ', condition='multiChain',
                       help='Specify the chains to model with symmetry. Comma separated chain pairs, in alphabetic '
                            'order corresponding to the input sequences. \ni.e: symmetry for chains first and third;'
                            ' second and fourth (hemoglobin), parameter must be: *A-C, B-D*')
        group.addParam('symAtom', params.StringParam, default='CA',
                       label='Symmetry atoms: ', condition='multiChain',
                       help='In order to accelerate the symmetry calculation, the symmetry distances will be '
                            'calculated over this type of atoms. Def: CA (alpha carbons)')

        group = form.addGroup('Scoring')
        for i, scoreName in enumerate(scoreChoices):
            defa = True if i == 0 else False
            group.addParam(f'score{scoreName}', params.BooleanParam, default=defa,
                          label=f"Report {scoreName} score: ")

        group = form.addGroup('Optimization')
        group.addParam('opt', params.EnumParam,
                       label='Optimization quality: ', default=1,
                       choices=['Low-Fast', 'Default', 'High-Slow'],
                       help="Modeller allows to modify how fast the optimization will be performed. Fast optimizations"
                            "will usually lead to lower quality.\nhttps://salilab.org/modeller/manual/node19.html")
        group.addParam('nReps', params.IntParam, default=1,
                       label="Number of optimization cycles: ",
                       help='Number of optimization cycles, including the energy optimization and Molecular Dynamics')
        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('modellerStep')
        self._insertFunctionStep('createOutputStep')

    def alignStep(self):
        alignFile = self.buildAlignFile()

    def modellerStep(self):
        pdbsFile = self.buildPDBsFile()
        Plugin.runScript(self, 'comparative_modelling.py', args=self._getModellerArgs(),
                         envDic=MODELLER_DIC, cwd=self._getPath())

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
            summary.append('GA341 score ranges from 0 to 1, the higher the better\n')
            summary.append('Rest of scores for the generated models are energy-like, the lower the better)\n')
            with open(scoresFile) as fSc:
              summary.append(fSc.read())
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        if self.adIni and not self.iniModel.get():
            errors.append('You have not specified the initial model')
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
    def getTargetID(self, seqObj=None):
        if not seqObj:
            if not self.multiChain:
                seqObj = self.inputSequence.get()
            else:
                for seqObj in self.inputSequences.get():
                    pass
        return seqObj.getId()

    def getTargetSequence(self, seqObj=None):
        if not seqObj:
            seqObj = self.inputSequence.get()
        return seqObj.getSequence()

    def _getModellerArgs(self):
        args = ''
        args += '-i {} '.format(self.getTargetID())
        args += '-af {} '.format(os.path.abspath(self.getAlignmentFile()))
        args += '-pf {} '.format(os.path.abspath(self.getPDBsFile()))
        args += '-pd {} '.format(os.path.abspath(self._getExtraPath()))
        if self.opt.get() != 0:
            args += '-n {} '.format(self.nModels.get())
        else:
            args += '-n 1 '
            print('With {} optimization, the initial model is not randomized so every output model is the same.\n'
                  'Therefore, only one model is output'.format(self.getEnumText('opt')))
        if self.getEnumText('alignMethod') == AUTOMODELLER:
            args += '--align '

        if self.adIni:
            args += '-im {} '.format(os.path.abspath(self.iniModel.get().getFileName()))
        if self.modelH:
            args += '--modelH '

        doScore = []
        for scoreName in scoreChoices:
            if getattr(self, f'score{scoreName}'):
                doScore.append(scoreName)
        if len(doScore) > 0:
            args += '-sc {} '.format(','.join(doScore))

        if self.opt.get() != 1:
            args += '-opt {} '.format(self.getEnumText('opt'))
        args += '-nr {} '.format(self.nReps.get())

        nChains = 1 if not self.multiChain.get() else len(self.inputSequences.get())

        renamed = False
        renameStr = self.renameChains.get().strip()
        if renameStr and len(renameStr.split(',')) == nChains:
          args += '-renam {} '.format(renameStr)
          renamed = True

        renumStr = self.renumberResidues.get().strip()
        if renumStr:
            renumStr = renumStr if len(renumStr.split(',')) == nChains else ','.join([renumStr] * nChains)
            args += '-renum {} '.format(renumStr)
            if not renamed:
                args += '-renam {} '.format(','.join(list(string.ascii_uppercase)[:nChains]))


        if self.symChains.get() and self.symChains.get().strip():
            symChains = self.symChains.get().replace(' ', '')
            args += '-sym {} '.format(symChains)
            args += '-symAtom {} '.format(self.symAtom.get())

        args += '-nj {} '.format(self.numberOfThreads.get())
        args += '-mPath {} '.format(Plugin.getPluginHome())

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
                    shutil.copy(tempJson['pdbFile'], self._getExtraPath(os.path.basename(tempJson['pdbFile']).lower()))
                else:
                    aSH = emconv.AtomicStructHandler()
                    aSH.readFromPDBDatabase(pdbCode, type='mmCif', dir=self._getExtraPath())

                f.write(pdbCode + '\n')
        return pdbsFile

    def buildAlignFile(self):
        alignFile = self.getAlignmentFile()
        programName = self.getEnumText('alignMethod')

        if not self.multiChain:
            if programName == CUSTOM:
                if not self.fromFile.get():
                    seqDic = self.parseInputAlignment()
                else:
                  shutil.copy(self.inputAlignFile.get(), alignFile)
                  return alignFile

            elif programName in [CLUSTALO, MUSCLE, MAFFT]:
                seqDic = self.makeScipionAlignment(programName)

            else:
                # todo: automatic alignment, meanwhile use mafft
                seqDic = self.makeScipionAlignment(MAFFT)

            with open(alignFile, 'w') as f:
                targetSeq = self.getTargetSequence()
                if self.getTargetID() in seqDic:
                    targetSeq = seqDic[self.getTargetID()]

                f.write('>P1;{}\nsequence:::A:::{}:::\n{}*\n'.
                        format(self.getTargetID(), self.inputSequence.get().getSeqName(), targetSeq))

                for tempLine in self.templateList.get().split('\n'):
                  if tempLine.strip():
                      tempJson = json.loads(tempLine.split(')')[1].strip())
                      pdbCode, chain = tempJson['pdbName'], tempJson['chain']
                      idxs = tempJson['index'].split('-')

                      tempSeq = seqDic[pdbCode]
                      f.write('>P1;{}\nstructureX:{}:{}:{}:{}:{}:{}:::\n{}*\n'.
                              format(pdbCode, pdbCode, idxs[0], chain, idxs[1], chain, pdbCode, tempSeq))

        else:
          if programName == CUSTOM:
              shutil.copy(self.inputAlignFile.get(), alignFile)
              return alignFile

          elif programName in [CLUSTALO, MUSCLE, MAFFT]:
              seqDic = self.makeScipionAlignment(programName)

          else:
              # todo: automatic alignment, meanwhile use mafft
              seqDic = self.makeScipionAlignment(MAFFT)

          with open(alignFile, 'w') as f:
              seqStr = ''
              for i, targetSeqObj in enumerate(self.inputSequences.get()):
                  targetSeq = self.getTargetSequence(targetSeqObj).strip()
                  seqStr += '{}/'.format(targetSeq)
              targetID = self.getTargetID(targetSeqObj)

              f.write('>P1;{}\nsequence:::A::{}:{}:::\n{}*\n'.
                      format(targetID, chainAlph[i], targetID, seqStr[:-1]))

              for tempLine in self.templateList.get().split('\n'):
                if tempLine.strip():
                  tempJson = json.loads(tempLine.split(')')[1].strip())
                  pdbCode = tempJson['pdbName']

                  allChains = tempJson['chains'].split(',')
                  chains = [allChains[0].split('-')[1], allChains[-1].split('-')[1]]

                  seqDic.pop(targetID)
                  chainsSeqs = '/'.join(seqDic.values())

                  f.write('>P1;{}\nstructureX:{}::{}::{}:{}:::\n{}*\n'.
                          format(pdbCode, pdbCode, chains[0], chains[1], pdbCode, chainsSeqs))

        return alignFile

    def getScipionAlignFile(self, idx=''):
        return os.path.abspath(self._getExtraPath('alignment{}.aln'.format(idx)))

    def parseInputAlignment(self):
        seqDic = {}
        for seq in self.inputAlign.get():
            seqDic[seq.getId()] = seq.getSequence()
        return seqDic

    def makeScipionAlignment(self, programName):
        if not self.multiChain:
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

            seqDic = self.performAlignment(inpSeqsFile, programName)
            return seqDic

        else:
            seqDic = {}
            for i, inSeq in enumerate(self.inputSequences.get()):
                inpSeqsFile = self._getTmpPath('inputSeqs_{}.fa'.format(i))
                with open(inpSeqsFile, 'w') as f:
                    f.write('>{}\n{}\n'.format(self.getTargetID(inSeq), self.getTargetSequence(inSeq)))

                    for tempLine in self.templateList.get().split('\n'):
                        if tempLine.strip():
                          tempJson = json.loads(tempLine.split(')')[1].strip())
                          seqFileTemplate = tempJson['seqFiles']
                          chainId = tempJson['chains'].split(',')[i].strip().split('-')[1]
                          seqFileTemplate = seqFileTemplate.replace('_*_*', '_{}_*'.format(chainId))
                          seqFile = glob.glob(seqFileTemplate)[0]
                          with open(seqFile) as fIn:
                            f.write(fIn.read().strip() + '\n')

                iSeqDic = self.performAlignment(inpSeqsFile, programName, idx=str(i))
                seqDic.update(iSeqDic)

            return seqDic


    def performAlignment(self, inpFile, programName, idx=''):
        alignFile = self.getScipionAlignFile(idx)
        cline = '%s && ' % (pwchemPlugin.getEnvActivationCommand(BIOCONDA_DIC))
        if programName == CLUSTALO:
          cline += 'clustalo -i {} --auto -o {} --outfmt=clu'.format(inpFile, alignFile)
        elif programName == MUSCLE:
          alignFile = alignFile.replace('.aln', '.fa')
          cline += 'muscle -align {} -output {}'.format(inpFile, alignFile)
        elif programName == MAFFT:
          cline += 'mafft --auto --clustalout {} > {}'.format(inpFile, alignFile)
        self.runJob(cline, '')

        seqIds = list(parseFasta(inpFile).keys())
        if programName == MUSCLE:
          seqDic = parseFasta(alignFile)
        else:
          seqDic = parseAlnFile(alignFile)

        nSeqDic = {}
        for i, old_key in enumerate(seqDic):
            nSeqDic[seqIds[i]] = seqDic[old_key]

        return nSeqDic
