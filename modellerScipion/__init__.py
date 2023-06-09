# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

import os, subprocess

from pyworkflow.utils import yellowStr
import pwchem

from .constants import *

_version_ = '0.1'
_logo = "modeller_logo.png"
_references = ['Webb2016']

class Plugin(pwchem.Plugin):
    _supportedVersions = [MODELLER_DIC['version']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. """
        cls._defineEmVar(MODELLER_DIC['home'], '{}-{}'.format(MODELLER_DIC['name'], MODELLER_DIC['version']))

    @classmethod
    def defineBinaries(cls, env):
        cls.addModellerPackage(env)
    
    @classmethod
    def addModellerPackage(cls, env):
        """ This function installs Modeller package. """
        MODELLER_INSTALLED = 'modeller_installed'

        installationCmd = cls.getCondaActivationCmd()
        installationCmd += ' conda create -y -c salilab -n {} modeller={} &&'.format(cls.getEnvName(MODELLER_DIC), MODELLER_DIC['version'])
        installationCmd += ' touch %s &&' % MODELLER_INSTALLED
        installationCmd += ' {} && export MODELLER_CONDA_PREFIX=$CONDA_PREFIX &&'.format(cls.getEnvActivationCommand(MODELLER_DIC, condaHook=False))

        asciiYellowCode = '\033[0;33m'
        asciiNoColorCode = '\033[0m'
        stringBefore = f'{asciiYellowCode}Once Modeller is downloaded and installed, remember to write the license key in file '
        stringAfter = f'/modlib/modeller/config.py in the conda environment. This key can be obtained by registration in https://salilab.org/modeller/registration.html{asciiNoColorCode}'
        installationCmd += ' echo "{}$MODELLER_CONDA_PREFIX{}" && echo $MODELLER_CONDA_PREFIX'.format(stringBefore, stringAfter)

        env.addPackage(MODELLER_DIC['name'],
                       version=MODELLER_DIC['version'],
                       tar='void.tgz',
                       commands=[(installationCmd, MODELLER_INSTALLED)],
                       neededProgs=cls.getDependencies(),
                       default=True)

    @classmethod
    def runScript(cls, protocol, scriptName, args, envDic, cwd=None, popen=False):
        """ Run modeller command from a given protocol. """
        scriptName = cls.getScriptsDir(scriptName)
        fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(envDic), 'python', scriptName)
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(fullProgram + args, cwd=cwd, shell=True)

    # ---------------------------------- Utils functions  -----------------------

    @classmethod
    def getPluginHome(cls, path=""):
        import modellerScipion
        fnDir = os.path.split(modellerScipion.__file__)[0]
        return os.path.join(fnDir, path)

    @classmethod
    def getScriptsDir(cls, scriptName=''):
        return cls.getPluginHome('scripts/%s' % scriptName)

    @classmethod
    def getDependencies(cls):
        # try to get CONDA activation command
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs
