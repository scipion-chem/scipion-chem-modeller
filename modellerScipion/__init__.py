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
import pwem

from .constants import *


_version_ = '0.1'
_logo = "modeller_logo.png"
_references = ['Webb2016']

MODELLER_DIC = {'name': 'modeller', 'version': '10.3', 'home': 'MODELLER_HOME'}

class Plugin(pwem.Plugin):
    _supportedVersions = [MODELLER_DIC['version']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineVar("MODELLER_ENV_ACTIVATION", 'conda activate modeller-env')

    @classmethod
    def defineBinaries(cls, env):
        MODELLER_INSTALLED = 'modeller_installed'

        installationCmd = cls.getCondaActivationCmd()
        installationCmd += ' conda create -y -c salilab -n modeller-env modeller={} &&'.format(MODELLER_DIC['version'])
        installationCmd += ' touch %s' % MODELLER_INSTALLED

        env.addPackage(MODELLER_DIC['name'],
                       version=MODELLER_DIC['version'],
                       tar='void.tgz',
                       commands=[(installationCmd, MODELLER_INSTALLED)],
                       neededProgs=cls.getDependencies(),
                       default=True)

        print(yellowStr('\nOnce Modeller is downloaded and installed, remember to write the license key '
                        'in file {}/envs/modeller-env/lib/modeller-{}/modlib/modeller/config.py in the conda '
                        'environment. This key can be obtained by registration in '
                        'https://salilab.org/modeller/registration.html\n'.
                        format(cls.getCondaEnvPath(), MODELLER_DIC['version'])))

    @classmethod
    def runScript(cls, protocol, scriptName, args, env, cwd=None, popen=False):
        """ Run modeller command from a given protocol. """
        scriptName = cls.getScriptsDir(scriptName)
        fullProgram = '%s %s && %s %s' % (cls.getCondaActivationCmd(), cls.getEnvActivation(env), 'python', scriptName)
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(fullProgram + args, cwd=cwd, shell=True)

    # ---------------------------------- Utils functions  -----------------------

    @classmethod
    def getEnvActivation(cls, env):
        activation = cls.getVar("{}_ENV_ACTIVATION".format(env.upper()))
        return activation

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

    @classmethod
    def getCondaEnvPath(cls):
        condaAct = cls.getCondaActivationCmd()
        try:
            condaPath = condaAct.split('$(')[1].split('/bin')[0]
        except:
            condaPath = 'YOUR_CONDA_PATH'
        return condaPath
