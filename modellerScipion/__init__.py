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

# General imports
import os, subprocess

# Scipion em imports
from pyworkflow.utils import yellowStr
import pwchem
from scipion.install.funcs import InstallHelper

# Plugin imports
from .constants import MODELLER_DIC

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
		# Instantiating install helper
		installer = InstallHelper(MODELLER_DIC['name'], packageHome=cls.getVar(MODELLER_DIC['home']), packageVersion=MODELLER_DIC['version'])

		# Defining message to print after installation
		innerPath = os.path.join('lib', cls.getEnvName(MODELLER_DIC), 'modlib', 'modeller', 'config.py')
		fullPath = cls.getEnvPath(packageDictionary=MODELLER_DIC, innerPath=innerPath)
		message = '\\n\\nOnce you have obtained a license key, remember to write it in file {}\\n'.format(fullPath)
		message += 'This key can be obtained by registration in https://salilab.org/modeller/registration.html\\n'
		message = yellowStr(message)

		# Installing package
		installer.getCondaEnvCommand().addCondaPackages([f'modeller={MODELLER_DIC["version"]}'], channel='salilab')\
			.addCommand('echo -e "{}"'.format(message), 'INFO_MESSAGE_PRINTED')\
			.addPackage(env, dependencies=['conda'])

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