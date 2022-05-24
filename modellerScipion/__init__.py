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

import pwem
from os.path import join
from .constants import *
from pyworkflow.utils import yellowStr

_version_ = '0.1'
_logo = "modeller_icon.png"
_references = ['Webb2016']

MODELLER_DIC = {'name': 'modeller', 'version': '10.1', 'home': 'MODELLER_HOME'}

class Plugin(pwem.Plugin):
    _homeVar = MODELLER_DIC['home']
    _pathVars = [MODELLER_DIC['home']]
    _supportedVersions = [MODELLER_DIC['version']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(MODELLER_DIC['home'], MODELLER_DIC['name'] + '-' + MODELLER_DIC['version'])

    @classmethod
    def defineBinaries(cls, env):
        tmpParams= 'installModellerParams.txt'

        installationCmd = 'wget %s -O %s && ' % (cls._getModellerDownloadUrl(), cls._getModellerTar())
        installationCmd += 'tar -xf %s --strip-components 1 && ' % cls._getModellerTar()
        installationCmd += 'rm %s && ' % cls._getModellerTar()

        installationCmd += 'printf "%s" > %s && ' % (FILE_INSTALLER.format(cls._pluginHome), tmpParams)
        installationCmd += './Install < %s > /dev/null 2>&1 && ' % tmpParams
        installationCmd += 'rm %s && ' % tmpParams

        # Creating validation file
        MODELLER_INSTALLED = '%s_installed' % MODELLER_DIC['name']
        installationCmd += 'touch %s' % MODELLER_INSTALLED  # Flag installation finished

        env.addPackage(MODELLER_DIC['name'],
                       version=MODELLER_DIC['version'],
                       tar='void.tgz',
                       commands=[(installationCmd, MODELLER_INSTALLED)],
                       neededProgs=["wget", "tar"],
                       default=True)

    @classmethod
    def runModeller(cls, protocol, program, args, cwd=None):
        """ Run Modeller command from a given protocol. """
        modellerArgs = ['python3', program, *args]
        protocol.runJob(join(cls.getVar(MODELLER_DIC['home']), 'bin/modpy.sh'), modellerArgs, cwd=cwd)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def _getModellerDownloadUrl(cls):
        print(yellowStr('\nOnce Modeller is downloaded and installed, remember to write the license key '
                        'in file {}/modlib/modeller/config.py. This key can be obtained by registration in '
                        'https://salilab.org/modeller/registration.html\n'.
                        format(MODELLER_DIC['home'])))
        return "\'https://salilab.org/modeller/10.1/modeller-10.1.tar.gz\'"

    @classmethod
    def _getModellerTar(cls):
        pluginHome = join(pwem.Config.EM_ROOT, MODELLER_DIC['name'] + '-' + MODELLER_DIC['version'])
        return pluginHome + '/' + MODELLER_DIC['name'] + '-' + MODELLER_DIC['version'] + '.tar.gz'
