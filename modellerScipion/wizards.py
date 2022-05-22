# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from .protocols import ModellerMutateResidue
from pwem.wizards import SelectChainWizard, SelectResidueWizard, pwobj, EmWizard

from .constants import AA_LIST

SelectChainWizard().addTarget(protocol=ModellerMutateResidue,
                              targets=['mutChain'],
                              inputs=['inputAtomStruct'],
                              outputs=['mutChain'])

SelectResidueWizard().addTarget(protocol=ModellerMutateResidue,
                                targets=['mutPosition'],
                                inputs=['inputAtomStruct', 'mutChain'],
                                outputs=['mutPosition'])

class AddMutationWizard(EmWizard):
  _targets = [(ModellerMutateResidue, ['addMutation'])]

  def show(self, form, *params):
    protocol = form.protocol
    chain, pos, resId = protocol.mutChain.get(), protocol.mutPosition.get(), protocol.mutResidue.get()
    res = AA_LIST[resId]
    form.setVar('toMutateList', protocol.toMutateList.get() +
                '{} | {} | {}\n'.format(chain, pos, res))
    
class ClearMutationList(EmWizard):
  _targets = [(ModellerMutateResidue, ['clearLabel'])]

  def show(self, form, *params):
    form.setVar('toMutateList', '')
