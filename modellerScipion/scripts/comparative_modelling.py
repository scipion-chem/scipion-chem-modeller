# A sample script for fully automated comparative modeling
# https://salilab.org/modeller/manual/node32.html

import os, argparse
from modeller import *
from modeller.automodel import *  # Load the AutoModel class
from modeller.parallel import *

def special_restraints(self, aln):
    # Constrain the A and B chains to be identical (but only restrain
    # the C-alpha atoms, to reduce the number of interatomic distances
    # that need to be calculated):
    if hasattr(self, 'symChains'):
        for chainPair in self.symChains:
            s1 = Selection(self.chains[chainPair[0]]).only_atom_types(self.symAtom)
            s2 = Selection(self.chains[chainPair[1]]).only_atom_types(self.symAtom)
            self.restraints.symmetry.append(Symmetry(s1, s2, 1.0))

def special_patches(self, aln):
    if hasattr(self, 'renam'):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=self.renam,
                             renumber_residues=self.renum)

def parsePDBCodes(pdbsFile):
    codes = []
    with open(pdbsFile) as f:
        for line in f:
            codes.append(line.strip())
    return codes

def parseScore(scoreStr):
    scores = scoreStr.split(',')
    ase, scoNames = [], []
    if 'DOPE' in scores:
        ase.append(assess.DOPE), scoNames.append('DOPE score')
    if 'DOPE-HR' in scores:
        ase.append(assess.DOPEHR), scoNames.append('DOPE-HR score')
    if 'Normalized_DOPE' in scores:
        ase.append(assess.normalized_dope), scoNames.append('Normalized DOPE score')
    if 'GA341' in scores:
        ase.append(assess.GA341), scoNames.append('GA341 score')

    return ase, scoNames

def parseSymmetries(symStr):
    chainPairs = []
    for cPair in symStr.split(','):
        chainPairs.append(tuple(cPair.strip().split('-')))
    return chainPairs

def comparativeModelling():
    setattr(AutoModel, 'special_restraints', special_restraints)
    setattr(AutoModel, 'special_patches', special_patches)
    setattr(AllHModel, 'special_restraints', special_restraints)
    setattr(AllHModel, 'special_patches', special_patches)

    parser = argparse.ArgumentParser(description='Mutate residue from a given chain of a pdb file')
    parser.add_argument('-i', '--inputSeqName', type=str, help='Name of the sequence to model in the alignment')
    parser.add_argument('-af', '--alignFile', type=str, help='Input alignment PIR file')
    parser.add_argument('-pf', '--pdbsFile', type=str, help='File containing the template PDBs specifications')
    parser.add_argument('-pd', '--pdbsDir', type=str, help='Directory containing atomic structure files')
    parser.add_argument('--align', default=False, action='store_true', help='Automatic align of the sequences')
    parser.add_argument('-n', '--nModels', type=int, default=1, required=False, help='Number of models')

    parser.add_argument('-im', '--iniModel', type=str, default='', help='File containing the initial PDB model')
    parser.add_argument('--modelH', default=False, action='store_true', help='Optimize also hydrogens')
    parser.add_argument('-sc', '--score', type=str, default='', help='Score of the finals models to save')
    parser.add_argument('-opt', '--optimization', type=str, default='', help='Quality of the optimization')
    parser.add_argument('-nr', '--nReps', type=int, default=1, required=False,
                        help='Number of optimization repetitions')
    parser.add_argument('-renum', type=str, default='', help='Renumber each chain first residue index')
    parser.add_argument('-renam', type=str, default='', help='Rename each chain name')
    parser.add_argument('-sym', '--symmetry', type=str, default='', help='Symmetry restrains by chain')
    parser.add_argument('-symAtom', '--symmetryAtom', type=str, default='',
                        help='Type of atoms to check the symmetry on')
    parser.add_argument('-nj', '--nCPUs', type=int, default=1, required=False,
                        help='Number of CPUs')
    parser.add_argument('-mPath', '--modellerPath', type=str, default='',
                        help='Path to modeller home')

    args = parser.parse_args()

    targetName, alignFile = args.inputSeqName, args.alignFile
    pdbCodes, pdbDir = parsePDBCodes(args.pdbsFile), args.pdbsDir
    align, modelH = args.align, args.modelH

    nModels, nReps = args.nModels, args.nReps
    score, optim = args.score, args.optimization

    if score != '':
        scoreFuncs, scoreKeys = parseScore(score)
    else:
        scoreFuncs, scoreKeys = None, ''

    iniModel = args.iniModel
    if iniModel == '':
        iniModel = None

    log.verbose()
    env = Environ()

    env.io.atom_files_directory = ['.', pdbDir]

    function = AutoModel if not modelH else AllHModel

    scF = tuple(scoreFuncs) if scoreFuncs else None
    a = function(env, alnfile=alignFile,
                 knowns=tuple(pdbCodes), sequence=targetName,
                 assess_methods=scF,
                 inifile=iniModel)

    ncpus, modellerPath = args.nCPUs, args.modellerPath
    if ncpus > 1:
        j = job()
        for i in range(ncpus):
            j.append(LocalWorker())
        a.use_parallel_job(j)

    if args.renam:
        nums = list(map(int, args.renum.split(','))) if args.renum else []
        names = args.renam.split(',')
        a.renam, a.renum = names, nums

    if args.symmetry != '':
        a.symChains = parseSymmetries(args.symmetry)
        a.symAtom = args.symmetryAtom

    a.starting_model = 1
    a.ending_model = nModels

    if align:
        a.auto_align()  # get an automatic alignment

    if optim == 'Low-Fast':
        a.very_fast()
    elif optim == 'High-Slow':
        a.library_schedule = autosched.slow
        a.max_var_iterations = 300
        a.md_level = refine.slow


    a.repeat_optimization = nReps
    a.make()  # do comparative modeling

    # Get a list of all successfully built models from a.outputs
    ok_models = [x for x in a.outputs if x['failure'] is None]

    # Rank the models by DOPE score
    if len(scoreKeys) > 0:
        scoreStr = ''
        for m in ok_models:
            outId = int(os.path.splitext(m['name'])[0][-2:])
            scoreStr += "Model: 'outputAtomStruct_{}".format(outId)
            for scoreKey in scoreKeys:
                if type(m[scoreKey]) == list:
                    m[scoreKey] = m[scoreKey][0]
                scoreStr += "(%s %.3f)" % (scoreKey, m[scoreKey])
            scoreStr += '\n'

        with open('scores.txt', 'w') as f:
            f.write(scoreStr)

if __name__ == '__main__':
    comparativeModelling()