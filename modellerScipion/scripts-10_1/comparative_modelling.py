# A sample script for fully automated comparative modeling
# https://salilab.org/modeller/manual/node32.html

import os, argparse
from modeller import *
from modeller.automodel import *  # Load the AutoModel class

def parsePDBCodes(pdbsFile):
    codes = []
    with open(pdbsFile) as f:
        for line in f:
            codes.append(line.strip())
    return codes

def parseScore(scoreStr):
    if scoreStr == 'molpdf':
        return None, scoreStr
    elif scoreStr == 'DOPE':
        return assess.DOPE, scoreStr + ' score'
    elif scoreStr == 'DOPE-HR':
        return assess.DOPEHR, scoreStr + ' score'
    elif scoreStr == 'Normalized_DOPE':
        return assess.normalized_dope, scoreStr.replace('_', ' ') + ' score'
    elif scoreStr == 'GA341':
        return assess.GA341, scoreStr + ' score'

def comparativeModelling():
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

    args = parser.parse_args()

    targetName, alignFile = args.inputSeqName, args.alignFile
    pdbCodes, pdbDir = parsePDBCodes(args.pdbsFile), args.pdbsDir
    align, modelH = args.align, args.modelH

    nModels, nReps = args.nModels, args.nReps
    score, optim = args.score, args.optimization
    if score != '':
        scoreFunc, scoreKey = parseScore(score)
    else:
        scoreFunc, scoreKey = None, ''

    iniModel = args.iniModel
    if iniModel == '':
        iniModel = None

    log.verbose()
    env = Environ()

    env.io.atom_files_directory = ['.', pdbDir]

    function = 'AutoModel' if not modelH else 'AllHModel'

    a = eval(function)(env, alnfile=alignFile,
                      knowns=tuple(pdbCodes), sequence=targetName,
                      assess_methods=(scoreFunc),
                      inifile=iniModel)

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
    if scoreKey:
        scoreStr = ''
        for m in ok_models:
            outId = int(os.path.splitext(m['name'])[0][-2:])
            if type(m[scoreKey]) == list:
                m[scoreKey] = m[scoreKey][0]
            scoreStr += "Model: %s (%s %.3f)\n" % ('outputAtomStruct_{}'.format(outId), scoreKey, m[scoreKey])

        with open('scores.txt', 'w') as f:
            f.write(scoreStr)

if __name__ == '__main__':
    comparativeModelling()