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


def comparativeModelling():
    parser = argparse.ArgumentParser(description='Mutate residue from a given chain of a pdb file')
    parser.add_argument('-i', '--inputSeqName', type=str, help='Name of the sequence to model in the alignment')
    parser.add_argument('-af', '--alignFile', type=str, help='Input alignment PIR file')
    parser.add_argument('-pf', '--pdbsFile', type=str, help='File containing the template PDBs specifications')
    parser.add_argument('-pd', '--pdbsDir', type=str, help='Directory containing atomic structure files')
    parser.add_argument('--align', default=False, action='store_true', help='Automatic align of the sequences')
    parser.add_argument('-n', '--nModels', type=int, default=1, required=False, help='Number of models')

    args = parser.parse_args()

    targetName = args.inputSeqName
    alignFile = args.alignFile
    pdbCodes = parsePDBCodes(args.pdbsFile)
    pdbDir = args.pdbsDir
    align = args.align
    nModels = args.nModels

    log.verbose()
    env = Environ()

    env.io.atom_files_directory = ['.', pdbDir]

    a = AutoModel(env,
                  # file with template codes and target sequence
                  alnfile=alignFile,
                  # PDB codes of the templates
                  knowns=tuple(pdbCodes),
                  # code of the target
                  sequence=targetName,
                  assess_methods=(assess.DOPE))

    a.starting_model = 1
    a.ending_model = nModels

    if align:
        a.auto_align()  # get an automatic alignment
    a.make()  # do comparative modeling

    # Get a list of all successfully built models from a.outputs
    ok_models = [x for x in a.outputs if x['failure'] is None]

    # Rank the models by DOPE score
    key = 'DOPE score'
    scoreStr = ''
    for m in ok_models:
        outId = int(os.path.splitext(m['name'])[0][-2:])
        scoreStr += "Model: %s (DOPE score %.3f)\n" % ('outputAtomStruct_{}'.format(outId), m[key])

    with open('scores.txt', 'w') as f:
        f.write(scoreStr)

if __name__ == '__main__':
    comparativeModelling()