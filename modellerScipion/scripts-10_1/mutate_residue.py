import os, argparse

from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched

#
#  mutate_model.py
#
#     Usage:   python mutate_model.py modelname respos resname chain > logfile
#
#     Example: python mutate_model.py 1t29 1699 LEU A > 1t29.log
#
#
#  Creates a single in silico point mutation to sidechain type and at residue position
#  input by the user, in the structure whose file is modelname.pdb
#  The conformation of the mutant sidechain is optimized by conjugate gradient and
#  refined using some MD.
#
#  Note: if the model has no chain identifier, specify "" for the chain argument.
#  Adaptation of modeller example by Scipion team
#


def optimize(atmsel, sched):
    #conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


#molecular dynamics
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False


#use homologs and dihedral library for dihedral angle restraints
def make_restraints(mdl1, aln):
   rsr = mdl1.restraints
   rsr.clear()
   s = Selection(mdl1)
   for typ in ('stereo', 'phi-psi_binormal'):
       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)

def mutateResidue():
    parser = argparse.ArgumentParser(description='Mutate residue from a given chain of a pdb file')
    parser.add_argument('-i', '--inputFilename', type=str, help='Input pdb file')
    parser.add_argument('-p', '--positions', type=str, help='Residue position to mutate')
    parser.add_argument('-r', '--newResidues',type=str, help='Residue to add in the substitution')
    parser.add_argument('-c', '--chains',type=str, help='Chain of the protein to mutate')
    parser.add_argument('-s', '--seed',type=int, default=-49837, required=False, help='Random seed')
    parser.add_argument('-o', '--outputFile',type=str, help='Output file')

    args = parser.parse_args()
    modelname, chains, respos, restypes = args.inputFilename, \
                                          eval(args.chains), eval(args.positions), eval(args.newResidues)
    seed, outputFile = args.seed, args.outputFile

    log.verbose()

    # Set a different value for rand_seed to get a different final model
    env = Environ(rand_seed=seed)

    env.io.hetatm = True
    #soft sphere potential
    env.edat.dynamic_sphere=False
    #lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')


    # Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = Model(env, file=modelname)
    ali = Alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

    for mutIdx in range(len(respos)):
        chain, resp, restyp = chains[mutIdx], respos[mutIdx], restypes[mutIdx]
        print(chain, resp, restyp)
        #set up the mutate residue selection segment
        s = Selection(mdl1.chains[chain].residues[resp])
        print(1)

        #perform the mutate residue operation
        s.mutate(residue_type=restyp)
        #get two copies of the sequence.  A modeller trick to get things set up
        ali.append_model(mdl1, align_codes=modelname)
        print(2)
        # Generate molecular topology for mutant
        mdl1.clear_topology()
        mdl1.generate_topology(ali[-1])

        # Transfer all the coordinates you can from the template native structure
        # to the mutant (this works even if the order of atoms in the native PDB
        # file is not standard):
        #here we are generating the model by reading the template coordinates
        mdl1.transfer_xyz(ali)
        print(3)
        # Build the remaining unknown coordinates
        mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

        #yes model2 is the same file as model1.  It's a modeller trick.
        mdl2 = Model(env, file=modelname)
        print(4)
        #required to do a transfer_res_numb
        #ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
        #transfers from "model 2" to "model 1"
        mdl1.res_num_from(mdl2,ali)

        #It is usually necessary to write the mutated sequence out and read it in
        #before proceeding, because not all sequence related information about MODEL
        #is changed by this command (e.g., internal coordinates, charges, and atom
        #types and radii are not updated).

        mdl1.write(file=modelname+restyp+resp+'.tmp')
        mdl1.read(file=modelname+restyp+resp+'.tmp')

        #set up restraints before computing energy
        #we do this a second time because the model has been written out and read in,
        #clearing the previously set restraints
        make_restraints(mdl1, ali)

        #a non-bonded pair has to have at least as many selected atoms
        mdl1.env.edat.nonbonded_sel_atoms=1

        sched = autosched.loop.make_for_model(mdl1)

        #only optimize the selected residue (in first pass, just atoms in selected
        #residue, in second pass, include nonbonded neighboring atoms)
        #set up the mutate residue selection segment
        s = Selection(mdl1.chains[chain].residues[resp])

        mdl1.restraints.unpick_all()
        mdl1.restraints.pick(s)

        s.energy()

        s.randomize_xyz(deviation=4.0)

        mdl1.env.edat.nonbonded_sel_atoms=2
        optimize(s, sched)

        #feels environment (energy computed on pairs that have at least one member
        #in the selected)
        mdl1.env.edat.nonbonded_sel_atoms=1
        optimize(s, sched)

        s.energy()

        # delete the temporary file
        os.remove(modelname + restyp + resp + '.tmp')

    #give a proper name
    tmpFile = outputFile.replace('.pdb', '_tmp.pdb')
    mdl1.write(file=tmpFile)
    fTmp = open(tmpFile, 'r').read()
    with open(outputFile, 'w') as f:
        f.write('#Mutations generated with modeller: \n')
        for i in range(len(respos)):
            f.write('#Chain: {} | Position: {} | New residue: {}\n'.format(chains[i], respos[i], restypes[i]))
        f.write(fTmp)
    os.remove(tmpFile)

if __name__ == '__main__':
    mutateResidue()