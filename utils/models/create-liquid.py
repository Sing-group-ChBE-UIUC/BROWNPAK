#!/usr/bin/env python

'''
Creates initial configuration for a system containing unbranched linear chains.

'''
import sys
import math
import numpy as np
from configuration import Configuration
from config_io import *

#-------------------------------------------------------------------------------

config = Configuration()

#Add simulation box and boundary condition. The box size can be updated later if
#necessary.
config.add_simbox(10.0, 10.0, 10.0, 1)

#Molecule details
num_mols_A = 60; num_mols_B = 60
mol_t_A = config.add_molecule_type('A', num_mols_A)
mol_t_B = config.add_molecule_type('B', num_mols_B)

#Total number of atoms in the molecule
natm = 1

#Atom types
# natm atoms as point particles with mass 1.0
atm_t_a = config.add_atom_type('A', 1, 1.0)
atm_t_b = config.add_atom_type('B', 1, 1.0)

#Vdw interaction
eps = 1.0; sigma = 2.0; rcut = 2.0**(1.0/6)*sigma
rcut_coul = 0; C = 2.0
config.add_ia_vdw('A', 'A', 5, np.array([eps, sigma, rcut, rcut_coul, C]))
config.add_ia_vdw('B', 'B', 5, np.array([eps, sigma, rcut, rcut_coul, C]))
config.add_ia_vdw('A', 'B', 5, np.array([eps, sigma, rcut, rcut_coul, C]))

#Add atoms, bonds, and angles
sep = 1.95 # Must be <= len_bond

#Loop over all molecules
imol = 1; iatm = 1; ibnd = 1; iang = 1 #Index of one past the last atom, bond, etc. added
for imol in range(1, num_mols_A+1):
    iatm_beg = iatm
    ibnd_beg = ibnd
    iang_beg = iang

    config.append_atom_unbonded( atm_t_a, 1.0, sep=sep )
    iatm += 1
    
    config.add_molecule(mol_t_A, natm, iatm_beg, 0, ibnd_beg, 0, iang_beg, 0, 0)

for imol in range(1, num_mols_B+1):
    iatm_beg = iatm
    ibnd_beg = ibnd
    iang_beg = iang

    config.append_atom_unbonded( atm_t_b, -1.0, sep=sep )
    iatm += 1
    
    config.add_molecule(mol_t_B, natm, iatm_beg, 0, ibnd_beg, 0, iang_beg, 0, 0)

config.apply_pbc()

write_ldf(config, 'tfr.txt', title='test_frame')
write_cfg(config, 'colloid.cfg', title='test_cfg')
