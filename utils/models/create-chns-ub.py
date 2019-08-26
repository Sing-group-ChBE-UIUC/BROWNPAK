#!/usr/bin/env python

'''
Creates initial configuration for a system containing unbranched linear
chains in a triply periodic domain.

'''
import sys
import math
import numpy as np
from configuration import Configuration
from config_io import *

#-------------------------------------------------------------------------------

config = Configuration()

#Add simulation box and boundary condition. The box size is fixed.
config.add_simbox(10.0, 10.0, 10.0, 1)

#Molecule details
num_mols = 10
mol_t = config.add_molecule_type('CHN-UB', num_mols)

#Total number of atoms in the molecule
natm = 10 #int(sys.argv[1])

#Total number of bonds
nbnd = natm - 1

#Total number of angles. If no angles set num_angles to zero.
nang = natm - 2 if natm > 2 else 0

#Atom types
# natm atoms as point particles with mass 1.0
atm_t_bb = config.add_atom_type('M', 1, 1.0)

#Vdw interaction
eps = 1.0; sigma = 2.0; rcut = 2.0**(1.0/6)*sigma
config.add_ia_vdw('M', 'M', 1, np.array([eps, sigma, rcut]))

#Bond types
bnd_t_bb = config.add_bond_type(3, np.array([30*eps/sigma**2, 1.5*sigma, eps, sigma]))

#Angle types
ang_t_bb = config.add_angle_type(1, np.array([2.0]))

#Add atoms, bonds, and angles
len_bond = 0.97*sigma #Equilibrium bond length
theta = (25.0/180.0)*math.pi
sep = 1.9 # Must be <= len_bond

#Loop over all molecules (here only a single molecule type)
imol = 1; iatm = 1; ibnd = 1; iang = 1 #Index of one past the last atom, bond, etc. added
for imol in range(1, num_mols+1):
    iatm_beg = iatm
    ibnd_beg = ibnd
    iang_beg = iang

    #Put the first atom anywhere within the simulation box
    config.append_atom_unbonded( atm_t_bb, 0.0 )
    iatm += 1
    #Successively add other atoms
    config.append_atom_bonded(1, 0.0, len_bond, 'efjc', iatm-1, sep=sep)
    iatm += 1
    for i in range(3,natm+1):
        config.append_atom_bonded(atm_t_bb, 0.0, len_bond, 'efjc',
                iatm-1, im2=iatm-2, theta=theta, sep=sep)
        iatm += 1
    
    #Specify bonds
    for i in range(1, nbnd+1):
        j = iatm_beg + i - 1
        config.add_bond(bnd_t_bb, j, j+1)
        ibnd += 1
    
    #Specify angles
    for i in range(1, nang+1):
        j = iatm_beg + i - 1
        config.add_angle(ang_t_bb, j, j+1, j+2)
        iang += 1

    #Specify the whole molecule 
    config.add_molecule(1, natm, iatm_beg, nbnd, ibnd_beg, nang, iang_beg, 0, 0)


config.apply_pbc()

write_ldf(config, 'tfr.txt', title='test_frame')
write_cfg(config, 'chns-%d.cfg'%natm, title='test_cfg')
