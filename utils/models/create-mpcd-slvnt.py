#!/usr/bin/env python

'''
Creates initial configuration for a system containing pply MPCD particles.

'''
import sys
import math
import numpy as np
from configuration import Configuration
from config_io import *

#-------------------------------------------------------------------------------

config = Configuration()

#Add simulation box with PBC.
config.add_simbox(20.0, 20.0, 20.0, 1)

'''
#The system consists of polymer chains and couterions.

#Atom types
#Monomer atoms as point particles with mass 1.0
atm_t_m = config.add_atom_type('M', 1, 1.0)
#Counterion atoms as point particles with mass 1.0.
atm_t_c = config.add_atom_type('C', 1, 1.0)

#Bond types
#Bonds along the chain. (Kremer-Grest)
eps = 1.0; sigma = 1.0; rcut = 2.0**(1.0/6)*sigma
bnd_t_chn = config.add_bond_type(3, np.array([30*eps/sigma**2, 1.5*sigma, eps, sigma]))

#Angle types
#Angles along the chain. (Cosine)
ang_t_chn = config.add_angle_type(1, np.array([2.0, 1.0]))

#Vdw (pairwise) interaction. 12-6 LJ+Coulomb (cut & shifted)
#Assume same for all pairs
eps = 1.0; sigma = 2.0; rcut = 2.0**(1.0/6)*sigma
rcut_coul = rcut + 3.0; C = 20.0
config.add_ia_vdw('M', 'M', 5, np.array([eps, sigma, rcut, rcut_coul, C]))
config.add_ia_vdw('C', 'C', 5, np.array([eps, sigma, rcut, rcut_coul, C]))
config.add_ia_vdw('M', 'C', 5, np.array([eps, sigma, rcut, rcut_coul, C]))

#Molecules
#Chain, let there be 5 such molecules.
num_mols = 5
mol_t_chn = config.add_molecule_type('CHN', num_mols)

#Number of atoms in a chain
natm = 30

#Total number of bonds in a chain
nbnd = natm - 1

#Number of angles in a chain.
nang = natm - 2 if natm > 2 else 0

#Add atoms, bonds, and angles of chains.
len_bond = 2.0  #Approx. equilibrium bond length
sep = 1.95 # For checking overlap, must be <= len_bond
theta = (45.0/180.0)*math.pi #Avg. initial bond angle

#Loop over all chains
iatm = 0; ibnd = 0; iang = 0 #Index of the last atom, bond, etc. added
for imol in range(1, num_mols+1):
    iatm_beg = iatm + 1
    ibnd_beg = ibnd + 1
    iang_beg = iang + 1
    #First atom can be at any location
    iatm = config.append_atom_unbonded( atm_t_m, 0.0, sep=sep)
    #Add second atom bonded to the first as a freely jointed link
    iatm = config.append_atom_bonded(atm_t_m, 0.0, len_bond, 'efjc',
            iatm, sep=sep)
    #Add rest of the atoms as freely rotating links
    for i in range(3,natm+1):
        iatm = config.append_atom_bonded(atm_t_m, 0.0, len_bond, 'efrc',
                iatm, im2=iatm-1, theta=theta, sep=sep)
    #Add bonds
    for i in range(1, nbnd+1):
        j = iatm_beg + i - 1
        ibnd = config.add_bond(bnd_t_chn, j, j+1)
    
    #Add angles
    for i in range(1, nang+1):
        j = iatm_beg + i - 1
        iang = config.add_angle(ang_t_chn, j, j+1, j+2)

    #Add the whole molecule 
    config.add_molecule(mol_t_chn, natm, iatm_beg, nbnd, ibnd_beg, nang, iang_beg, 0, 0)

#Now add counterions. These are monoatomic molecules. Let's keep their name same
#as the atomic species.
#Let there be 50 counterions.
num_mols = 50
mol_t_ci = config.add_molecule_type('C', num_mols)
#No bonds & angles. Let's keep sep same as before. Put these counterions
#randomly within the box.
for imol in range(1, num_mols+1):
    iatm_beg = iatm + 1
    iatm = config.append_atom_unbonded( atm_t_c, 0.0, sep=sep)
    #Add the whole molecule 
    config.add_molecule(mol_t_ci, 1, iatm_beg, 0, 0, 0, 0, 0, 0)

#Randomly assign charge, while keeping the whole system neutral.
config.add_random_charges()

#Apply periodic boundary conditions
config.apply_pbc()
'''

#Add MPCD atoms
config.add_mpcd_atoms(5)

#write_ldf(config, 'test.txt', title='test_frame')
write_cfg(config, 'sample.cfg', title='Sample')
