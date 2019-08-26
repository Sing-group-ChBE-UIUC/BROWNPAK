#!/usr/bin/env python

'''
Creates initial configuration for a system containing single unbranched linear
chain in an unbounded domain.

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
config.add_simbox(20.0, 20.0, 20.0, 0)

#Molecule details
num_mols = 1
config.add_molecule_type('CHN-UB', num_mols)

#Total number of atoms in the molecule
natm = 30 #int(sys.argv[1])

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

#Tether types (No tethers)
#teth_t = config.add_tether_type(1, np.array([2.0, 1e-6]))

#Add atoms, bonds, and angles
len_bond = 2.0 #0.97*sigma #Equilibrium bond length
theta = (25.0/180.0)*math.pi
sep = 2.0 # Must be <= len_bond

#Loop over all molecules (here only a single molecule type)
imol = 1; iatm = 1; ibnd = 1; iang = 1 #Index of one past the last atom, bond, etc. added
for imol in range(1, num_mols+1):
    iatm_beg = iatm
    ibnd_beg = ibnd
    iang_beg = iang

    config.append_atom_unbonded( atm_t_bb, 0.0 )
    iatm += 1
    config.append_atom_bonded(1, 0.0, len_bond, 'efjc', iatm-1, sep=sep)
    iatm += 1
    for i in range(3,natm+1):
        config.append_atom_bonded(atm_t_bb, 0.0, len_bond, 'efjc',
                iatm-1, im2=iatm-2, theta=theta, sep=sep)
        iatm += 1
    
    for i in range(1, nbnd+1):
        j = iatm_beg + i - 1
        config.add_bond(bnd_t_bb, j, j+1)
        ibnd += 1
    
    for i in range(1, nang+1):
        j = iatm_beg + i - 1
        config.add_angle(ang_t_bb, j, j+1, j+2)
        iang += 1
    
    config.add_molecule(1, natm, iatm_beg, nbnd, ibnd_beg, nang, iang_beg, 0, 0)

#config.add_tether(teth_t, 1, np.zeros((3,)) )

config.update_simbox()
config.to_center(1)

config.add_flow_field( 1, np.array([0.1]) )

write_ldf(config, 'tfr.txt', title='test_frame')
write_cfg(config, 'lin-%d.cfg'%natm, title='test_cfg')
