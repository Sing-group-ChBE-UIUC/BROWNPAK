#!/usr/bin/env python

'''
Creates initial configuration for a system containing single branched linear 
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
config.add_simbox(40.0, 40.0, 40.0, 0)

#Number of backbone atoms
na_bbone = 30 #int(sys.argv[1])
#Number of atoms on each side chain
na_sc = 4 #int(sys.argv[2])
#Number of atoms between consecutive branch points
na_sp = 0
#Number of side chains growing from a branch point
f = 1

#Atom types
#All atoms are point particles with unit mass
atm_t_bb = config.add_atom_type('BB', 1, 1.0)
atm_t_sc = config.add_atom_type('SC', 1, 1.0)

#Vdw interaction
eps = 1.0; sigma = 2.0; rcut = 2.0**(1.0/6)*sigma
config.add_ia_vdw('BB', 'BB', 1, np.array([eps, sigma, rcut]))
config.add_ia_vdw('BB', 'SC', 1, np.array([eps, sigma, rcut]))
config.add_ia_vdw('SC', 'SC', 1, np.array([eps, sigma, rcut]))

#Bond types
bnd_t_bb = config.add_bond_type(3, np.array([30*eps/sigma**2, 1.5*sigma, eps, sigma]))
bnd_t_bs = config.add_bond_type(3, np.array([30*eps/sigma**2, 1.5*sigma, eps, sigma]))
bnd_t_ss = config.add_bond_type(3, np.array([30*eps/sigma**2, 1.5*sigma, eps, sigma]))

#Angle types (No angles)

#Dihedral types (No dihedrals)

#Tether types (No tethers)

#Molecule type
num_mols = 1
config.add_molecule_type('BTLBRS', num_mols)

len_bond = 1.94 #0.97*sigma #Equilibrium bond length
theta_bb = (45.0/180.0)*math.pi
sep = 1.9 # Must be <= len_bond

#Add atoms, bonds, etc
#Loop over all molecules (here only a single molecule type)
imol = 1; iatm = 1; ibnd = 1 #Index of atom, bond, etc. of the next entity to be added
for imol in range(1, num_mols+1):
    iatm_beg = iatm
    ibnd_beg = ibnd

    #Backbone
    #Put first atom at the origin
    config.add_atom(atm_t_bb, 0.0, np.zeros((3,)))
    iatm += 1
    #Put second atom along x-axis
    config.append_atom_bonded(atm_t_bb, 0.0, len_bond, 'alignx', iatm-1)
    iatm += 1
    #Put the rest of the backbone atoms
    for i in range(3,na_bbone+1):
        config.append_atom_bonded(atm_t_bb, 0.0, len_bond, 'alignx',
                iatm-1, im2=iatm-2, theta=theta_bb, sep=sep)
        iatm += 1
    #Backbone bonds 
    nbnd = na_bbone - 1
    for i in range(1, nbnd+1):
        config.add_bond(bnd_t_bb, i, i+1)
        ibnd += 1

    #Set backbone branch
    config.set_branch(-1, 1, na_bbone)

    #Side chains (arranged as efjc)
    for i in range(1, na_bbone+1, na_sp+1):
        ibp = i #Branch point
        for i_sc in range(1, f+1):
            for ia_sc in range(1, na_sc+1):
                if ia_sc == 1:
                    config.append_atom_bonded(atm_t_sc, 0.0, len_bond, 'efjc',
                            ibp, sep=sep)
                    iatm += 1
                    config.add_bond(bnd_t_bs, ibp, iatm-1)
                    ibnd +=1
                    ia_br_beg = iatm - 1
                else:
                    config.append_atom_bonded(atm_t_sc, 0.0, len_bond, 'efjc',
                            iatm-1, sep=sep)
                    iatm += 1
                    config.add_bond(bnd_t_ss, iatm-2, iatm-1)
                    ibnd +=1
            #Set branch
            config.set_branch(ibp, ia_br_beg, na_sc)

    #Add molecule
    config.add_molecule(1, iatm-1, iatm_beg, ibnd-1, ibnd_beg, 0, 0, 0, 0)

config.update_simbox()
config.to_center(1)
config.add_flow_field( 0, np.array([]) )

write_ldf(config, 'tfr.txt', title='test_brush')
write_cfg(config, 'bb-sus-%d-%d.cfg'%(na_bbone, na_sc), title='brush')
