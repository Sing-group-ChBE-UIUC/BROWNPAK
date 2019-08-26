#!/usr/bin/env python

import sys
import warnings
import math
import numpy as np
from geom_utils import rotate_vectors_random, get_ransphere, get_frc_link


class Configuration(object):

    def __init__(self):
        self.simbox = np.zeros((3,))
        self.imcon = 0

        self.atom_types = {}
        self.bond_types = {}
        self.angle_types = {}
        self.dihedral_types = {}
        self.tether_types = {}
        self.molecule_types = {}
        self.vdw = {}
        self.external = {}
        self.flow_field = None

        self.atoms = {}
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.molecules = {}
        self.tethers = {}
        self.branches = {}
        self.mpcd_avnc = 0
        self.num_mpcd_atoms = 0
        

    def add_simbox(self, lx, ly, lz, imcon):
        '''
        Adds simulation box with dimensions lx x ly x lz. The box is placed in
        the first octant with its left, bottom, back corner at the origin. x-dir
        is toward the right, y is up, and z points out of the screen.
        imcon = 0: Unbounded domain
        imcon = 1: PBC along all directions.

        '''
        assert imcon in [0, 1]
        self.simbox = np.array([lx, ly, lz])
        self.imcon = imcon


    def add_atom_type(self, name, style, mass):
        '''
        Adds a new atom type. The type defines an atom with name `name`, style `style`,
        and mass `mass`.

        Returns an integer identifier for the atom type added.

        Style 1: Point particle (no orientation)
        Style 2: Body particle (orientation must be specified)

        '''
        if len(name) > 8:
            print('Atom name `%s` longer than 8 characters.'%name)
            print('Truncating name to the first 8 characters.')
            atm_name = name[0:8]
        else:
            atm_name = name
        assert mass > 0
        iat = len(self.atom_types) + 1
        self.atom_types[iat] = {'name': atm_name, 'style': style, 'mass': mass}
        return iat


    def add_bond_type(self, style, params):
        ibt = len(self.bond_types) + 1
        self.bond_types[ibt] = {'style': style, 'params': np.copy(params)}
        return ibt


    def add_angle_type(self, style, params):
        iant = len(self.angle_types) + 1
        self.angle_types[iant] = {'style': style, 'params': np.copy(params)}
        return iant


    def add_dihedral_type(self, style, params):
        idt = len(self.dihedral_types) + 1
        self.dihedral_types[idt] = {'style': style, 'params': np.copy(params)}
        return idt


    def add_tether_type(self, style, params):
        itt = len(self.tether_types) + 1
        self.tether_types[itt] = {'style': style, 'params': np.copy(params)}
        return itt


    def add_molecule_type(self, name, nmol):
        if len(name) > 8:
            print('Molecule name `%s` longer than 8 characters.'%name)
            print('Truncating name to the first 8 characters.')
            mol_name = name[0:8]
        else:
            mol_name = name
        assert nmol > 0
        imt = len(self.molecule_types) + 1
        self.molecule_types[imt] = {'name': mol_name, 'nmol': nmol}
        return imt


    def add_atom(self, typ, charge, coords, orientation=None):
        assert typ in self.atom_types
        atm_id = len(self.atoms) + 1
        self.atoms[atm_id] = {'type': typ, 'charge': charge,
                'coords': np.copy(coords[0:3])}
        if orientation:
            self.atoms[atm_id]['orientation'] = np.copy(orientation[0:4])
        return atm_id


    def append_atom_bonded(self, typ, charge, len_bond, method, im1,
            im2=None, theta=None, uhat=None, sep=None, orientation=None):
        assert typ in self.atom_types
        atm_id = len(self.atoms) + 1

        itr = 0; maxitr = 100 #For method = 'efjc'|'efrc'

        if method == 'alignx':
            rim1 = self.atoms[im1]['coords']
            ri = rim1 + np.array([len_bond, 0.0, 0.0])
        elif method == 'aligny':
            rim1 = self.atoms[im1]['coords']
            ri = rim1 + np.array([0.0, len_bond, 0.0])
        elif method == 'alignz':
            rim1 = self.atoms[im1]['coords']
            ri = rim1 + np.array([0.0, 0.0, len_bond])
        elif method == 'align':
            rim1 = self.atoms[im1]['coords']
            ri = rim1 + len_bond*uhat
        elif method == 'fjc':
            rim1 = self.atoms[im1]['coords']
            ri = rim1 + get_ransphere(len_bond)
        elif method == 'frc':
            rim2 = self.atoms[im2]['coords']
            rim1 = self.atoms[im1]['coords']
            ri = rim1 + len_bond*get_frc_link(rim1-rim2, theta)
        elif method == 'efjc':
            rim1 = self.atoms[im1]['coords']
            while itr < maxitr:
                ri = rim1 + get_ransphere(len_bond)
                if not self.has_overlap(ri, sep):
                    break
                itr += 1
            else:
                warnings.warn('Maximum iteration reached for appending bonded atom')
        elif method == 'efrc':
            rim2 = self.atoms[im2]['coords']
            rim1 = self.atoms[im1]['coords']
            while itr < maxitr:
                ri = rim1 + len_bond*get_frc_link(rim1-rim2, theta)
                if not self.has_overlap(ri, sep):
                    break
                itr += 1
            else:
                warnings.warn('Maximum iteration reached for appending bonded atom')
        else:
            raise ValueError('Unknown method "%s"'%method)

        self.atoms[atm_id] = {'type': typ, 'charge': charge, 'coords': ri}
        if orientation:
            self.atoms[atm_id]['orientation'] = np.copy(orientation[0:4])
        return atm_id


    def append_atom_unbonded(self, typ, charge, sep=None, orientation=None):
        assert typ in self.atom_types
        atm_id = len(self.atoms) + 1
        itr = 0; maxitr = 100
        if not sep:
            ri = np.random.random_sample((3,))
            ri[0] *= self.simbox[0]
            ri[1] *= self.simbox[1]
            ri[2] *= self.simbox[2]
        else:
            while itr < maxitr:
                #Choose any random position within the box. 
                ri = np.random.random_sample((3,))
                ri[0] *= self.simbox[0]
                ri[1] *= self.simbox[1]
                ri[2] *= self.simbox[2]
                if not self.has_overlap(ri, sep):
                    break
                itr += 1
            else:
                warnings.warn('Maximum iteration reached for appending unbonded atom')
        self.atoms[atm_id] = {'type': typ, 'charge': charge, 'coords': ri}
        if orientation:
            self.atoms[atm_id]['orientation'] = np.copy(orientation[0:4])
        return atm_id


    def add_bond(self, typ, atom_i, atom_j):
        assert typ in self.bond_types
        assert atom_i in self.atoms
        assert atom_j in self.atoms
        bnd_id = len(self.bonds) + 1
        self.bonds[bnd_id] = {'type': typ, 'atm_i': atom_i, 'atm_j': atom_j}
        return bnd_id


    def add_angle(self, typ, atom_i, atom_j, atom_k):
        assert typ in self.angle_types
        assert atom_i in self.atoms
        assert atom_j in self.atoms
        assert atom_k in self.atoms
        ang_id = len(self.angles) + 1
        self.angles[ang_id] = {'type': typ, 'atm_i': atom_i, 'atm_j': atom_j,
                'atm_k': atom_k}
        return ang_id


    def add_dihedral(self, typ, atom_i, atom_j, atom_k, atom_l):
        assert typ in self.dihedral_types
        assert atom_i in self.atoms
        assert atom_j in self.atoms
        assert atom_k in self.atoms
        assert atom_l in self.atoms
        dhd_id = len(self.dihedrals) + 1
        self.dihedrals[dhd_id] = {'type': typ, 'atm_i': atom_i, 'atm_j': atom_j,
                'atm_k': atom_k, 'atm_l': atom_l}
        return dhd_id


    def add_molecule(self, typ, natm, iatm_beg, nbnd, ibnd_beg, nang, iang_beg,
            ndhd, idhd_beg):
        assert typ in self.molecule_types
        assert iatm_beg in self.atoms
        assert (iatm_beg+natm-1) in self.atoms
        if nbnd > 0:
            assert ibnd_beg in self.bonds
            assert (ibnd_beg+nbnd-1) in self.bonds
        if nang > 0:
            assert iang_beg in self.angles
            assert (iang_beg+nang-1) in self.angles
        if ndhd > 0:
            assert idhd_beg in self.dihedrals
            assert (idhd_beg+ndhd-1) in self.dihedrals

        imol = len(self.molecules) + 1
        self.molecules[imol] = {'type': typ,
                'natm': natm, 'iatm_beg': iatm_beg,
                'nbnd': nbnd, 'ibnd_beg': ibnd_beg,
                'nang': nang, 'iang_beg': iang_beg,
                'ndhd': ndhd, 'idhd_beg': idhd_beg }
        return imol


    def add_tether(self, typ, iatm, tp):
        '''
        Before calling this, all atoms must be added.

        '''
        assert typ in self.tether_types
        assert iatm in self.atoms
        iteth = len(self.tethers) + 1
        self.tethers[iteth] = {'type': typ, 'iatm': iatm, 'tp': np.copy(tp)}
        return iteth


    def set_branch(self, iatm_teth, iatm_beg, na_br):
        '''
        iatm_teth: index of the atom to which the branch is tethered to.
        iatm_beg: index of the first atom of this branch
        na_br: Number of atoms in this branch
        Branch with id 1 denotes the main branch. For the main branch, iatm_teth is -1.
        '''
        assert iatm_beg in self.atoms
        #Index of the last atom in this branch
        iatm_end = iatm_beg + na_br - 1
        assert iatm_end in self.atoms
        #Number of existing branches
        br_id = len(self.branches) + 1
        if br_id == 1:
            self.branches[br_id] = {'iatm_teth': -1, 'na_br': na_br,
                    'iatm_beg': iatm_beg}
        else:
            assert iatm_teth in self.atoms
            self.branches[br_id] = {'iatm_teth': iatm_teth, 'na_br': na_br,
                    'iatm_beg': iatm_beg}
        return br_id


    def add_ia_vdw(self, atm_nam_i, atm_nam_j, style, params):
        '''
        All atom types must be specified before calling this function.

        '''
        num_atom_types = len(self.atom_types)

        if self.vdw is None:
            self.vdw = {}
            ivdwt = 1
            nat = len(self.atom_types)
            for jat in range(1, nat+1):
                for iat in range(jat, nat+1):
                    self.vdw[ivdwt] = {'type_i': at_i, 'type_j': at_j, 'style': 0}
                    ivdwt += 1

        for iat in self.atom_types.keys():
            if self.atom_types[iat]['name'] == atm_nam_i:
                at_i = iat
                break

        for jat in self.atom_types.keys():
            if self.atom_types[jat]['name'] == atm_nam_j:
                at_j = jat
                break

        if at_i < at_j:
            ivdwt = at_j + (2*num_atom_types-at_i)*(at_i-1)//2
        else:
            ivdwt = at_i + (2*num_atom_types-at_j)*(at_j-1)//2

        self.vdw[ivdwt] = {'type_i': at_i, 'type_j': at_j,
                'style': style, 'params': np.copy(params)}
        return ivdwt


    def add_external(self, style, params):
        '''
        All atom types must be specified before calling this function.

        '''
        iext = len(self.external) + 1
        self.external[iext] = {'style': style, 'params': np.copy(params)}
        return iext


    def add_flow_field(self, style, params):
        '''
        Style 0: no flow; params: empty
        Style 1: shear; params[0]: shear rate
        Style 2: planar extension; params[0]: extension rate
        Style 3: uniaxial extension; params[0]: extension rate

        '''
        self.flow_field = {'style': style, 'params': np.copy(params)}


    def add_mpcd_atoms(self, avnc):
        '''
        Final box size must be set before calling this function.

        '''
        vol = self.simbox.prod()
        self.num_mpcd_atoms = math.ceil(avnc*vol)
        self.mpcd_avnc = avnc


    def to_center(self, imol):
        '''
        Bring the c.o.m of molecule imol to the box center. All atoms must be
        added before calling this function. This may be useful for single
        untethered molecules in an unbounded domain.

        '''
        com = np.zeros((3,))
        n = self.molecules[imol]['natm'] 
        ibeg = self.molecules[imol]['iatm_beg'] 
        iend = ibeg + n

        for i in range(ibeg, iend):
            coords = self.atoms[i]['coords']
            com += coords
        com /= n
        for i in range(ibeg, iend):
            self.atoms[i]['coords'] += (-com + 0.5*self.simbox)


    def apply_pbc(self):
        '''
        Applies PBC along all directions.

        '''
        assert self.imcon == 1
        for i in range(1,len(self.atoms)+1):
            coords = self.atoms[i]['coords']
            coords -= self.simbox*np.floor(coords/self.simbox)


    def add_random_charges(self):
        '''
        Add random charges to all atoms while keeping the system neutral.
        Charges are in units of +/-1.0.

        '''
        chge = np.zeros( (len(self.atoms),), dtype=np.int32 )
        n = len(self.atoms)
        hafn = n//2
        if n%2 == 0:
            #For even number of atoms
            chge[:hafn] = np.random.randint(-1, 2, hafn)
            chge[hafn:] = -chge[:hafn]
        else:
            #For odd number of atoms
            chge[:hafn] = np.random.randint(-1, 2, hafn)
            chge[hafn+1:] = -chge[:hafn]

        assert chge.sum() == 0

        for i in range(1,len(self.atoms)+1):
            self.atoms[i]['charge'] = chge[i-1]


    def update_simbox(self):
        '''
        Updates the simulation box dimensions such that all atoms are contained
        within it. Does not change boundary conditions as set before.
        `self.add_simbox` must be called and all atoms already added 
        before calling this function.
        Useful for setting the box size when it is hard to guess or calculate
        box size, e.g. for single molecules in unbounded domain.
        Note: If tethered molecules are not entirely within the first octant,
        they will be translated as well(may result in incorrect connectivity).

        '''
        #Move each molecule to the first octant.
        for imol in range(1, len(self.molecules)+1):
            natm = self.molecules[imol]['natm']
            ibeg = self.molecules[imol]['iatm_beg']
            iend = ibeg + natm

            minx = 0.0; miny = 0.0; minz = 0.0
            for iatm in range(ibeg, iend):
                coords = self.atoms[iatm]['coords']
                minx = min(minx, coords[0])
                miny = min(miny, coords[1])
                minz = min(minz, coords[2])
            dr = np.zeros((3,))
            if minx < 0.0: dr[0] = abs(minx)
            if miny < 0.0: dr[1] = abs(miny)
            if minz < 0.0: dr[2] = abs(minz)
            for iatm in range(ibeg, iend):
                coords = self.atoms[iatm]['coords']
                coords += dr
        #Update box dimensions
        maxx = 0.0; maxy = 0.0; maxz = 0.0
        for iatm in range(1, len(self.atoms)+1):
            coords = self.atoms[iatm]['coords']
            maxx = max(maxx, coords[0])
            maxy = max(maxy, coords[1])
            maxz = max(maxz, coords[2])
        Lmax = max(maxx, maxy, maxz)
        if Lmax > self.simbox.max():
            self.simbox[:] = Lmax


    def has_overlap(self, ri, sep):
        '''
        Checks if an atom at position `ri` overlaps with existing atoms within a
        separation distance of `sep`.
        '''
        out = False
        sepsq = sep*sep
        for j in self.atoms.keys():
            rj = self.atoms[j]['coords']
            rij = rj - ri
            rijmsq = np.dot(rij, rij)
            if rijmsq < sepsq:
                out = True
                break
        return out
