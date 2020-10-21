!********************************************************************************!
!                                                                                !
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2020 Sarit Dutta                                                 !
!                                                                                !
! Permission is hereby granted, free of charge, to any person obtaining a copy   !
! of this software and associated documentation files (the "Software"), to deal  !
! in the Software without restriction, including without limitation the rights   !
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      !
! copies of the Software, and to permit persons to whom the Software is          !
! furnished to do so, subject to the following conditions:                       !
!                                                                                !
! The above copyright notice and this permission notice shall be included in all !
! copies or substantial portions of the Software.                                !
!                                                                                !
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     !
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       !
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    !
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         !
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  !
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  !
! SOFTWARE.                                                                      !
!                                                                                !
!********************************************************************************!

module atmcfg_m

use constants_m

implicit none

integer, parameter :: mxparam = 12
    !! Maximum number of parameters for bonds, angles, etc.

type atmcfg_t
    !Particle configuration: Atoms
    integer :: num_atom_types = 0
        !! Number of *atom_type*s
    character(len=8), dimension(:), allocatable :: atom_names
        !! (*num_atom_types*,) array. Name of atoms of each type.
    integer, dimension(:), allocatable :: atom_styles
        !! (*num_atom_types*,) array. Style of atoms of each type.
    real(rp), dimension(:), allocatable :: atom_mass
        !! (*num_atom_types*,) array. Mass of atoms of each type.
    integer :: num_atoms = 0
        !!  Number of atoms
    integer , dimension(:), allocatable:: atoms
        !! (*num_atoms*,) array. 
        !!
        !! For atom *i*, its type *at = atoms(i)*, with style
        !! *atom_styles(at)*, name *atom_names(at)*, mass *atom_mass(at)*, charge
        !! *charge(i)*, position *coordinates(:,i)*, velocity *velocities(:,i)*,
        !! orientation (if the style requires) *orientations(:,i)*. The force 
        !! acting on atom *i* is *forces(:,i)*.
    real(rp), dimension(:), allocatable :: charge
        !! (*num_atoms*,) array.
    real(rp), dimension(:,:), allocatable :: coordinates
        !!  (3, *num_atoms*) array
    real(rp), dimension(:,:), allocatable :: forces
        !!  (3, *num_atoms*) array
    
    !Particle configuration: Bonds
    integer :: num_bond_types = 0
        !!  Number of *bond_type*s
    integer, dimension(:), allocatable :: bond_styles
        !!  (*num_bond_types*,) array.
    real(rp), dimension(:,:), allocatable :: bond_params
        !!  (*mxparam*,*num_bond_types*) array.
    integer :: num_bonds = 0
        !!  Total number of bonds.
    integer, dimension(:,:), allocatable :: bonds
        !! (3, *num_bonds*) array. Bond *i* is of type *bt = bonds(1,i)*,  directed from
        !! atom *bonds(2,i)* to *bonds(3,i)*. Its style is *bond_styles(bt)* with
        !! parameters *bond_params(:,bt)*.
    
    !Particle configuration: Angles
    integer :: num_angle_types = 0
        !!  Number of *angle_type*s
    integer, dimension(:), allocatable :: angle_styles
        !!  (*num_angle_types*,) array
    real(rp), dimension(:,:), allocatable :: angle_params
        !!  (*mxparam*, *num_angle_types*) array
    integer :: num_angles = 0
        !!  Number of angles
    integer, dimension(:,:), allocatable :: angles
        !! (4, *num_angles*) array. Angle *i* is of type *ant = angles(1,i)*, incident
        !! to atoms *angles(2,i)*, *angles(3,i)*, and *angles(4,i)*. Its style is
        !! *angle_styles(ant)* with parameters *angle_params(:,ant)*.
    
    !Particle configuration: Dihedrals
    integer :: num_dihedral_types = 0
        !!  Number of *dihedral_type*s
    integer, dimension(:), allocatable :: dihedral_styles
        !!  (*num_dihedral_types*,) array
    real(rp), dimension(:,:), allocatable :: dihedral_params
        !!  (*mxparam*, *num_dihedral_types*) array
    integer :: num_dihedrals = 0
        !!  Number of dihedrals
    integer, dimension(:,:), allocatable :: dihedrals
        !! (5, *num_dihedrals*) array. Dihedral *i* is of type *dt = dihedrals(1,i)*, incident
        !! to atoms *dihedrals(2,i)*, *dihedrals(3,i)*, *dihedrals(4,i)*, and *dihedrals(5,i)*.
        !! Its style is *dihedral_styles(dt)* with parameters *dihedral_params(:,dt)*.
    
    !Particle configuration: Branches
    integer :: num_branches = 0
        !! Total number of branches (including the backbone)
    integer , dimension(:,:), allocatable:: branches
        !! (3,*num_branches*) array. Branch *i* is tethered to atom *branches(1,i)*,
        !! contains *branches(2,i)* atoms, with the beginning atom index *branches(3,i)*.
    
    !Particle configuration: Molecules
    integer :: num_molecule_types = 0
        !!  Number of *molecule_type*s
    character(len=8), dimension(:), allocatable :: molecule_names
        !! (*num_molecule_types*,) array
    integer, dimension(:), allocatable :: molecule_pop
        !! (*num_molecule_types*,) array
    integer :: num_molecules = 0
        !!  Number of molecules
    integer, dimension(:,:), allocatable :: molecules
        !! (9,*num_molecules*) array. For molecule *i*, its type *mt = molecules(1,i)*, 
        !! containing *molecules(2,i)* atoms with beginning index *molecules(3,i)*,
        !! *molecules(4,i)* bonds with beginning index *molecules(5,i)*,
        !! *molecules(6,i)* angles with beginning index *molecules(7,i)*, and
        !! *molecules(8,i)* dihedrals with beginning index *molecules(9,i)*.
    real(rp), dimension(3) :: molecule_com = 0.0_rp
        !! Center of mass of the molecule. This is used only when imcon == 0, i.e.
        !! for a single molecule without periodic boundaries.
    
    !Particle configuration: Tethers
    integer :: num_tether_types = 0
        !!  Number of *tether_type*s
    integer, dimension(:), allocatable :: tether_styles
        !!  (*num_tether_types*,) array
    real(rp), dimension(:,:), allocatable :: tether_params
        !!  (*mxparam*, *num_tether_types*) array
    integer :: num_tethers = 0
        !!  Number of tethers
    integer, dimension(:,:), allocatable :: tethers
        !! (2, *num_tethers*) array. Tether *i* is of type *tt = tethers(1,i)*, tethering
        !! atom *tethers(2,i)* to a point *tether_points(:,i)*.
        !! Its style is *tether_styles(tt)* with parameters *tether_params(:,tt)*.
    real(rp), dimension(:,:), allocatable :: tether_points
        !!  (3, *num_tethers*) array
    
    !Particle configuration: VDW (pair) interactions
    integer :: num_vdw_types = 0
        !!  Number of *vdw_type*s
    integer, dimension(:), allocatable :: vdw_styles
        !!  (*num_vdw_types*,) array
    real(rp), dimension(:,:), allocatable :: vdw_params
        !!  (*mxparam*, *num_vdw_types*) array
    integer, dimension(:,:), allocatable :: vdw_pairs
        !!  (2, *num_vdw_types*) array. Stores atom type of interacting pairs, such
        !! that at_i >= at_j.
    
    !Particle configuration: External force field
    integer :: num_externals = 0
        !!  Number of external fields
    integer, dimension(:), allocatable :: external_styles
        !!  (*num_external*,) array
    real(rp), dimension(:,:), allocatable :: external_params
        !!  (*mxparam*, *num_external*) array
    
    !Particle configuration: Flow field
    integer :: flow_style = 0
    real(rp), dimension(:), allocatable :: flow_params
        !!  (*mxparam*,) array
end type atmcfg_t

contains

!******************************************************************************

subroutine atmcfg_delete(this)
    !! Deallocates all memory acquired by a `configuration_t` object and resets
    !! all other components to zero. Exception: `num_coeffs` is not reset to
    !! zero.

    type(atmcfg_t), intent(in out) :: this

    this%num_atom_types = 0; this%num_atoms = 0
    if (allocated(this%atom_names))  deallocate(this%atom_names)
    if (allocated(this%atom_styles)) deallocate(this%atom_styles)
    if (allocated(this%atom_mass  )) deallocate(this%atom_mass)
    if (allocated(this%atoms      )) deallocate(this%atoms)
    if (allocated(this%coordinates)) deallocate(this%coordinates)
    if (allocated(this%forces     )) deallocate(this%forces)
    if (allocated(this%charge     )) deallocate(this%charge)

    this%num_bond_types = 0; this%num_bonds = 0
    if (allocated(this%bond_styles)) deallocate(this%bond_styles)
    if (allocated(this%bond_params)) deallocate(this%bond_params)
    if (allocated(this%bonds      )) deallocate(this%bonds)

    this%num_angle_types = 0; this%num_angles = 0
    if (allocated(this%angle_styles)) deallocate(this%angle_styles)
    if (allocated(this%angle_params)) deallocate(this%angle_params)
    if (allocated(this%angles      )) deallocate(this%angles)

    this%num_dihedral_types = 0; this%num_dihedrals = 0
    if (allocated(this%dihedral_styles)) deallocate(this%dihedral_styles)
    if (allocated(this%dihedral_params)) deallocate(this%dihedral_params)
    if (allocated(this%dihedrals      )) deallocate(this%dihedrals)

    this%num_branches = 0
    if (allocated(this%branches)) deallocate(this%branches)

    this%num_molecule_types = 0; this%num_molecules = 0
    if (allocated(this%molecule_names)) deallocate(this%molecule_names)
    if (allocated(this%molecule_pop  )) deallocate(this%molecule_pop)
    if (allocated(this%molecules     )) deallocate(this%molecules)

    this%num_tether_types = 0; this%num_tethers = 0
    if (allocated(this%tether_styles)) deallocate(this%tether_styles)
    if (allocated(this%tether_params)) deallocate(this%tether_params)
    if (allocated(this%tethers      )) deallocate(this%tethers)
    if (allocated(this%tether_points)) deallocate(this%tether_points)

    this%num_vdw_types = 0
    if (allocated(this%vdw_styles)) deallocate(this%vdw_styles)
    if (allocated(this%vdw_params)) deallocate(this%vdw_params)
    if (allocated(this%vdw_pairs )) deallocate(this%vdw_pairs)

    this%num_externals = 0
    if (allocated(this%external_styles)) deallocate(this%external_styles)
    if (allocated(this%external_params)) deallocate(this%external_params)

    if (allocated(this%flow_params)) deallocate(this%flow_params)

    end subroutine

!******************************************************************************

end module atmcfg_m
