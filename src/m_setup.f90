module m_setup
    !! Routines for doing allocation, etc. in preparation for simulation run.

use m_precision
use m_ran_num
use m_simbox
use m_globals
use m_config_io
use m_interaction
use m_bd_solver
use m_mpcd
use m_relax
use m_stats_io

implicit none

contains

!*******************************************************************************

subroutine setup()

    integer :: ierr

    ierr = 0

    !Initialize random number generator
    if (read_seed) then
        call init_stream('random_seed.txt'//trim(adjustl(job_tag)))
    else 
        call init_stream('')
    end if
    if (write_seed) call save_seed('random_seed.txt'//trim(adjustl(job_tag)))

    !Create the simulation box. Its attributes will change later based on input
    !data.
    call smbx_init(simbox)

    !Read configuration from files
    if ( lrevive ) then
        !Restarting simulation: Read revive file
        call read_dump(fn_revive//trim(adjustl(job_tag)))
    else 
        !New simulation: Read initial configuration file
        call read_config(fn_cfg//trim(adjustl(job_tag)))
        nts = 0
    end if

    !Allocate memory for forces. Forces are not saved in revive file or
    !in config file as they can be calculated from position data.
    !MPCD particles do not have any forces acting on them.
    allocate( forces(3,num_atoms) )
    forces = 0.0_rp

    !Initialize stats collection
    call stats_init()

    !Set up interactions
    call ia_setup()

    !Set up solver. No need to set up for structure relaxation.
    if (sim_style == 1) then
        call bds_init(ierr)
    else if (sim_style == 2) then
        call mpcd_init(ierr)
    end if
    if (ierr /= 0) call finish()

    end subroutine

!*******************************************************************************

subroutine run()

    if (sim_style == 0) then
        call rlx_run()
    else if (sim_style == 1) then
        call bds_run()
    else if (sim_style == 2) then
        call mpcd_run()
    end if

    end subroutine

!*******************************************************************************

subroutine config_clear()
    !! Clears out all configuration related variables in module `m_globals`.

    num_atom_types = 0; num_atoms = 0
    if (allocated(atom_names))  deallocate(atom_names)
    if (allocated(atom_styles)) deallocate(atom_styles)
    if (allocated(atom_mass  )) deallocate(atom_mass)
    if (allocated(atoms      )) deallocate(atoms)
    if (allocated(coordinates)) deallocate(coordinates)
    if (allocated(velocities )) deallocate(velocities)
    if (allocated(forces     )) deallocate(forces)
    if (allocated(charge     )) deallocate(charge)

    num_bond_types = 0; num_bonds = 0
    if (allocated(bond_styles)) deallocate(bond_styles)
    if (allocated(bond_params)) deallocate(bond_params)
    if (allocated(bonds      )) deallocate(bonds)

    num_angle_types = 0; num_angles = 0
    if (allocated(angle_styles)) deallocate(angle_styles)
    if (allocated(angle_params)) deallocate(angle_params)
    if (allocated(angles      )) deallocate(angles)

    num_dihedral_types = 0; num_dihedrals = 0
    if (allocated(dihedral_styles)) deallocate(dihedral_styles)
    if (allocated(dihedral_params)) deallocate(dihedral_params)
    if (allocated(dihedrals      )) deallocate(dihedrals)

    num_branches = 0
    if (allocated(branches)) deallocate(branches)

    num_molecule_types = 0; num_molecules = 0
    if (allocated(molecule_names)) deallocate(molecule_names)
    if (allocated(molecule_pop  )) deallocate(molecule_pop)
    if (allocated(molecules     )) deallocate(molecules)

    num_tether_types = 0; num_tethers = 0
    if (allocated(tether_styles)) deallocate(tether_styles)
    if (allocated(tether_params)) deallocate(tether_params)
    if (allocated(tethers      )) deallocate(tethers)
    if (allocated(tether_points)) deallocate(tether_points)

    num_vdw_types = 0
    if (allocated(vdw_styles)) deallocate(vdw_styles)
    if (allocated(vdw_params)) deallocate(vdw_params)
    if (allocated(vdw_pairs )) deallocate(vdw_pairs)

    num_externals = 0
    if (allocated(external_styles)) deallocate(external_styles)
    if (allocated(external_params)) deallocate(external_params)

    if (allocated(flow_params)) deallocate(flow_params)

    end subroutine

!*******************************************************************************

subroutine finish()

    if (sim_style == 1) then
        call bds_finish()
    else if (sim_style == 2) then
        call mpcd_finish()
    end if

    call ia_finish()
    call stats_finish()
    call config_clear()

    end subroutine

!*******************************************************************************

end module m_setup
