module m_config_io
    !! Routines for IO of config and dump files.

use m_precision
use m_strings
use m_simbox
use m_globals

implicit none

contains

!******************************************************************************

subroutine read_dump(fn)
    !! Reads from DUMP file

    character(len=*), intent(in) :: fn
    real(rp), dimension(3,3) :: mat
    integer :: fu
    logical :: lvel

    open(newunit=fu, file=fn, access='stream', form='unformatted', &
        action='read', status='old')

    read(fu) leql, nts

    read(fu) mat
    call simbox%set_basis(mat)
    read(fu) imcon

    read(fu) num_atom_types
    allocate(atom_names (num_atom_types))
    allocate(atom_styles(num_atom_types))
    allocate(atom_mass  (num_atom_types))
    read(fu) atom_names, atom_styles, atom_mass

    read(fu) mpcd_avnc, num_mpcd_atoms
    read(fu) num_atoms, num_atoms_tot

    allocate(atoms(num_atoms))
    read(fu) atoms

    allocate(coordinates(3,num_atoms_tot))
    read(fu) coordinates

    read(fu) lvel
    if (lvel) then
        allocate(velocities(3,num_atoms_tot))
        read(fu) velocities
    end if

    allocate(charge(num_atoms))
    read(fu) charge

    read(fu) num_bond_types
    if (num_bond_types > 0) then
        allocate(bond_styles(num_bond_types))
        allocate(bond_params(mxparam,num_bond_types))
        read(fu) bond_styles, bond_params
    end if

    read(fu) num_bonds
    if (num_bonds > 0) then
        allocate(bonds(3,num_bonds))
        read(fu) bonds
    end if

    read(fu) num_angle_types
    if (num_angle_types > 0) then
        allocate(angle_styles(num_angle_types))
        allocate(angle_params(mxparam,num_angle_types))
        read(fu) angle_styles, angle_params
    end if

    read(fu) num_angles
    if (num_angles > 0) then
        allocate(angles(4,num_angles))
        read(fu) angles
    end if

    read(fu) num_dihedral_types
    if (num_dihedral_types > 0) then
        allocate(dihedral_styles(num_dihedral_types))
        allocate(dihedral_params(mxparam,num_dihedral_types))
        read(fu) dihedral_styles, dihedral_params
    end if

    read(fu) num_dihedrals
    if (num_dihedrals > 0) then
        allocate(dihedrals(5,num_dihedrals))
        read(fu) dihedrals
    end if

    read(fu) num_branches
    if (num_branches > 0) then
        allocate(branches(3,num_branches))
        read(fu) branches
    end if

    read(fu) num_molecule_types
    allocate(molecule_names(num_molecule_types))
    allocate(molecule_pop  (num_molecule_types))
    read(fu) molecule_names, molecule_pop
    read(fu) num_molecules
    allocate(molecules(9,num_molecules))
    read(fu) molecules

    read(fu) num_tether_types
    if (num_tether_types > 0) then
        allocate(tether_styles(num_tether_types))
        allocate(tether_params(mxparam, num_tether_types))
        read(fu) tether_styles, tether_params
    end if

    read(fu) num_tethers
    if (num_tethers > 0) then
        allocate(tethers(2,num_tethers))
        allocate(tether_points(3,num_tethers))
        read(fu) tethers, tether_points
    end if

    read(fu) num_vdw_types
    if (num_vdw_types > 0) then
        allocate(vdw_styles(num_vdw_types))
        allocate(vdw_pairs(2,num_vdw_types))
        allocate(vdw_params(mxparam,num_vdw_types))
        read(fu) vdw_pairs, vdw_styles, vdw_params
    end if

    read(fu) num_externals
    if (num_externals > 0) then
        allocate(external_styles(num_externals))
        allocate(external_params(mxparam,num_externals))
        read(fu) external_styles, external_params
    end if

    read(fu) flow_style
    if (flow_style /= 0) then
        allocate(flow_params(mxparam))
        read(fu) flow_params
    end if

    if ( (num_tethers==0) .and. (imcon==0) ) read(fu) molecule_com

    close(fu)

    end subroutine

!******************************************************************************

subroutine write_dump(fn)
    !! Writes to DUMP file.

    character(len=*), intent(in) :: fn
    integer :: fu
    logical :: lvel

    open(newunit=fu, file=fn, access='stream', form='unformatted', &
        action='write', status='replace')

    write(fu) leql, nts

    write(fu) simbox%basis
    write(fu) imcon

    write(fu) num_atom_types, atom_names, atom_styles, atom_mass
    write(fu) mpcd_avnc, num_mpcd_atoms
    write(fu) num_atoms, num_atoms_tot

    write(fu) atoms
    write(fu) coordinates
    if (allocated(velocities)) then
        !A logical indicating whether velocities are written
        lvel = .true.
        write(fu) lvel, velocities
    else
        lvel = .false.
        write(fu) lvel
    end if
    write(fu) charge

    write(fu) num_bond_types
    if (num_bond_types > 0) write(fu) bond_styles, bond_params
    write(fu) num_bonds
    if (num_bonds > 0) write(fu) bonds

    write(fu) num_angle_types
    if (num_angle_types > 0) write(fu) angle_styles, angle_params
    write(fu) num_angles
    if (num_angles > 0) write(fu) angles

    write(fu) num_dihedral_types
    if (num_dihedral_types > 0) write(fu) dihedral_styles, dihedral_params
    write(fu) num_dihedrals
    if (num_dihedrals > 0) write(fu) dihedrals

    write(fu) num_branches
    if (num_branches > 0) write(fu) branches

    write(fu) num_molecule_types, molecule_names, molecule_pop
    write(fu) num_molecules, molecules

    write(fu) num_tether_types
    if (num_tether_types > 0) write(fu) tether_styles, tether_params
    write(fu) num_tethers
    if (num_tethers > 0) write(fu) tethers, tether_points

    write(fu) num_vdw_types
    if (num_vdw_types > 0) write(fu) vdw_pairs, vdw_styles, vdw_params

    write(fu) num_externals
    if (num_externals > 0) write(fu) external_styles, external_params

    write(fu) flow_style
    if (flow_style /= 0) write(fu) flow_params

    if ( (num_tethers==0) .and. (imcon==0) ) write(fu) molecule_com

    close(fu)

    end subroutine

!******************************************************************************

subroutine read_config(fn)
    !! Read from CONFIG file

    character(len=*), intent(in) :: fn
    character(len=mxrdln) :: line
    character(len=:), allocatable :: word
    real(rp), dimension(3,3) :: mat
    integer :: i, it, npar, ibr
    integer :: fu, ios

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '#', ios)
        if (ios /= 0) exit

        line = adjustl(line)

        if (str_startswith(line, 'SIMBOX')) then
            !Read simulation box size into lattice vectors defining the box
            read(fu,*) mat(:,1)
            read(fu,*) mat(:,2)
            read(fu,*) mat(:,3)
            call simbox%set_basis(mat)
        end if

        if (str_startswith(line, 'IMCON')) then
            !Read box boundary condition
            call str_split(line, ' ', word)
            imcon = str_to_i(line)
        end if

        if (str_startswith(line, 'ATOM_TYPES')) then
            call str_split(line, ' ', word)
            num_atom_types = str_to_i(line)

            allocate(atom_names (num_atom_types))
            allocate(atom_styles(num_atom_types))
            allocate(atom_mass  (num_atom_types))

            do it = 1, num_atom_types
                read(fu,*) atom_names(it), atom_styles(it), atom_mass(it)
            end do
        end if

        if (str_startswith(line, 'MPCD_ATOMS')) then
            call str_split(line, ' ', word)
            call str_split(line, ' ', word)
            num_mpcd_atoms = str_to_i(word)
            mpcd_avnc = str_to_i(line)
        end if

        if (str_startswith(line, 'ATOMS')) then
            call str_split(line, ' ', word)
            num_atoms = str_to_i(line)

            num_atoms_tot = num_mpcd_atoms + num_atoms
            allocate(atoms (num_atoms))
            allocate(charge (num_atoms))
            allocate(coordinates(3,num_atoms_tot))

            do i = 1, num_atoms
                read(fu,*) atoms(i), charge(i), coordinates(:,i)
            end do
        end if

        if (str_startswith(line, 'BOND_TYPES')) then
            call str_split(line, ' ', word)
            num_bond_types = str_to_i(line)

            allocate(bond_styles(num_bond_types))
            allocate(bond_params(mxparam, num_bond_types))
            bond_params = 0.0_rp

            do it = 1, num_bond_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                bond_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) bond_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'BONDS')) then
            call str_split(line, ' ', word)
            num_bonds = str_to_i(line)
            allocate(bonds(3,num_bonds))
            do i = 1, num_bonds
                read(fu,*) bonds(:,i)
            end do
        end if

        if (str_startswith(line, 'ANGLE_TYPES')) then
            call str_split(line, ' ', word)
            num_angle_types = str_to_i(line)

            allocate(angle_styles(num_angle_types))
            allocate(angle_params(mxparam, num_angle_types))
            angle_params = 0.0_rp

            do it = 1, num_angle_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                angle_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) angle_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'ANGLES')) then
            call str_split(line, ' ', word)
            num_angles = str_to_i(line)
            allocate(angles(4,num_angles))
            do i = 1, num_angles
                read(fu,*) angles(:,i)
            end do
        end if

        if (str_startswith(line, 'DIHEDRAL_TYPES')) then
            call str_split(line, ' ', word)
            num_dihedral_types = str_to_i(line)

            allocate(dihedral_styles(num_dihedral_types))
            allocate(dihedral_params(mxparam, num_dihedral_types))
            dihedral_params = 0.0_rp

            do it = 1, num_dihedral_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                dihedral_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) dihedral_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'DIHEDRALS')) then
            call str_split(line, ' ', word)
            num_dihedrals = str_to_i(line)
            allocate(dihedrals(5,num_dihedrals))
            do i = 1, num_dihedrals
                read(fu,*) dihedrals(:,i)
            end do
        end if

        if (str_startswith(line, 'BRANCHES')) then
            call str_split(line, ' ', word)
            num_branches = str_to_i(line)
            allocate(branches(3,num_branches))
            do ibr = 1, num_branches
                read(fu,*) branches(:,ibr)
            end do
        end if

        if (str_startswith(line, 'MOLECULE_TYPES')) then
            call str_split(line, ' ', word)
            num_molecule_types = str_to_i(line)

            allocate(molecule_names(num_molecule_types))
            allocate(molecule_pop(num_molecule_types))

            do it = 1, num_molecule_types
                read(fu, *) molecule_names(it), molecule_pop(it)
            end do
        end if

        if (str_startswith(line, 'MOLECULES')) then
            call str_split(line, ' ', word)
            num_molecules = str_to_i(line)
            allocate(molecules(9,num_molecules))
            do i = 1, num_molecules
                read(fu,*) molecules(:,i)
            end do
        end if

        if (str_startswith(line, 'TETHER_TYPES')) then
            call str_split(line, ' ', word)
            num_tether_types = str_to_i(line)

            allocate(tether_styles(num_tether_types))
            allocate(tether_params(mxparam, num_tether_types))
            tether_params = 0.0_rp

            do it = 1, num_tether_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                tether_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) tether_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'TETHERS')) then
            call str_split(line, ' ', word)
            num_tethers = str_to_i(line)
            allocate(tethers(2,num_tethers))
            allocate(tether_points(3,num_tethers))
            do i = 1, num_tethers
                read(fu,*) tethers(:,i), tether_points(:,i)
            end do
        end if

        if (str_startswith(line, 'VDW')) then
            call str_split(line, ' ', word)
            num_vdw_types = str_to_i(line)

            allocate(vdw_styles(num_vdw_types))
            allocate(vdw_pairs(2,num_vdw_types))
            allocate(vdw_params(mxparam,num_vdw_types))
            vdw_params = 0.0_rp

            do it = 1, num_vdw_types
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                vdw_pairs(1,it) = str_to_i(word)
                call str_split(line, ' ', word)
                vdw_pairs(2,it) = str_to_i(word)
                call str_split(line, ' ', word)
                vdw_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) vdw_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'EXTERNAL')) then
            call str_split(line, ' ', word)
            num_externals = str_to_i(line)

            allocate(external_styles(num_externals))
            allocate(external_params(mxparam,num_externals))
            external_params = 0.0_rp

            do it = 1, num_externals
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                external_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) external_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'FLOW_FIELD')) then
            allocate(flow_params(mxparam))
            flow_params = 0.0_rp
            call readline(fu, line, '#', ios)
            call str_split(line, ' ', word)
            flow_style = str_to_i(word)
            if (flow_style > 0) then
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) flow_params(1:npar)
                end if
            end if
        end if
    end do

    close(fu)

    end subroutine

!*******************************************************************************

subroutine write_config(fn, title)
    !! Write to cfg file

    character(len=*), intent(in) :: fn
    character(len=*), intent(in) :: title
    integer :: fu
    integer :: i, it

    open(newunit=fu, file=fn, action='write', status='replace')

    write(fu, '(a)') '#'//trim(adjustl(title))
    write(fu, '(a)') 'version 1.0'

    write(fu, *)
    write(fu, '(a)') 'SIMBOX'
    write(fu, '(3(g0.8,2x))') simbox%basis(:,1)
    write(fu, '(3(g0.8,2x))') simbox%basis(:,2)
    write(fu, '(3(g0.8,2x))') simbox%basis(:,3)
    write(fu, *)
    write(fu, '(a,2x,i0)') 'IMCON', imcon

    write(fu, *)
    write(fu, '(a,2x,i0)') 'ATOM_TYPES', num_atom_types
    do it = 1, num_atom_types
        write(fu, '(a,2x,i0,2x,g0.6)') trim(adjustl(atom_names(it))), &
            atom_styles(it), atom_mass(it)
    end do

    if (num_mpcd_atoms > 0) then
        write(fu, *)
        write(fu, '(a,2x, i0,2x,i0)') 'MPCD_ATOMS', num_mpcd_atoms, mpcd_avnc
    end if

    write(fu, *)
    write(fu, '(a,2x,i0)') 'ATOMS', num_atoms
    do i = 1, num_atoms
        write(fu, '(i0,2x,*(g0.14,2x))') atoms(i), charge(i), coordinates(:,i)
    end do

    if (num_bonds > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'BOND_TYPES', num_bond_types
        do it = 1, num_bond_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') bond_styles(it), &
                size(bond_params,1), bond_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'BONDS', num_bonds
        do i = 1, num_bonds
            write(fu, '(*(i0,2x))') bonds(:,i)
        end do
    end if

    if (num_angles > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'ANGLE_TYPES', num_angle_types
        do it = 1, num_angle_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') angle_styles(it), &
                size(angle_params,1), angle_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'ANGLES', num_angles
        do i = 1, num_angles
            write(fu, '(*(i0,2x))') angles(:,i)
        end do
    end if

    if (num_dihedrals > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'DIHEDRAL_TYPES', num_dihedral_types
        do it = 1, num_dihedral_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') dihedral_styles(it), &
                size(dihedral_params,1), dihedral_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'DIHEDRALS', num_dihedrals
        do i = 1, num_dihedrals
            write(fu, '(*(i0,2x))') dihedrals(:,i)
        end do
    end if

    if (num_branches > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'BRANCHES', num_branches
        do i = 1, num_branches
            write(fu, '(*(i0,2x))') branches(:,i)
        end do
    end if

    write(fu, *)
    write(fu, '(a,2x,i0)') 'MOLECULE_TYPES', num_molecule_types
    do it = 1, num_molecule_types
        write(fu, '(a,2x,i0)') trim(adjustl(molecule_names(it))), molecule_pop(it)
    end do

    write(fu, *)
    write(fu, '(a,2x,i0)') 'MOLECULES', num_molecules
    do i = 1, num_molecules
        write(fu, '(*(i0,2x))') molecules(:,i)
    end do

    if (num_tethers > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'TETHER_TYPES', num_tether_types
        do it = 1, num_tether_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') tether_styles(it), &
                size(tether_params,1), tether_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'TETHERS', num_tethers
        do i = 1, num_tethers
            write(fu, '(2(i0,2x),3(g0.6,2x))') tethers(:,i), tether_points(:,i)
        end do
    end if

    if (num_vdw_types > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'VDW', num_vdw_types
        do it = 1, num_vdw_types
            write(fu,'(2(i0,2x),i0,2x,i0,2x,*(g0.6,2x))') vdw_pairs(:,it), &
                vdw_styles(it), size(vdw_params,1), vdw_params(:,it)
        end do
    end if

    if (num_externals > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'EXTERNAL', num_externals
        do it = 1, num_externals
            write(fu,'(i0,2x,i0,2x,*(g0.6,2x))') external_styles(it), &
                size(external_params,1), external_params(:,it)
        end do
    end if

    write(fu, *)
    write(fu, '(a)') 'FLOW_FIELD'
    write(fu,'(i0,2x,i0,2x,*(g0.6,2x))', advance='no') flow_style
    if (flow_style > 0) then 
        write(fu,'(i0,2x,*(g0.6,2x))') size(flow_params), flow_params
    else
        write(fu,*)
    end if

    close(fu)

    end subroutine

!******************************************************************************

subroutine write_ldf(fn_ld, title, with_mpcd_atoms)
    !! Write to a LAMMPS data file.

    character(len=*),       intent(in) :: fn_ld
        !! Name of the file
    character(len=*),       intent(in) :: title
        !! Title of the configuation
    logical, intent(in), optional :: with_mpcd_atoms
        !! Include MPCD atoms in the file? {T, F}
    real(rp), dimension(3) :: tilt_factors
    real(rp) :: chge
    logical :: with_mpcd_atoms_
    integer :: fu_ld
    integer :: cntr_atm, iatm_beg, natm
    integer :: i, iatm, imol

    if (present(with_mpcd_atoms)) then
        with_mpcd_atoms_ = with_mpcd_atoms
    else
        with_mpcd_atoms_ = .false.
    end if

    open(newunit=fu_ld, file=fn_ld, action='write')

    !Header
    write(fu_ld,'(a)') '#'//trim(adjustl(title))

    if (with_mpcd_atoms_) then
        write(fu_ld,'(i0,2x,a)') num_atoms_tot, 'atoms'
        write(fu_ld,'(i0,2x,a)') num_atom_types+1, 'atom types'
    else
        write(fu_ld,'(i0,2x,a)') num_atoms, 'atoms'
        write(fu_ld,'(i0,2x,a)') num_atom_types, 'atom types'
    end if

    if (num_bonds > 0) then
        write(fu_ld,'(i0,2x,a)') num_bonds, 'bonds'
        write(fu_ld,'(i0,2x,a)') num_bond_types, 'bond types'
    end if

    if (num_angles > 0) then
        write(fu_ld,'(i0,2x,a)') num_angles, 'angles'
        write(fu_ld,'(i0,2x,a)') num_angle_types, 'angle types'
    end if

    if (num_dihedrals > 0) then
        write(fu_ld,'(i0,2x,a)') num_dihedrals, 'dihedrals'
        write(fu_ld,'(i0,2x,a)') num_dihedral_types, 'dihedral types'
    end if

    !Simulation box & tilt factors
    write(fu_ld,'(a,2x,g0.6,2x,a)') '0.0', simbox%basis(1,1), 'xlo xhi'
    write(fu_ld,'(a,2x,g0.6,2x,a)') '0.0', simbox%basis(2,2), 'ylo yhi'
    write(fu_ld,'(a,2x,g0.6,2x,a)') '0.0', simbox%basis(3,3), 'zlo zhi'
    write(fu_ld,'(a)') '0.0 0.0 0.0 xy xz yz'

    !TODO: Triclinic box
    !   tilt_factors(1) = dot_product(simbox%basis(:,1)/norm2(simbox%basis(:,1)), &
    !       simbox%basis(:,2)/norm2(simbox%basis(:,2)))
    !   tilt_factors(2) = dot_product(simbox%basis(:,1)/norm2(simbox%basis(:,1)), &
    !       simbox%basis(:,3)/norm2(simbox%basis(:,3)))
    !   tilt_factors(3) = dot_product(simbox%basis(:,2)/norm2(simbox%basis(:,2)), &
    !       simbox%basis(:,3)/norm2(simbox%basis(:,3)))
    !   tilt_factors = acos(tilt_factors)
    !   write(fu_ld,'(3(g0.6,2x),a)') tilt_factors, 'xy xz yz'

    !Body: Atoms
    write(fu_ld,*)
    write(fu_ld,'(a)') 'Atoms # full'
    write(fu_ld,*)

    cntr_atm = 1
    do imol = 1, num_molecules
        natm = molecules(2,imol)
        iatm_beg = molecules(3,imol)
        do i = 1, natm
            iatm = iatm_beg + i -1
            write(fu_ld,'(i0,2x,i0,2x,i0,2x,4(g0.8,2x))') cntr_atm, imol, &
                atoms(iatm), charge(iatm), coordinates(:,iatm)
            cntr_atm = cntr_atm + 1
        end do
    end do

    if (with_mpcd_atoms_) then
        imol = 0; chge = 0.0_rp !No charge on MPCD atoms
        do i = num_atoms+1, num_atoms_tot
            write(fu_ld,'(i0,2x,i0,2x,i0,2x,4(g0.8,2x))') cntr_atm, imol, &
                (num_atom_types+1), chge, coordinates(:,i)
            cntr_atm = cntr_atm + 1
        end do
    end if

    !Body: Bonds
    if (num_bonds > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Bonds'
        write(fu_ld,*)
        do i = 1, num_bonds
            write(fu_ld,'(i0,2x,3(i0,2x))') i, bonds(:,i)
        end do
    end if

    !Body: Angles
    if (num_angles > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Angles'
        write(fu_ld,*)
        do i = 1, num_angles
            write(fu_ld,'(i0,2x,4(i0,2x))') i, angles(:,i)
        end do
    end if

    !Body: Dihedrals
    if (num_dihedrals > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Dihedrals'
        write(fu_ld,*)
        do i = 1, num_dihedrals
            write(fu_ld,'(i0,2x,5(i0,2x))') i, dihedrals(:,i)
        end do
    end if

    close(fu_ld)

    end subroutine

!*******************************************************************************

subroutine write_xyz(fn_xyz, title, with_mpcd_atoms)
    !! Write to an XYZ file.

    character(len=*),       intent(in) :: fn_xyz
        !! Name of the XYZ file
    character(len=*),       intent(in) :: title
        !! Title (for the configuration)
    logical, intent(in), optional :: with_mpcd_atoms
        !! Include MPCD atoms in the file? {T, F}
    integer :: iatm
    integer :: na
    integer :: fu_xyz

    if (present(with_mpcd_atoms)) then
        if (with_mpcd_atoms) then
            na = num_atoms_tot
        else
            na = num_atoms
        end if
    else
        na = num_atoms
    end if

    open(newunit=fu_xyz, file=fn_xyz, action='write')

    write(fu_xyz, '(i0)') na
    write(fu_xyz, '(a)') title

    do iatm = 1, na
        write(fu_xyz,'(3(g0.7,2x))') coordinates(:,iatm)
    end do

    close(fu_xyz)

    end subroutine

!*******************************************************************************

end module m_config_io
