module config_io_m
    !! Routines for IO of config and dump files.

use constants_m
use strings_m
use simbox_m
use atmcfg_m

implicit none

contains

!******************************************************************************

subroutine read_dump(nts, simbox, atc, fn)
    !! Reads from DUMP file

    integer(ip_long), intent(out) :: nts
    type(smbx_t),     intent(out) :: simbox
    type(atmcfg_t),   intent(out) :: atc
    character(len=*), intent(in)  :: fn
    real(rp), dimension(3,3) :: mat
    integer :: fu, nt, n

    open(newunit=fu, file=fn, access='stream', form='unformatted', &
        action='read', status='old')

    read(fu) nts

    read(fu) n, mat
    call simbox%init(n)
    call simbox%set_basis(mat)

    !Atom types data
    read(fu) nt
    atc%num_atom_types = nt
    allocate(atc%atom_names(nt))
    allocate(atc%atom_styles(nt))
    allocate(atc%atom_mass(nt))
    read(fu) atc%atom_names, atc%atom_styles, atc%atom_mass

    !Atom data
    read(fu) n
    atc%num_atoms = n
    allocate(atc%atoms(n));
    allocate(atc%coordinates(3,n))
    allocate(atc%charge(n))
    read(fu) atc%atoms, atc%coordinates, atc%charge

    !Bond types data
    read(fu) nt
    atc%num_bond_types = nt
    if (nt > 0) then
        allocate(atc%bond_styles(nt))
        allocate(atc%bond_params(mxparam,nt))
        read(fu) atc%bond_styles, atc%bond_params
    end if

    !Bond data
    read(fu) n
    atc%num_bonds = n
    if (n > 0) then
        allocate(atc%bonds(3,n))
        read(fu) atc%bonds
    end if

    !Angle types data
    read(fu) nt
    atc%num_angle_types = nt
    if (nt > 0) then
        allocate(atc%angle_styles(nt))
        allocate(atc%angle_params(mxparam,nt))
        read(fu) atc%angle_styles, atc%angle_params
    end if

    !Angles data
    read(fu) n
    atc%num_angles = n
    if (n > 0) then
        allocate(atc%angles(4,n))
        read(fu) atc%angles
    end if

    !Dihedral types data
    read(fu) nt
    atc%num_dihedral_types = nt
    if (nt > 0) then
        allocate(atc%dihedral_styles(nt))
        allocate(atc%dihedral_params(mxparam,nt))
        read(fu) atc%dihedral_styles, atc%dihedral_params
    end if

    !Dihedral data
    read(fu) n
    atc%num_dihedrals = n
    if (n > 0) then
        allocate(atc%dihedrals(5,n))
        read(fu) atc%dihedrals
    end if

    !Branches data
    read(fu) n
    atc%num_branches = n
    if (n > 0) then
        allocate(atc%branches(3,n))
        read(fu) atc%branches
    end if

    !Molecule types data
    read(fu) nt
    atc%num_molecule_types = nt
    allocate(atc%molecule_names(nt))
    allocate(atc%molecule_pop  (nt))
    read(fu) atc%molecule_names, atc%molecule_pop

    !Molecule data
    read(fu) n
    atc%num_molecules = n
    allocate(atc%molecules(9,n))
    read(fu) atc%molecules

    !Tether types data
    read(fu) nt
    atc%num_tether_types = nt
    if (nt > 0) then
        allocate(atc%tether_styles(nt))
        allocate(atc%tether_params(mxparam, nt))
        read(fu) atc%tether_styles, atc%tether_params
    end if

    !Tether data
    read(fu) n
    atc%num_tethers = n
    if (n > 0) then
        allocate(atc%tethers(2,n))
        allocate(atc%tether_points(3,n))
        read(fu) atc%tethers, atc%tether_points
    end if

    !vdw types data
    read(fu) nt
    atc%num_vdw_types = nt
    if (nt > 0) then
        allocate(atc%vdw_styles(nt))
        allocate(atc%vdw_pairs(2,nt))
        allocate(atc%vdw_params(mxparam,nt))
        read(fu) atc%vdw_pairs, atc%vdw_styles, atc%vdw_params
    end if

    !External field data
    read(fu) n
    atc%num_externals = n
    if (n > 0) then
        allocate(atc%external_styles(n))
        allocate(atc%external_params(mxparam,n))
        read(fu) atc%external_styles, atc%external_params
    end if

    !Flow field data
    read(fu) atc%flow_style
    if (atc%flow_style /= 0) then
        allocate(atc%flow_params(mxparam))
        read(fu) atc%flow_params
    end if

    if ( (atc%num_tethers==0) .and. (simbox%imcon==0) ) then
        read(fu) atc%molecule_com
    end if

    close(fu)

    end subroutine

!******************************************************************************

subroutine write_dump(nts, simbox, atc, fn)
    !! Writes to DUMP file.

    integer(ip_long), intent(in) :: nts
    type(smbx_t),     intent(in) :: simbox
    type(atmcfg_t),   intent(in) :: atc
    character(len=*), intent(in) :: fn
    integer :: fu

    open(newunit=fu, file=fn, access='stream', form='unformatted', &
        action='write', status='replace')

    write(fu) nts

    write(fu) simbox%imcon
    write(fu) simbox%basis

    write(fu) atc%num_atom_types
    write(fu) atc%atom_names, atc%atom_styles, atc%atom_mass

    write(fu) atc%num_atoms
    write(fu) atc%atoms, atc%coordinates, atc%charge

    write(fu) atc%num_bond_types
    if (atc%num_bond_types > 0) write(fu) atc%bond_styles, atc%bond_params
    write(fu) atc%num_bonds
    if (atc%num_bonds > 0) write(fu) atc%bonds

    write(fu) atc%num_angle_types
    if (atc%num_angle_types > 0) write(fu) atc%angle_styles, atc%angle_params
    write(fu) atc%num_angles
    if (atc%num_angles > 0) write(fu) atc%angles

    write(fu) atc%num_dihedral_types
    if (atc%num_dihedral_types > 0) &
        write(fu) atc%dihedral_styles, atc%dihedral_params
    write(fu) atc%num_dihedrals
    if (atc%num_dihedrals > 0) write(fu) atc%dihedrals

    write(fu) atc%num_branches
    if (atc%num_branches > 0) write(fu) atc%branches

    write(fu) atc%num_molecule_types
    write(fu) atc%molecule_names, atc%molecule_pop

    write(fu) atc%num_molecules
    write(fu) atc%molecules

    write(fu) atc%num_tether_types
    if (atc%num_tether_types > 0) write(fu) atc%tether_styles, atc%tether_params
    write(fu) atc%num_tethers
    if (atc%num_tethers > 0) write(fu) atc%tethers, atc%tether_points

    write(fu) atc%num_vdw_types
    if (atc%num_vdw_types > 0) write(fu) atc%vdw_pairs, atc%vdw_styles, &
                                atc%vdw_params

    write(fu) atc%num_externals
    if (atc%num_externals > 0) write(fu) atc%external_styles, atc%external_params

    write(fu) atc%flow_style
    if (atc%flow_style /= 0) write(fu) atc%flow_params

    if ( (atc%num_tethers==0) .and. (simbox%imcon==0) ) write(fu) atc%molecule_com

    close(fu)

    end subroutine

!******************************************************************************

subroutine read_config(simbox, atc, fn)
    !! Read from CONFIG file

    type(smbx_t),     intent(out) :: simbox
    type(atmcfg_t),   intent(out) :: atc
    character(len=*), intent(in)  :: fn
    integer, parameter :: mxrdln = 128
    character(len=mxrdln) :: line
    character(len=:), allocatable :: word
    real(rp), dimension(3,3) :: mat
    integer :: i, it, npar, ibr, n, nt
    integer :: fu, ios

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '#', ios)
        if (ios /= 0) exit

        line = adjustl(line)

        if (str_startswith(line, 'SIMBOX')) then
            !Read box boundary condition
            call str_split(line, ' ', word)
            if (len_trim(line) == 0) then
                !Read simulation box size into lattice vectors defining the box
                read(fu,*) mat(:,1)
                read(fu,*) mat(:,2)
                read(fu,*) mat(:,3)
                read(fu,*)
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                n = str_to_i(line)
                call simbox%init(n)
                call simbox%set_basis(mat)
            else
                n = str_to_i(line)
                call simbox%init(n)
                !Read simulation box size into lattice vectors defining the box
                read(fu,*) mat(:,1)
                read(fu,*) mat(:,2)
                read(fu,*) mat(:,3)
                call simbox%set_basis(mat)
            end if
        end if

        if (str_startswith(line, 'ATOM_TYPES')) then
            call str_split(line, ' ', word)
            nt = str_to_i(line)
            atc%num_atom_types = nt
            allocate(atc%atom_names (nt))
            allocate(atc%atom_styles(nt))
            allocate(atc%atom_mass  (nt))

            do it = 1, nt
                read(fu,*) atc%atom_names(it), atc%atom_styles(it), &
                            atc%atom_mass(it)
            end do
        end if

        if (str_startswith(line, 'ATOMS')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_atoms = n

            allocate(atc%atoms (n))
            allocate(atc%charge (n))
            allocate(atc%coordinates(3,n))

            do i = 1, n
                read(fu,*) atc%atoms(i), atc%charge(i), atc%coordinates(:,i)
            end do
        end if

        if (str_startswith(line, 'BOND_TYPES')) then
            call str_split(line, ' ', word)
            nt = str_to_i(line)
            atc%num_bond_types = nt

            allocate(atc%bond_styles(nt))
            allocate(atc%bond_params(mxparam,nt))
            atc%bond_params = 0.0_rp

            do it = 1, nt
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                atc%bond_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) atc%bond_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'BONDS')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_bonds = n
            allocate(atc%bonds(3,n))
            do i = 1, n
                read(fu,*) atc%bonds(:,i)
            end do
        end if

        if (str_startswith(line, 'ANGLE_TYPES')) then
            call str_split(line, ' ', word)
            nt = str_to_i(line)
            atc%num_angle_types = nt

            allocate(atc%angle_styles(nt))
            allocate(atc%angle_params(mxparam, nt))
            atc%angle_params = 0.0_rp

            do it = 1, nt
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                atc%angle_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) atc%angle_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'ANGLES')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_angles = n
            allocate(atc%angles(4,n))
            do i = 1, n
                read(fu,*) atc%angles(:,i)
            end do
        end if

        if (str_startswith(line, 'DIHEDRAL_TYPES')) then
            call str_split(line, ' ', word)
            nt = str_to_i(line)
            atc%num_dihedral_types = nt

            allocate(atc%dihedral_styles(nt))
            allocate(atc%dihedral_params(mxparam, nt))
            atc%dihedral_params = 0.0_rp

            do it = 1, nt
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                atc%dihedral_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) atc%dihedral_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'DIHEDRALS')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_dihedrals = n
            allocate(atc%dihedrals(5,n))
            do i = 1, n
                read(fu,*) atc%dihedrals(:,i)
            end do
        end if

        if (str_startswith(line, 'BRANCHES')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_branches = n
            allocate(atc%branches(3,n))
            do ibr = 1, n
                read(fu,*) atc%branches(:,ibr)
            end do
        end if

        if (str_startswith(line, 'MOLECULE_TYPES')) then
            call str_split(line, ' ', word)
            nt = str_to_i(line)
            atc%num_molecule_types = nt

            allocate(atc%molecule_names(nt))
            allocate(atc%molecule_pop(nt))

            do it = 1, nt
                read(fu,*) atc%molecule_names(it), atc%molecule_pop(it)
            end do
        end if

        if (str_startswith(line, 'MOLECULES')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_molecules = n
            allocate(atc%molecules(9,n))
            do i = 1, n
                read(fu,*) atc%molecules(:,i)
            end do
        end if

        if (str_startswith(line, 'TETHER_TYPES')) then
            call str_split(line, ' ', word)
            nt = str_to_i(line)
            atc%num_tether_types = nt

            allocate(atc%tether_styles(nt))
            allocate(atc%tether_params(mxparam, nt))
            atc%tether_params = 0.0_rp

            do it = 1, nt
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                atc%tether_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) atc%tether_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'TETHERS')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_tethers = n
            allocate(atc%tethers(2,n))
            allocate(atc%tether_points(3,n))
            do i = 1, n
                read(fu,*) atc%tethers(:,i), atc%tether_points(:,i)
            end do
        end if

        if (str_startswith(line, 'VDW')) then
            call str_split(line, ' ', word)
            nt = str_to_i(line)
            atc%num_vdw_types = nt

            allocate(atc%vdw_styles(nt))
            allocate(atc%vdw_pairs(2,nt))
            allocate(atc%vdw_params(mxparam,nt))
            atc%vdw_params = 0.0_rp

            do it = 1, nt
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                atc%vdw_pairs(1,it) = str_to_i(word)
                call str_split(line, ' ', word)
                atc%vdw_pairs(2,it) = str_to_i(word)
                call str_split(line, ' ', word)
                atc%vdw_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) atc%vdw_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'EXTERNAL')) then
            call str_split(line, ' ', word)
            n = str_to_i(line)
            atc%num_externals = n

            allocate(atc%external_styles(n))
            allocate(atc%external_params(mxparam,n))
            atc%external_params = 0.0_rp

            do it = 1, n
                call readline(fu, line, '#', ios)
                call str_split(line, ' ', word)
                atc%external_styles(it) = str_to_i(word)
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) atc%external_params(1:npar,it)
                end if
            end do
        end if

        if (str_startswith(line, 'FLOW_FIELD')) then
            allocate(atc%flow_params(mxparam))
            atc%flow_params = 0.0_rp
            call readline(fu, line, '#', ios)
            call str_split(line, ' ', word)
            atc%flow_style = str_to_i(word)
            if (atc%flow_style > 0) then
                call str_split(line, ' ', word)
                npar = str_to_i(word)
                if (npar > 0) then
                    read(line, *) atc%flow_params(1:npar)
                end if
            end if
        end if
    end do

    close(fu)

    end subroutine

!*******************************************************************************

subroutine write_config(simbox, atc, fn, title)
    !! Write to cfg file

    type(smbx_t),     intent(in) :: simbox
    type(atmcfg_t),   intent(in) :: atc
    character(len=*), intent(in) :: fn
    character(len=*), intent(in) :: title
    integer :: fu
    integer :: i, it

    open(newunit=fu, file=fn, action='write', status='replace')

    write(fu, '(a)') '#'//trim(adjustl(title))
    write(fu, '(a)') 'version 2.0'

    write(fu, *)
    write(fu, '(a,2x,i0)') 'SIMBOX', simbox%imcon
    write(fu, '(3(g0.8,2x))') simbox%basis(:,1)
    write(fu, '(3(g0.8,2x))') simbox%basis(:,2)
    write(fu, '(3(g0.8,2x))') simbox%basis(:,3)

    write(fu, *)
    write(fu, '(a,2x,i0)') 'ATOM_TYPES', atc%num_atom_types
    do it = 1, atc%num_atom_types
        write(fu, '(a,2x,i0,2x,g0.6)') trim(adjustl(atc%atom_names(it))), &
            atc%atom_styles(it), atc%atom_mass(it)
    end do

    write(fu, *)
    write(fu, '(a,2x,i0)') 'ATOMS', atc%num_atoms
    do i = 1, atc%num_atoms
        write(fu, '(i0,2x,*(g0.14,2x))') atc%atoms(i), atc%charge(i), &
            atc%coordinates(:,i)
    end do

    if (atc%num_bonds > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'BOND_TYPES', atc%num_bond_types
        do it = 1, atc%num_bond_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') atc%bond_styles(it), &
                size(atc%bond_params,1), atc%bond_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'BONDS', atc%num_bonds
        do i = 1, atc%num_bonds
            write(fu, '(*(i0,2x))') atc%bonds(:,i)
        end do
    end if

    if (atc%num_angles > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'ANGLE_TYPES', atc%num_angle_types
        do it = 1, atc%num_angle_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') atc%angle_styles(it), &
                size(atc%angle_params,1), atc%angle_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'ANGLES', atc%num_angles
        do i = 1, atc%num_angles
            write(fu, '(*(i0,2x))') atc%angles(:,i)
        end do
    end if

    if (atc%num_dihedrals > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'DIHEDRAL_TYPES', atc%num_dihedral_types
        do it = 1, atc%num_dihedral_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') atc%dihedral_styles(it), &
                size(atc%dihedral_params,1), atc%dihedral_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'DIHEDRALS', atc%num_dihedrals
        do i = 1, atc%num_dihedrals
            write(fu, '(*(i0,2x))') atc%dihedrals(:,i)
        end do
    end if

    if (atc%num_branches > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'BRANCHES', atc%num_branches
        do i = 1, atc%num_branches
            write(fu, '(*(i0,2x))') atc%branches(:,i)
        end do
    end if

    write(fu, *)
    write(fu, '(a,2x,i0)') 'MOLECULE_TYPES', atc%num_molecule_types
    do it = 1, atc%num_molecule_types
        write(fu, '(a,2x,i0)') trim(adjustl(atc%molecule_names(it))), &
            atc%molecule_pop(it)
    end do

    write(fu, *)
    write(fu, '(a,2x,i0)') 'MOLECULES', atc%num_molecules
    do i = 1, atc%num_molecules
        write(fu, '(*(i0,2x))') atc%molecules(:,i)
    end do

    if (atc%num_tethers > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'TETHER_TYPES', atc%num_tether_types
        do it = 1, atc%num_tether_types
            write(fu, '(i0,2x,i0,2x,*(g0.6,2x))') atc%tether_styles(it), &
                size(atc%tether_params,1), atc%tether_params(:,it)
        end do

        write(fu, *)
        write(fu, '(a,2x,i0)') 'TETHERS', atc%num_tethers
        do i = 1, atc%num_tethers
            write(fu, '(2(i0,2x),3(g0.6,2x))') atc%tethers(:,i), &
                atc%tether_points(:,i)
        end do
    end if

    if (atc%num_vdw_types > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'VDW', atc%num_vdw_types
        do it = 1, atc%num_vdw_types
            write(fu,'(2(i0,2x),i0,2x,i0,2x,*(g0.6,2x))') atc%vdw_pairs(:,it), &
                atc%vdw_styles(it), size(atc%vdw_params,1), atc%vdw_params(:,it)
        end do
    end if

    if (atc%num_externals > 0) then
        write(fu, *)
        write(fu, '(a,2x,i0)') 'EXTERNAL', atc%num_externals
        do it = 1, atc%num_externals
            write(fu,'(i0,2x,i0,2x,*(g0.6,2x))') atc%external_styles(it), &
                size(atc%external_params,1), atc%external_params(:,it)
        end do
    end if

    write(fu, *)
    write(fu, '(a)') 'FLOW_FIELD'
    write(fu,'(i0,2x,i0,2x,*(g0.6,2x))', advance='no') atc%flow_style
    if (atc%flow_style > 0) then 
        write(fu,'(i0,2x,*(g0.6,2x))') size(atc%flow_params), atc%flow_params
    else
        write(fu,*)
    end if

    close(fu)

    end subroutine

!******************************************************************************

subroutine write_ldf(simbox, atc, fn_ld, title)
    !! Write to a LAMMPS data file.

    type(smbx_t),     intent(in) :: simbox
    type(atmcfg_t),   intent(in) :: atc
    character(len=*), intent(in) :: fn_ld
        !! Name of the file
    character(len=*), intent(in) :: title
        !! Title of the configuation
    real(rp), dimension(3) :: tilt_factors
    integer :: fu_ld
    integer :: cntr_atm, iatm_beg, natm
    integer :: i, iatm, imol

    open(newunit=fu_ld, file=fn_ld, action='write')

    !Header
    write(fu_ld,'(a)') '#'//trim(adjustl(title))

    write(fu_ld,'(i0,2x,a)') atc%num_atoms, 'atoms'
    write(fu_ld,'(i0,2x,a)') atc%num_atom_types, 'atom types'

    if (atc%num_bonds > 0) then
        write(fu_ld,'(i0,2x,a)') atc%num_bonds, 'bonds'
        write(fu_ld,'(i0,2x,a)') atc%num_bond_types, 'bond types'
    end if

    if (atc%num_angles > 0) then
        write(fu_ld,'(i0,2x,a)') atc%num_angles, 'angles'
        write(fu_ld,'(i0,2x,a)') atc%num_angle_types, 'angle types'
    end if

    if (atc%num_dihedrals > 0) then
        write(fu_ld,'(i0,2x,a)') atc%num_dihedrals, 'dihedrals'
        write(fu_ld,'(i0,2x,a)') atc%num_dihedral_types, 'dihedral types'
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
    do imol = 1, atc%num_molecules
        natm = atc%molecules(2,imol)
        iatm_beg = atc%molecules(3,imol)
        do i = 1, natm
            iatm = iatm_beg + i -1
            write(fu_ld,'(i0,2x,i0,2x,i0,2x,4(g0.8,2x))') cntr_atm, imol, &
                atc%atoms(iatm), atc%charge(iatm), atc%coordinates(:,iatm)
            cntr_atm = cntr_atm + 1
        end do
    end do

    !Body: Bonds
    if (atc%num_bonds > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Bonds'
        write(fu_ld,*)
        do i = 1, atc%num_bonds
            write(fu_ld,'(i0,2x,3(i0,2x))') i, atc%bonds(:,i)
        end do
    end if

    !Body: Angles
    if (atc%num_angles > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Angles'
        write(fu_ld,*)
        do i = 1, atc%num_angles
            write(fu_ld,'(i0,2x,4(i0,2x))') i, atc%angles(:,i)
        end do
    end if

    !Body: Dihedrals
    if (atc%num_dihedrals > 0) then
        write(fu_ld,*)
        write(fu_ld,'(a)') 'Dihedrals'
        write(fu_ld,*)
        do i = 1, atc%num_dihedrals
            write(fu_ld,'(i0,2x,5(i0,2x))') i, atc%dihedrals(:,i)
        end do
    end if

    close(fu_ld)

    end subroutine

!*******************************************************************************

subroutine write_xyz(atc, fn_xyz, title)
    !! Write to an XYZ file.

    type(atmcfg_t),   intent(in) :: atc
    character(len=*), intent(in) :: fn_xyz
        !! Name of the XYZ file
    character(len=*), intent(in) :: title
        !! Title (for the configuration)
    integer :: iatm
    integer :: na
    integer :: fu_xyz

    na = atc%num_atoms

    open(newunit=fu_xyz, file=fn_xyz, action='write')

    write(fu_xyz, '(i0)') na
    write(fu_xyz, '(a)') title

    do iatm = 1, na
        write(fu_xyz,'(3(es24.15,2x))') atc%coordinates(:,iatm)
    end do

    close(fu_xyz)

    end subroutine

!*******************************************************************************

end module config_io_m
