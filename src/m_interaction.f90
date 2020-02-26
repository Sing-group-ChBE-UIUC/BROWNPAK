module m_interaction
!! Driver routines for force & energy calculation.

use omp_lib
use m_precision
use m_constants_math
use m_globals
use m_connectivity
use m_verlet
use m_cell_list
use m_ia_bond
use m_ia_angle
use m_ia_dihedral
use m_ia_vdw
use m_ia_tether
use m_ia_external

implicit none

private

public :: ia_setup, ia_finish, ia_calc_forces

contains

!******************************************************************************

subroutine ia_setup()
    !! Sets up parameters for potentials 

    !Initialize and build atom -> bond, atom -> atom and next nearest bonded
    !neighbor tables
    call atbo_build()

    !Debug statements
    !write(*,*) 'ATBO_TAB'
    !call atbo_tab%print()

    !Build excluded atoms table
    call exat_build()

    !Debug statements
    !write(*,*) 'EXAT_TAB'
    !call exat_tab%print()

    !Delete atom -> bond table, it is no longer needed.
    call atbo_tab%delete()

    !Check
    if (num_vdw_types == 0) lvdw = .false.

    if ((imcon==0) .and. use_cell_list) use_cell_list = .false.
    if (use_verlet_tab .and. use_cell_list) use_verlet_tab = .false.

    !Initialize verlet list
    if (use_verlet_tab) call verlet_init(rcutoff+tskin, tskin)

    !Initialize cell list
    if ( (sim_style == 0) .or. (sim_style == 1) ) then
        ! For relaxation/BD simulation the cell list is used for short range
        ! force calculation.
        if (use_cell_list) then
            call cl_init(num_atoms, rcutoff)
            call cl_set_cell_size(rcutoff)
            call cl_build_cell_nbrs()
        end if
    else if (sim_style == 2) then
        ! For MPCD simulation the cell list is always used for sorting. If
        ! use_cell_list == .true., the cell list will also be used for force
        ! calculation.
        call cl_init(num_atoms_tot, 1.0_rp) !Collision cell size is 1.0.
        if (use_cell_list) then
            call cl_set_cell_size(rcutoff) !Set cell size to global cutoff 
            call cl_build_cell_nbrs()
        end if
    end if

    if (num_bonds > 0) call ia_bond_setup()
    if (num_angles > 0) call ia_angle_setup()
    if (num_dihedrals > 0) call ia_dihedral_setup()
    if (num_tethers > 0) call ia_tether_setup()
    if (lvdw) call ia_vdw_setup()
    if (num_externals > 0) call ia_external_setup()

    end subroutine

!******************************************************************************

subroutine ia_finish()
    !! Releases memory allocated in `ia_setup`.

    call exat_tab%delete()
    if (use_verlet_tab) call verlet_delete()
    if (use_cell_list) call cl_delete()

    end subroutine

!******************************************************************************

subroutine ia_calc_forces(ierr)
    !! Calculates total forces and energies

    integer, intent(out) :: ierr

    ierr = 0

    !Zeroing out force, energies, & bond length
    forces = 0.0_rp

    energy_bond = 0.0_rp
    energy_angle = 0.0_rp
    energy_dihedral = 0.0_rp
    energy_tether = 0.0_rp
    energy_vdw = 0.0_rp
    energy_external = 0.0_rp
    energy_tot = 0.0_rp
    stress = 0.0_rp

    if (num_bonds > 0) then
        bndlen = 0.0_rp
        bndlen_min = huge(0.0_rp)
        bndlen_max = 0.0_rp
    end if

    !Calculation of bonded interactions
    if (num_bonds > 0) then
        call ia_add_bond_forces(ierr)
        if (ierr /= 0) return
    end if

    !Calculation of angular interactions
    if (num_angles > 0) call ia_add_angle_forces()

    !Calculation of dihedral interactions
    if (num_dihedrals > 0) call ia_add_dihedral_forces()

    !Calculation of tether interactions
    if (num_tethers > 0) then
        call ia_add_tether_forces(ierr)
        if (ierr /= 0) return
    end if

    !Calculation for pairwise interactions
    if (lvdw) then
        if (use_verlet_tab) then
            call ia_add_vdw_forces_vl(ierr)
        else if (use_cell_list) then
            call ia_add_vdw_forces_cl(ierr)
        else
            call ia_add_vdw_forces(ierr)
        end if
        if (ierr /= 0) return
    end if

    !Calculation of external interactions
    if (num_externals > 0) then
        call ia_add_external_forces(ierr)
        if (ierr /= 0) return
    end if

    !Calculate total energy
    energy_tot = energy_bond + energy_angle + energy_dihedral + energy_vdw &
        + energy_tether + energy_external

    !Update stress considering volume for finite concentration
    if (imcon /= 0) stress = stress/simbox%volume

    end subroutine

!********************************************************************************

subroutine ia_add_vdw_forces(ierr)
    !! Calculates force and energy due to all short-ranged non-bonded pairwise
    !! interactions and adds to `energy_vdw` & 'forces` in module `m_globals`.
    !! Uses direct N^2 calculation.

    integer, intent(out) :: ierr
    real(rp), dimension(3) :: ri, rj, rij, fi
    real(rp) :: rij_mag
    real(rp) :: qi, qj
    real(rp) :: enrg, frc
    integer :: i, j, at_i, at_j, typ

    ierr = 0

    do i = 1, (num_atoms-1)
        ri = coordinates(:,i)
        qi = charge(i)
        at_i = atoms(i)
        do j = i+1, num_atoms
            !Check if atom j is an excluded atom for atom i
            if ( exat_tab%is_in(i,j) ) cycle
            rj = coordinates(:,j)
            qj = charge(j)
            at_j = atoms(j)
            if (at_i < at_j) then
                typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
            else
                typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
            end if
            rij = rj - ri
            if (imcon /= 0) call simbox%get_image(rij)
            rij_mag = norm2(rij)
            call ia_get_vdw_force(rij_mag, qi, qj, typ, enrg, frc, ierr)
            if (ierr /= 0) return
            energy_vdw = energy_vdw + enrg
            fi = frc*rij/rij_mag

            !Update forces
            forces(:,i) = forces(:,i) + fi
            forces(:,j) = forces(:,j) - fi

            !Update stress
            if (update_stress /= 0) then
                stress(1,1) = stress(1,1) - rij(1)*fi(1)
                stress(2,1) = stress(2,1) - rij(2)*fi(1)
                stress(3,1) = stress(3,1) - rij(3)*fi(1)

                stress(1,2) = stress(1,2) - rij(1)*fi(2)
                stress(2,2) = stress(2,2) - rij(2)*fi(2)
                stress(3,2) = stress(3,2) - rij(3)*fi(2)

                stress(1,3) = stress(1,3) - rij(1)*fi(3)
                stress(2,3) = stress(2,3) - rij(2)*fi(3)
                stress(3,3) = stress(3,3) - rij(3)*fi(3)
            end if
        end do
    end do

    end subroutine

!********************************************************************************

subroutine ia_add_vdw_forces_vl(ierr)
    !! Calculates force and energy due to all short-ranged non-bonded pairwise
    !! interactions and adds to `energy_vdw` & 'forces` in module `m_globals`.
    !! Uses Verlet table.

    integer, intent(out) :: ierr
    integer, dimension(:), pointer :: nbrs => null()
    real(rp), dimension(3) :: ri, rj, rij, fi
    real(rp) :: rij_mag
    real(rp) :: qi, qj
    real(rp) :: enrg, frc
    integer :: i, j, k, at_i, at_j, typ

    ierr = 0

    call verlet_build()

    !call verlet_tab%print()

    do i = 1, (num_atoms-1)
        ri = coordinates(:,i)
        qi = charge(i)
        at_i = atoms(i)

        !Getting list of neighbors of particle i using pointer nbrs
        call verlet_tab%get_row(i, nbrs)

        do k = 1, size(nbrs)
            j = nbrs(k)
            !Check if atom j is an excluded atom for atom i
            if ( exat_tab%is_in(i,j) ) cycle

            rj = coordinates(:,j)
            qj = charge(j)
            at_j = atoms(j)
            if (at_i < at_j) then
                typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
            else
                typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
            end if
            rij = rj - ri
            if (imcon /= 0) call simbox%get_image(rij)
            rij_mag = norm2(rij)
            call ia_get_vdw_force(rij_mag, qi, qj, typ, enrg, frc, ierr)
            if (ierr /= 0) return

            energy_vdw = energy_vdw + enrg
            fi = frc*rij/rij_mag
            !Update forces
            forces(:,i) = forces(:,i) + fi
            forces(:,j) = forces(:,j) - fi

            !Update stress
            if (update_stress /= 0) then
                stress(1,1) = stress(1,1) - rij(1)*fi(1)
                stress(2,1) = stress(2,1) - rij(2)*fi(1)
                stress(3,1) = stress(3,1) - rij(3)*fi(1)

                stress(1,2) = stress(1,2) - rij(1)*fi(2)
                stress(2,2) = stress(2,2) - rij(2)*fi(2)
                stress(3,2) = stress(3,2) - rij(3)*fi(2)

                stress(1,3) = stress(1,3) - rij(1)*fi(3)
                stress(2,3) = stress(2,3) - rij(2)*fi(3)
                stress(3,3) = stress(3,3) - rij(3)*fi(3)
            end if
        end do
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_add_vdw_forces_cl(ierr)
    !! Calculates force and energy due to all short-ranged non-bonded pairwise
    !! interactions and adds to `energy_vdw` & 'forces` in module `m_globals`.
    !! Uses cell list.

    integer, intent(out) :: ierr
    integer, dimension(:), pointer :: aic => null()
    integer, dimension(:), pointer :: ainc => null()
    integer, dimension(:), pointer :: nbr_cells => null()
    real(rp), dimension(3) :: ri, rj, rij, fi
    real(rp), dimension(3,num_atoms) :: forces_temp
    real(rp) :: enrg_temp
    real(rp) :: rij_mag
    real(rp) :: qi, qj
    real(rp) :: enrg, frc
    integer :: num_cells
    integer :: icell, jcell, iatm, jatm
    integer :: i, j, k, at_i, at_j, typ
    logical :: is_parallel = .false.

    ierr = 0

    call cl_build(coordinates(:,1:num_atoms))
    num_cells = cl_get_num_cells()

    !$omp parallel &
    !$omp default(shared) &
    !$omp private(icell,aic,nbr_cells,i,iatm,ri,qi,at_i,j,jatm,rj,qj,at_j,typ,rij,rij_mag,frc,fi,forces_temp,enrg_temp,jcell,ainc) 

    forces_temp = 0.0_rp
    enrg_temp = 0.0_rp
    ! is_parallel = omp_in_parallel()
    ! if (nts <= 10 .and. is_parallel) then
    !     write(*,*) '1'
    ! end if
    !$omp do schedule(static)

    do icell = 0, (num_cells-1)
        call cl_get_contents(icell, aic)
        call cl_get_nbr_cells(icell, nbr_cells)

        !Interaction with particles within the cell
        do i = 1, size(aic) - 1
            iatm = aic(i)
            ri = coordinates(:,iatm)
            qi = charge(iatm)
            at_i = atoms(iatm)

            do j = i+1, size(aic)
                jatm = aic(j)
                !Check if atom jatm is an excluded atom for atom iatm
                if ( exat_tab%is_in(iatm,jatm) ) cycle

                rj = coordinates(:,jatm)
                qj = charge(jatm)
                at_j = atoms(jatm)
                if (at_i < at_j) then
                    typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
                else
                    typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
                end if
                rij = rj - ri
                if (imcon /= 0) call simbox%get_image(rij)
                rij_mag = norm2(rij)
                call ia_get_vdw_force(rij_mag, qi, qj, typ, enrg, frc, ierr)
                ! if (ierr /= 0) return

                enrg_temp = enrg_temp + enrg
                fi = frc*rij/rij_mag

                !Update forces
                forces_temp(:,iatm) = forces_temp(:,iatm) + fi
                forces_temp(:,jatm) = forces_temp(:,jatm) - fi
                ! !$omp atomic update
                !     forces(1,i) = forces(1,i) + fi(1)
                !     forces(2,i) = forces(2,i) + fi(2)
                !     forces(3,i) = forces(3,i) + fi(3)
                ! !$omp atomic update 
                !     forces(1,j) = forces(1,j) - fi(1)
                !     forces(2,j) = forces(2,j) - fi(2)
                !     forces(3,j) = forces(3,j) - fi(3)

                !Update stress
                if (update_stress /= 0) then
                    stress(1,1) = stress(1,1) - rij(1)*fi(1)
                    stress(2,1) = stress(2,1) - rij(2)*fi(1)
                    stress(3,1) = stress(3,1) - rij(3)*fi(1)

                    stress(1,2) = stress(1,2) - rij(1)*fi(2)
                    stress(2,2) = stress(2,2) - rij(2)*fi(2)
                    stress(3,2) = stress(3,2) - rij(3)*fi(2)

                    stress(1,3) = stress(1,3) - rij(1)*fi(3)
                    stress(2,3) = stress(2,3) - rij(2)*fi(3)
                    stress(3,3) = stress(3,3) - rij(3)*fi(3)
                end if
            end do
        end do

        !Interaction with particles belonging to neighboring cells
        do i = 1, size(aic)
            iatm = aic(i)
            ri = coordinates(:,iatm)
            qi = charge(iatm)
            at_i = atoms(iatm)

            do j = 1, size(nbr_cells)
                jcell = nbr_cells(j)
                call cl_get_contents(jcell, ainc)

                do k = 1, size(ainc)
                    jatm = ainc(k)
                    !Check if atom jatm is an excluded atom for atom iatm
                    if ( exat_tab%is_in(iatm,jatm) ) cycle

                    rj = coordinates(:,jatm)
                    qj = charge(jatm)
                    at_j = atoms(jatm)

                    if (at_i < at_j) then
                        typ = at_j + (2*num_atom_types-at_i)*(at_i-1)/2
                    else
                        typ = at_i + (2*num_atom_types-at_j)*(at_j-1)/2
                    end if
                    rij = rj - ri
                    if (imcon /= 0) call simbox%get_image(rij)
                    rij_mag = norm2(rij)
                    call ia_get_vdw_force(rij_mag, qi, qj, typ, enrg, frc, ierr)
                    ! if (ierr /= 0) return

                    enrg_temp = enrg_temp + enrg
                    fi = frc*rij/rij_mag

                    !Update forces
                    forces_temp(:,iatm) = forces_temp(:,iatm) + fi
                    forces_temp(:,jatm) = forces_temp(:,jatm) - fi

                    !Update stress
                    if (update_stress /= 0) then
                        stress(1,1) = stress(1,1) - rij(1)*fi(1)
                        stress(2,1) = stress(2,1) - rij(2)*fi(1)
                        stress(3,1) = stress(3,1) - rij(3)*fi(1)

                        stress(1,2) = stress(1,2) - rij(1)*fi(2)
                        stress(2,2) = stress(2,2) - rij(2)*fi(2)
                        stress(3,2) = stress(3,2) - rij(3)*fi(2)

                        stress(1,3) = stress(1,3) - rij(1)*fi(3)
                        stress(2,3) = stress(2,3) - rij(2)*fi(3)
                        stress(3,3) = stress(3,3) - rij(3)*fi(3)
                    end if
                end do

            end do

        end do

    end do
    !$omp end do

    !$omp critical
    forces = forces + forces_temp

    energy_vdw = energy_vdw + enrg_temp

    !$omp end critical

    !$omp end parallel

    end subroutine

!*******************************************************************************

subroutine ia_add_bond_forces(ierr)
    !! Calculates forces & energy due to all bonds. Will add to 
    !! `energy_bond` & and 'forces` in module `m_globals`.

    integer, intent(out) :: ierr
    real(rp), dimension(3) :: ri, rj, rij, fi
    real(rp) :: rij_mag
    real(rp) :: enrg, frc
    integer  :: ibnd, typ
    integer  :: i, j

    ierr = 0

    ! abandoned paralization: not acceleration
    ! $omp parallel &
    ! $omp default(private) &
    ! $omp shared (forces,stress,ierr,imcon,energy_bond,bndlen_min,bndlen_max) 
    ! $omp do
    do ibnd = 1, num_bonds
        typ = bonds(1,ibnd)
        i = bonds(2,ibnd)
        j = bonds(3,ibnd)
        ri = coordinates(:,i)
        rj = coordinates(:,j)
        rij = rj - ri
        if (imcon /= 0) call simbox%get_image(rij)
        rij_mag = norm2(rij)

        bndlen = bndlen + rij_mag
        ! $omp atomic update
        bndlen_min = min(bndlen_min, rij_mag)
        bndlen_max = max(bndlen_max, rij_mag)

        call ia_get_bond_force(rij_mag, typ, enrg, frc, ierr)
        if (ierr /= 0) return
        ! this should be commented out because Openmp does not
        ! support jumping in/out of parallel regions
        fi = frc*rij/rij_mag 
        energy_bond = energy_bond + enrg
        ! $omp atomic update

        forces(:,i) = forces(:,i) + fi
        forces(:,j) = forces(:,j) - fi

        !Update stress
        ! $omp atomic update             
        if (update_stress /= 0) then
            stress(1,1) = stress(1,1) - rij(1)*fi(1)
            stress(2,1) = stress(2,1) - rij(2)*fi(1)
            stress(3,1) = stress(3,1) - rij(3)*fi(1)

            stress(1,2) = stress(1,2) - rij(1)*fi(2)
            stress(2,2) = stress(2,2) - rij(2)*fi(2)
            stress(3,2) = stress(3,2) - rij(3)*fi(2)

            stress(1,3) = stress(1,3) - rij(1)*fi(3)
            stress(2,3) = stress(2,3) - rij(2)*fi(3)
            stress(3,3) = stress(3,3) - rij(3)*fi(3)
        end if
    end do

    ! $omp end do
    ! $omp end parallel
    bndlen = bndlen/num_bonds

    end subroutine
    
!********************************************************************************

subroutine ia_add_angle_forces()
    !! Calculates forces & energy due to all angles. Will add to 
    !! `energy_angle` & 'forces` in module `m_globals`.
     
    real(rp), dimension(3) :: rim1, ri, rip1, q1, q2
    real(rp), dimension(3) :: fim1, fi, fip1
    real(rp) :: enrg
    integer  :: iang
    integer  :: typ
    integer  :: i, im1, ip1

    do iang = 1, num_angles
        typ = angles(1,iang)
        im1 = angles(2,iang)
        i   = angles(3,iang)
        ip1 = angles(4,iang)

        rim1 = coordinates(:,im1)
        ri   = coordinates(:,i)
        rip1 = coordinates(:,ip1)
        q1 = ri - rim1; q2 = rip1 - ri
        if (imcon /= 0) call simbox%get_image(q1)
        if (imcon /= 0) call simbox%get_image(q2)

        call ia_get_angle_force(q1, q2, typ, enrg, fim1, fi, fip1)
        energy_angle = energy_angle + enrg

        !Update forces
        forces(:,im1) = forces(:,im1) + fim1
        forces(:,i)   = forces(:,i)   + fi
        forces(:,ip1) = forces(:,ip1) + fip1

        !Update stress
        if (update_stress /= 0) then
            stress(1,1) = stress(1,1) - q1(1)*fim1(1) + q2(1)*fip1(1)
            stress(2,1) = stress(2,1) - q1(2)*fim1(1) + q2(2)*fip1(1)
            stress(3,1) = stress(3,1) - q1(3)*fim1(1) + q2(3)*fip1(1)

            stress(1,2) = stress(1,2) - q1(1)*fim1(2) + q2(1)*fip1(2)
            stress(2,2) = stress(2,2) - q1(2)*fim1(2) + q2(2)*fip1(2)
            stress(3,2) = stress(3,2) - q1(3)*fim1(2) + q2(3)*fip1(2)

            stress(1,3) = stress(1,3) - q1(1)*fim1(3) + q2(1)*fip1(3)
            stress(2,3) = stress(2,3) - q1(2)*fim1(3) + q2(2)*fip1(3)
            stress(3,3) = stress(3,3) - q1(3)*fim1(3) + q2(3)*fip1(3)
        end if
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_add_dihedral_forces()
    !! Calculates forces & energy due to all dihedrals. Will add to 
    !! `energy_dihedral` & 'forces` in module `m_globals`.

    real(rp), dimension(3) :: ri, rj, rk, rl
    real(rp), dimension(3) :: q1, q2, q3
    real(rp), dimension(3) :: fi, fj, fk, fl
    real(rp) :: enrg
    integer :: idhd
    integer :: typ
    integer :: i, j, k, l

    do idhd = 1, num_dihedrals
        typ = dihedrals(1,idhd)
        i = dihedrals(2, idhd)
        j = dihedrals(3, idhd)
        k = dihedrals(4, idhd)
        l = dihedrals(5, idhd)

        ri = coordinates(:,i)
        rj = coordinates(:,j)
        rk = coordinates(:,k)
        rl = coordinates(:,l)

        q1 = rj - ri; q2 = rk - rj; q3 = rl - rk
        if (imcon /= 0) call simbox%get_image(q1)
        if (imcon /= 0) call simbox%get_image(q2)
        if (imcon /= 0) call simbox%get_image(q3)

        call ia_get_dihedral_force(q1, q2, q3, typ, enrg, fi, fj, fk, fl)
        energy_dihedral = energy_dihedral + enrg
        forces(:,i) = forces(:,i) + fi
        forces(:,j) = forces(:,j) + fj
        forces(:,k) = forces(:,k) + fk
        forces(:,l) = forces(:,l) + fl

        !TODO: Add the contribution to stress here. 
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_add_tether_forces(ierr)
    !! Calculates forces & energy due to all tethers. Will add to `energy_tether` &
    !! 'forces` in module `m_globals`. Tether forces cannot be subject to
    !! periodic boundary conditions.

    integer, intent(out) :: ierr
    real(rp), dimension(3) :: r, tp, q, fj
    real(rp) :: qmag
    real(rp) :: enrg, frc
    integer :: iteth, typ, teth_iatm

    ierr = 0

    do iteth = 1, num_tethers
        typ = tethers(1,iteth)
        teth_iatm = tethers(2,iteth) !Index of the tethered atom
        tp = tether_points(:,iteth)
        r = coordinates(:,teth_iatm)
        q = r - tp
        qmag = norm2(q)
        call ia_get_tether_force(qmag, typ, enrg, frc, ierr)
        if (ierr /= 0) return

        fj = -frc*q/qmag
        energy_tether = energy_tether + enrg
        forces(:,teth_iatm) = forces(:,teth_iatm) + fj

        !Update stress
        !Sign flipped since fj, not fi is involved
        if (update_stress /= 0) then
            stress(1,1) = stress(1,1) + q(1)*fj(1)
            stress(2,1) = stress(2,1) + q(2)*fj(1)
            stress(3,1) = stress(3,1) + q(3)*fj(1)

            stress(1,2) = stress(1,2) + q(1)*fj(2)
            stress(2,2) = stress(2,2) + q(2)*fj(2)
            stress(3,2) = stress(3,2) + q(3)*fj(2)

            stress(1,3) = stress(1,3) + q(1)*fj(3)
            stress(2,3) = stress(2,3) + q(2)*fj(3)
            stress(3,3) = stress(3,3) + q(3)*fj(3)
        end if
    end do

    end subroutine
    
!******************************************************************************

end module m_interaction
