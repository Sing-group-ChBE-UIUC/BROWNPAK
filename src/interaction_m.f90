module interaction_m
!! Driver routines for force & energy calculation.

use constants_m
use table_m
use simbox_m
use pairtab_m
use ia_bond_m
use ia_angle_m
use ia_dihedral_m
use ia_vdw_m
use ia_tether_m
use ia_external_m
use control_m
use atmcfg_m
use stats_m

implicit none

private
public :: ia_setup, ia_finish, ia_calc_forces

type(itable_t) :: pair_tab
logical :: lvdw

contains

!******************************************************************************

subroutine ia_setup(cpar, simbox, atc)
    !! Builds necessary neighbor tables and sets up parameters for potentials.

    type(ctrlpar_t), intent(in) :: cpar
    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t),  intent(in out) :: atc

    lvdw = .true.
    if ( (.not. cpar%lvdw) .or. (atc%num_vdw_types==0) ) lvdw = .false.

    if (lvdw) then
        call pt_init(cpar%mth_ptgen, atc%num_atoms, cpar%excluded_atoms, &
                cpar%rcutoff, cpar%tskin, atc%bonds, simbox, pair_tab)

        call ia_vdw_setup(atc%num_vdw_types, atc%vdw_styles, atc%vdw_params)
    end if

    if (atc%num_bonds > 0) then
        call ia_bond_setup(atc%num_bond_types, atc%bond_styles, atc%bond_params)
    end if

    if (atc%num_angles > 0) then
        call ia_angle_setup(atc%num_angle_types, atc%angle_styles, &
            atc%angle_params)
    end if

    if (atc%num_dihedrals > 0) then
        call ia_dihedral_setup(atc%num_dihedral_types, atc%dihedral_styles, &
            atc%dihedral_params)
    end if

    if (atc%num_tethers > 0) then
        call ia_tether_setup(atc%num_tether_types, atc%tether_styles, &
            atc%tether_params)
    end if

    if (atc%num_externals > 0) then
        call ia_external_setup(atc%num_externals, atc%external_styles, &
            atc%external_params)
    end if

    end subroutine

!******************************************************************************

subroutine ia_finish()
    !! Cleanup routine for interaction calculation.

    call pt_delete(pair_tab)

    end subroutine

!******************************************************************************

subroutine ia_calc_forces(simbox, atc, ierr)
    !! Calculates total forces and energies

    type(smbx_t),       intent(in) :: simbox
    type(atmcfg_t), intent(in out) :: atc
    integer, intent(out) :: ierr

    ierr = 0

    !Zeroing out force, energies, & bond length
    atc%forces = 0.0_rp

    energy_bond = 0.0_rp
    energy_angle = 0.0_rp
    energy_dihedral = 0.0_rp
    energy_tether = 0.0_rp
    energy_vdw = 0.0_rp
    energy_external = 0.0_rp
    energy_tot = 0.0_rp
    stress = 0.0_rp

    !Calculation of bonded interactions
    if (atc%num_bonds > 0) then
        bndlen = 0.0_rp
        bndlen_min = huge(0.0_rp)
        bndlen_max = 0.0_rp
        call ia_add_bond_forces(simbox, atc, ierr)
        if (ierr /= 0) return
    end if

    !Calculation of angular interactions
    if (atc%num_angles > 0) call ia_add_angle_forces(simbox, atc)

    !Calculation of dihedral interactions
    if (atc%num_dihedrals > 0) call ia_add_dihedral_forces(simbox, atc)

    !Calculation of tether interactions
    if (atc%num_tethers > 0) then
        call ia_add_tether_forces(atc, ierr)
        if (ierr /= 0) return
    end if

    !Calculation for pairwise interactions
    if (lvdw) then
        call ia_add_vdw_forces(simbox, atc, ierr)
        if (ierr /= 0) return
    end if

    !Calculation of external interactions
    if (atc%num_externals > 0) then
        call ia_add_external_forces(atc%num_externals, atc%external_styles, &
            atc%external_params, atc%coordinates, energy_external,      &
            atc%forces, stress, ierr)
        if (ierr /= 0) return
    end if

    !Calculate total energy
    energy_tot = energy_bond + energy_angle + energy_dihedral + energy_vdw &
        + energy_tether + energy_external

    !Update stress considering volume for finite concentration
    if (simbox%imcon /= 0) stress = stress/simbox%volume

    end subroutine

!********************************************************************************

subroutine ia_add_vdw_forces(simbox, atc, ierr)
    !! Calculates force and energy due to all short-ranged non-bonded pairwise
    !! interactions based on `pair_tab`.

    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in out) :: atc
    integer, intent(out) :: ierr
    integer, dimension(:), pointer :: pnbrs => null()
    real(rp), dimension(3) :: ri, rj, rij, fi
    real(rp) :: rij_mag
    real(rp) :: qi, qj
    real(rp) :: enrg, frc
    integer :: i, j, k, at_i, at_j, typ

    ierr = 0

    call pt_build(simbox, atc%coordinates, pair_tab)

    do i = 1, atc%num_atoms
        ri = atc%coordinates(:,i)
        qi = atc%charge(i)
        at_i = atc%atoms(i)

        !Getting list of neighbors of atom i using pointer pnbrs
        call pair_tab%get_row(i, pnbrs)

        do k = 1, size(pnbrs)
            j = pnbrs(k)
            rj = atc%coordinates(:,j)
            qj = atc%charge(j)
            at_j = atc%atoms(j)
            if (at_i < at_j) then
                typ = at_j + (2*atc%num_atom_types-at_i)*(at_i-1)/2
            else
                typ = at_i + (2*atc%num_atom_types-at_j)*(at_j-1)/2
            end if
            rij = rj - ri
            if (simbox%imcon /= 0) call simbox%get_image(rij)
            rij_mag = norm2(rij)
            call ia_get_vdw_force(rij_mag, qi, qj, atc%vdw_styles(typ), &
                atc%vdw_params(:,typ), enrg, frc, ierr)
            if (ierr /= 0) return

            energy_vdw = energy_vdw + enrg
            fi = frc*rij/rij_mag
            !Update forces
            atc%forces(:,i) = atc%forces(:,i) + fi
            atc%forces(:,j) = atc%forces(:,j) - fi

            !Update stress
            stress(1,1) = stress(1,1) - rij(1)*fi(1)
            stress(2,1) = stress(2,1) - rij(2)*fi(1)
            stress(3,1) = stress(3,1) - rij(3)*fi(1)

            stress(1,2) = stress(1,2) - rij(1)*fi(2)
            stress(2,2) = stress(2,2) - rij(2)*fi(2)
            stress(3,2) = stress(3,2) - rij(3)*fi(2)

            stress(1,3) = stress(1,3) - rij(1)*fi(3)
            stress(2,3) = stress(2,3) - rij(2)*fi(3)
            stress(3,3) = stress(3,3) - rij(3)*fi(3)
        end do
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_add_bond_forces(simbox, atc, ierr)
    !! Calculates forces & energy due to all bonds. Will add to 
    !! `energy_bond` & and 'forces` in module `m_globals`.

    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in out) :: atc
    integer, intent(out) :: ierr
    real(rp), dimension(3) :: ri, rj, rij, fi
    real(rp) :: rij_mag
    real(rp) :: enrg, frc
    integer  :: ibnd, typ
    integer  :: i, j

    ierr = 0

    do ibnd = 1, atc%num_bonds
        typ = atc%bonds(1,ibnd)
        i = atc%bonds(2,ibnd)
        j = atc%bonds(3,ibnd)
        ri = atc%coordinates(:,i)
        rj = atc%coordinates(:,j)
        rij = rj - ri
        if (simbox%imcon /= 0) call simbox%get_image(rij)
        rij_mag = norm2(rij)

        bndlen = bndlen + rij_mag
        bndlen_min = min(bndlen_min, rij_mag)
        bndlen_max = max(bndlen_max, rij_mag)

        call ia_get_bond_force(rij_mag, atc%bond_styles(typ), &
                atc%bond_params(:,typ), enrg, frc, ierr)
        if (ierr /= 0) return
        energy_bond = energy_bond + enrg
        fi = frc*rij/rij_mag
        atc%forces(:,i) = atc%forces(:,i) + fi
        atc%forces(:,j) = atc%forces(:,j) - fi

        !Update stress
        stress(1,1) = stress(1,1) - rij(1)*fi(1)
        stress(2,1) = stress(2,1) - rij(2)*fi(1)
        stress(3,1) = stress(3,1) - rij(3)*fi(1)

        stress(1,2) = stress(1,2) - rij(1)*fi(2)
        stress(2,2) = stress(2,2) - rij(2)*fi(2)
        stress(3,2) = stress(3,2) - rij(3)*fi(2)

        stress(1,3) = stress(1,3) - rij(1)*fi(3)
        stress(2,3) = stress(2,3) - rij(2)*fi(3)
        stress(3,3) = stress(3,3) - rij(3)*fi(3)
    end do

    bndlen = bndlen/atc%num_bonds

    end subroutine
    
!********************************************************************************

subroutine ia_add_angle_forces(simbox, atc)
    !! Calculates forces & energy due to all angles. Will add to 
    !! `energy_angle` & 'forces`.
     
    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in out) :: atc
    real(rp), dimension(3) :: rim1, ri, rip1, q1, q2
    real(rp), dimension(3) :: fim1, fi, fip1
    real(rp) :: enrg
    integer  :: iang
    integer  :: typ
    integer  :: i, im1, ip1

    do iang = 1, atc%num_angles
        typ = atc%angles(1,iang)
        im1 = atc%angles(2,iang)
        i   = atc%angles(3,iang)
        ip1 = atc%angles(4,iang)

        rim1 = atc%coordinates(:,im1)
        ri   = atc%coordinates(:,i)
        rip1 = atc%coordinates(:,ip1)
        q1 = ri - rim1; q2 = rip1 - ri
        if (simbox%imcon /= 0) call simbox%get_image(q1)
        if (simbox%imcon /= 0) call simbox%get_image(q2)

        call ia_get_angle_force(q1, q2, atc%angle_styles(typ), &
            atc%angle_params(:,typ), enrg, fim1, fi, fip1)
        energy_angle = energy_angle + enrg

        !Update forces
        atc%forces(:,im1) = atc%forces(:,im1) + fim1
        atc%forces(:,i)   = atc%forces(:,i)   + fi
        atc%forces(:,ip1) = atc%forces(:,ip1) + fip1

        !Update stress
        stress(1,1) = stress(1,1) - q1(1)*fim1(1) + q2(1)*fip1(1)
        stress(2,1) = stress(2,1) - q1(2)*fim1(1) + q2(2)*fip1(1)
        stress(3,1) = stress(3,1) - q1(3)*fim1(1) + q2(3)*fip1(1)

        stress(1,2) = stress(1,2) - q1(1)*fim1(2) + q2(1)*fip1(2)
        stress(2,2) = stress(2,2) - q1(2)*fim1(2) + q2(2)*fip1(2)
        stress(3,2) = stress(3,2) - q1(3)*fim1(2) + q2(3)*fip1(2)

        stress(1,3) = stress(1,3) - q1(1)*fim1(3) + q2(1)*fip1(3)
        stress(2,3) = stress(2,3) - q1(2)*fim1(3) + q2(2)*fip1(3)
        stress(3,3) = stress(3,3) - q1(3)*fim1(3) + q2(3)*fip1(3)
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_add_dihedral_forces(simbox, atc)
    !! Calculates forces & energy due to all dihedrals. Will add to 
    !! `energy_dihedral` & 'forces`.

    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in out) :: atc
    real(rp), dimension(3) :: ri, rj, rk, rl
    real(rp), dimension(3) :: q1, q2, q3
    real(rp), dimension(3) :: fi, fj, fk, fl
    real(rp) :: enrg
    integer :: idhd
    integer :: typ
    integer :: i, j, k, l

    do idhd = 1, atc%num_dihedrals
        typ = atc%dihedrals(1,idhd)
        i   = atc%dihedrals(2, idhd)
        j   = atc%dihedrals(3, idhd)
        k   = atc%dihedrals(4, idhd)
        l   = atc%dihedrals(5, idhd)

        ri = atc%coordinates(:,i)
        rj = atc%coordinates(:,j)
        rk = atc%coordinates(:,k)
        rl = atc%coordinates(:,l)

        q1 = rj - ri; q2 = rk - rj; q3 = rl - rk
        if (simbox%imcon /= 0) call simbox%get_image(q1)
        if (simbox%imcon /= 0) call simbox%get_image(q2)
        if (simbox%imcon /= 0) call simbox%get_image(q3)

        call ia_get_dihedral_force(q1, q2, q3, atc%dihedral_styles(typ), &
            atc%dihedral_params(:,typ), enrg, fi, fj, fk, fl)
        energy_dihedral = energy_dihedral + enrg
        atc%forces(:,i) = atc%forces(:,i) + fi
        atc%forces(:,j) = atc%forces(:,j) + fj
        atc%forces(:,k) = atc%forces(:,k) + fk
        atc%forces(:,l) = atc%forces(:,l) + fl

        !TODO: Add the contribution to stress here. 
    end do

    end subroutine
    
!********************************************************************************

subroutine ia_add_tether_forces(atc, ierr)
    !! Calculates forces & energy due to all tethers. Will add to `energy_tether` &
    !! 'forces`. Tether forces are not subject to periodic boundary conditions.

    type(atmcfg_t), intent(in out) :: atc
    integer, intent(out) :: ierr
    real(rp), dimension(3) :: r, tp, q, fj
    real(rp) :: qmag
    real(rp) :: enrg, frc
    integer :: iteth, typ, teth_iatm

    ierr = 0

    do iteth = 1, atc%num_tethers
        typ = atc%tethers(1,iteth)
        teth_iatm = atc%tethers(2,iteth) !Index of the tethered atom
        tp = atc%tether_points(:,iteth)
        r = atc%coordinates(:,teth_iatm)
        q = r - tp
        qmag = norm2(q)
        call ia_get_tether_force(qmag, atc%tether_styles(typ), &
            atc%tether_params(:,typ), enrg, frc, ierr)
        if (ierr /= 0) return

        fj = -frc*q/qmag
        energy_tether = energy_tether + enrg
        atc%forces(:,teth_iatm) = atc%forces(:,teth_iatm) + fj

        !Update stress
        !Sign flipped since fj, not fi is involved
        stress(1,1) = stress(1,1) + q(1)*fj(1)
        stress(2,1) = stress(2,1) + q(2)*fj(1)
        stress(3,1) = stress(3,1) + q(3)*fj(1)

        stress(1,2) = stress(1,2) + q(1)*fj(2)
        stress(2,2) = stress(2,2) + q(2)*fj(2)
        stress(3,2) = stress(3,2) + q(3)*fj(2)

        stress(1,3) = stress(1,3) + q(1)*fj(3)
        stress(2,3) = stress(2,3) + q(2)*fj(3)
        stress(3,3) = stress(3,3) + q(3)*fj(3)
    end do

    end subroutine
    
!******************************************************************************

end module interaction_m
