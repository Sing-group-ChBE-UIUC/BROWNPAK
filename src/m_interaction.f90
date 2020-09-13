module m_interaction
!! Driver routines for force & energy calculation.

use m_precision
use m_constants_math
use m_globals
use m_pairtab
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
    !! Builds necessary neighbor tables and sets up parameters for potentials.

    if (num_vdw_types == 0) lvdw = .false.
    if (lvdw) then
        call pt_init()
        call ia_vdw_setup()
    end if

    if (num_bonds > 0) call ia_bond_setup()
    if (num_angles > 0) call ia_angle_setup()
    if (num_dihedrals > 0) call ia_dihedral_setup()
    if (num_tethers > 0) call ia_tether_setup()
    if (num_externals > 0) call ia_external_setup()

    end subroutine

!******************************************************************************

subroutine ia_finish()
    !! Cleanup routine for interaction calculation.

    call pt_finish()

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
        call ia_add_vdw_forces(ierr)
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
    !! interactions based on `pair_tab` and adds to `energy_vdw` & `forces`
    !! in module `m_globals`.

    integer, intent(out) :: ierr
    integer, dimension(:), pointer :: pnbrs => null()
    real(rp), dimension(3) :: ri, rj, rij, fi
    real(rp) :: rij_mag
    real(rp) :: qi, qj
    real(rp) :: enrg, frc
    integer :: i, j, k, at_i, at_j, typ

    ierr = 0

    call pt_build()

    do i = 1, num_atoms
        ri = coordinates(:,i)
        qi = charge(i)
        at_i = atoms(i)

        !Getting list of neighbors of atom i using pointer pnbrs
        call pair_tab%get_row(i, pnbrs)

        do k = 1, size(pnbrs)
            j = pnbrs(k)
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
        bndlen_min = min(bndlen_min, rij_mag)
        bndlen_max = max(bndlen_max, rij_mag)

        call ia_get_bond_force(rij_mag, typ, enrg, frc, ierr)
        if (ierr /= 0) return
        energy_bond = energy_bond + enrg
        fi = frc*rij/rij_mag
        forces(:,i) = forces(:,i) + fi
        forces(:,j) = forces(:,j) - fi

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

end module m_interaction
