module m_ia_bond
!! This module contains routines to evaluate bond potentials and their
!! derivative.
!! 
!! The following styles are available:  
!!
!! * Style 0. None (only topology)
!! * Style 1. Harmonic. See [[bond_harm_set]].
!! * Style 2. FENE. See [[bond_fene_set]].
!! * Style 3. Kremer-Grest. See [[bond_kg_set]]. 
!! * Style 4. Marko-Siggia. See [[bond_ms_set]].  

use m_precision
use m_constants_math
use m_strings
use m_globals
use m_logger

implicit none

private

public :: ia_bond_setup, ia_get_bond_force

contains

!******************************************************************************

subroutine ia_bond_setup()
    !! Sets up parameters for bond potentials 

    integer :: i
    integer :: sty

    !Set bond interactions
    do i = 1, num_bond_types
        sty = bond_styles(i)
        select case(sty)
        case(1)
            call bond_harm_set(bond_params(:,i))
        case(2)
            call bond_fene_set(bond_params(:,i))
        case(3)
            call bond_kg_set(bond_params(:,i))
        case(4)
            call bond_ms_set(bond_params(:,i))
        case default
            continue
        end select
    end do

    end subroutine

!******************************************************************************

subroutine ia_get_bond_force(rij_mag, bnd_typ, enrg, frc, ierr)
    !! Calculates the energy & its derivative due to a bond.

    real(rp), intent(in) :: rij_mag
        !! Distance between bonded atoms
    integer, intent(in) :: bnd_typ
        !! Type of the bond
    real(rp), intent(out) :: enrg
        !! Bond energy
    real(rp), intent(out) :: frc
        !! Derivative of the potential. This is the magnitude of the force due
        !! to this potential.
    integer, intent(out)  :: ierr
        !! Error flag
    integer :: bnd_styl

    ierr = 0; enrg = 0.0_rp; frc = 0.0_rp
    bnd_styl = bond_styles(bnd_typ)

    select case (bnd_styl)
    case (1)
        call bond_harm(rij_mag, bond_params(:,bnd_typ), enrg, frc) 
    case (2)
        call bond_fene(rij_mag, bond_params(:,bnd_typ), enrg, frc, ierr) 
    case (3)
        call bond_kg(rij_mag, bond_params(:,bnd_typ), enrg, frc, ierr) 
    case (4)
        call bond_ms(rij_mag, bond_params(:,bnd_typ), enrg, frc, ierr) 
    case default
        continue
    end select

    end subroutine

!********************************************************************************

subroutine bond_harm_set(params, k, r0)
    !! Setter for harmonic bond interaction.
    !!
    !!```
    !!  U = (1/2)*k*(r - r0)^2
    !!```
    !! User-set parameters:
    !!
    !! * params(1) = `k` (spring constant)
    !! * params(2) = `r0` (equilibrium distance)

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: k
    real(rp), intent(in), optional :: r0

    if (present(k)) params(1) = k
    if (present(r0)) params(2) = r0

    end subroutine

!********************************************************************************

subroutine bond_harm(r, params, enrg, frc)
    !! Calculates energy & its derivative for harmonic bond. See [[bond_harm_set]].

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc

    enrg = 0.5_rp*params(1)*(r-params(2))*(r-params(2))
    frc = params(1)*(r-params(2))

    end subroutine

!********************************************************************************

subroutine bond_fene_set(params, k, rmax, r0)
    !! Setter for FENE bond.
    !!
    !!```
    !!  U = -0.5 k rmax^2 log [1 - ((r - r0)/rmax)^2]
    !!```
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `k`
    !! * params(2) = `rmax`
    !! * params(3) = `r0`
    !!
    !! @Note The bond cannot extend beyond (rmax+r0), where r0 is the
    !! equilibrium bond length. If r0 = 0, this reduces to the standard definition
    !! of FENE bonds.
    !!
    !! Internally stored parameters:
    !!
    !! * params(4) = `rmax^2`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: k
    real(rp), intent(in), optional :: rmax
    real(rp), intent(in), optional :: r0

    if (present(k))    params(1) = k
    if (present(rmax)) params(2) = rmax
    if (present(r0))   params(3) = r0

    params(4) = params(2)*params(2)

    end subroutine

!********************************************************************************

subroutine bond_fene(r, params, enrg, frc, ierr)
    !! Calculates energy & its derivative for FENE bond. See [[bond_fene_set]].
    !!
    !! If bond length exceeds maximum extensible spring length, an error will be
    !! reported by `ierr = 1`.

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    integer,  intent(out) :: ierr
    real(rp) :: k, rmax, r0
    real(rp) :: extn, extnsq, rmaxsq

    ierr = 0; enrg = 0.0_rp; frc = 0.0_rp

    k = params(1); rmax = params(2); r0 = params(3); rmaxsq = params(4)
    extn = r - r0; extnsq = extn*extn

    if ( r >= (rmax+r0) ) then
        ierr = 1
        call logger%log_msg('<bond_fene> bondlength too large')
        call logger%log_msg('<bond_fene> r = '//str_from_num(r))
        return
    else
        enrg = -0.5_rp*k*rmaxsq*log(1.0_rp - extnsq/rmaxsq)
        frc = k*extn/(1.0_rp - extnsq/rmaxsq)
    end if

    end subroutine

!******************************************************************************

subroutine bond_kg_set(params, k, rmax, eps, sigma)
    !! Setter for FENE bond interaction.
    !!
    !!```
    !!   V = 4*eps*[(r/sigma)^12 - (r/sigma)^6] + eps
    !!   W = -0.5 k rmax^2 log [1 - (r/rmax)^2]
    !!   U = W + V, r < 2^(1/6)*sigma
    !!       W, r >= 2^(1/6)*sigma
    !!```
    !! User-set parameters:
    !!
    !! * params(1) = `k`
    !! * params(2) = `rmax`
    !! * params(3) = `eps`
    !! * params(4) = `sigma`
    !!
    !! Internally stored parameters:
    !!
    !! * params(5) = `rmax^2`
    !! * params(6) = `2^(1/6)*sigma`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: k
    real(rp), intent(in), optional :: rmax
    real(rp), intent(in), optional :: eps
    real(rp), intent(in), optional :: sigma
    real(rp) :: k_, rmax_, eps_, sigma_

    if (present(k))     params(1) = k
    if (present(rmax))  params(2) = rmax
    if (present(eps))   params(3) = eps
    if (present(sigma)) params(4) = sigma

    k_ = params(1); rmax_ = params(2); eps_ = params(3); sigma_ = params(4)

    params(5) = rmax_**2
    params(6) = math_sxrt2*sigma_

    end subroutine

!********************************************************************************

subroutine bond_kg(r, params, enrg, frc, ierr)
    !! Calculates energy & its derivative for Kremer-Grest bond. See [[bond_kg_set]].
    !!
    !! If bond length exceeds maximum extensible spring length, an error will be
    !! reported as `ierr = 1`.

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    integer, intent(out)  :: ierr
    real(rp) :: k, rmax, eps, sigma
    real(rp) :: rmaxsq, rcut
    real(rp) :: rsq, sir, sir2, sir12, sir6

    ierr = 0; enrg = 0.0_rp; frc = 0.0_rp

    k = params(1); rmax = params(2); eps = params(3); sigma = params(4)
    rmaxsq = params(5); rcut = params(6)
    rsq = r*r

    if ( r >= rmax ) then
        ierr = 1
        call logger%log_msg('<bond_kg> r > rmax')
        call logger%log_msg('<bond_kg> r = '//str_from_num(r))
        return
    else if ( (r >= rcut) .and. (r < rmax) ) then
        enrg = -0.5_rp*k*rmaxsq*log(1.0_rp - rsq/rmaxsq)
        frc = k*r/(1.0_rp - rsq/rmaxsq)
    else
        sir = sigma/r
        sir2 = sir*sir; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
        enrg = -0.5_rp*k*rmaxsq*log(1.0_rp - rsq/rmaxsq) &
                + 4*eps*(sir12 - sir6) + eps
        frc = k*r/(1.0_rp - rsq/rmaxsq) - 24*eps*(2*sir12 - sir6)/r
    end if

    end subroutine

!******************************************************************************

subroutine bond_ms_set(params, lp, rmax)
    !! Setter for Marko-Siggia bond.
    !!
    !!```
    !!   U = [-(1/2)*rtilde^2 + 0.25/(1-rtilde)^2 + 0.25*rtilde]*(rmax/lp), r < rmax
    !!   where rtilde = r/rmax.
    !!```
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `lp` (persistence length)
    !! * params(2) = `rmax`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: lp
    real(rp), intent(in), optional :: rmax

    if (present(lp))    params(1) = lp
    if (present(rmax))  params(2) = rmax

    end subroutine

!********************************************************************************

subroutine bond_ms(r, params, enrg, frc, ierr)
    !! Evaluates the potential & its derivative for Marko-Siggia bond.
    !! See [[bond_ms_set]].
    !!
    !! If bond length exceeds maximum extensible spring length, an error will be
    !! reported as `ierr = 1`.

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    integer, intent(out)  :: ierr
    real(rp) :: lp, rmax, rtilde, rrtilde

    ierr = 0; enrg = 0.0_rp; frc = 0.0_rp
    lp = params(1); rmax = params(2)

    if ( r >= rmax ) then
        ierr = 1
        call logger%log_msg('<bond_ms> r > rmax')
        call logger%log_msg('<bond_ms> r = '//str_from_num(r))
        return
    else
        rtilde = r/rmax; rrtilde = 1.0_rp/(1.0_rp-rtilde)
        enrg = (-0.5*rtilde*rtilde + 0.25_rp*rrtilde*rrtilde + 0.25*rtilde)*(rmax/lp)
        frc = ( -rtilde - 0.25_rp*rrtilde*rrtilde + 0.25_rp )/lp
    end if

    end subroutine

!******************************************************************************

end module m_ia_bond
