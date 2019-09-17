module m_ia_vdw
!! Routines to evaulate pairwise potentials and their derivative.
!!
!! The following styles are available:
!!
!! * Style 1. 12-6 LJ. See [[vdw_lj_set]].
!! * Style 2. Gaussian. See [[vdw_gaussian_set]].
!! * Style 3. Cosine. See [[vdw_cosine_set]].
!! * Style 4. Screened Coulomb + LJ. See [[vdw_lj_coul_debye_set]].
!! * Style 5. Coulomb + LJ. See [[vdw_lj_coul_set]].
!! * Style 6. Standard DPD. See [[vdw_dpd_set]].
!! * Style 7. expn/rx. See [[vdw_expnrx_set]].

use m_precision
use m_constants_math
use m_globals
use m_logger, only: logger_init, logger => master_logger

implicit none

private

public :: ia_vdw_setup, ia_get_vdw_force

contains

!******************************************************************************

subroutine ia_vdw_setup()
    !! Sets up parameters for vdw potentials 

    integer :: i
    integer :: sty

    !Set vdw interactions
    do i = 1, num_vdw_types
        sty = vdw_styles(i)

        select case(sty)
        case(1)
            call vdw_lj_set(vdw_params(:,i))
        case(2)
            call vdw_gaussian_set(vdw_params(:,i))
        case(3)
            call vdw_cosine_set(vdw_params(:,i))
        case(4)
            call vdw_lj_coul_debye_set(vdw_params(:,i))
        case(5)
            call vdw_lj_coul_set(vdw_params(:,i))
        case(6)
            call vdw_dpd_set(vdw_params(:,i))
        case(7)
            call vdw_expnrx_set(vdw_params(:,i))
        case default
            continue
        end select
    end do

    end subroutine

!******************************************************************************

subroutine ia_get_vdw_force(rij_mag, qi, qj, typ, enrg, frc, ierr)
    !! Calculates the energy & its derivative due to a single interacting pair of atoms.

    real(rp), intent(in) :: rij_mag
        !! Distance between two atoms
    real(rp), intent(in) :: qi
        !! Charge on atom i
    real(rp), intent(in) :: qj
        !! Charge on atom j
    integer, intent(in) :: typ
        !! Type of vdw interaction
    real(rp), intent(out) :: enrg
        !! Energy
    real(rp), intent(out) :: frc
        !! Derivative of the potential. This is the magnitude of force due to
        !! this potential.
    integer, intent(out) :: ierr
        !! Error flag
    integer :: styl

    ierr = 0; enrg = 0.0_rp; frc = 0.0_rp

    styl = vdw_styles(typ)
    select case(styl)
    case(1)
        call vdw_lj(rij_mag, vdw_params(:,typ), enrg, frc)
    case(2)
        call vdw_gaussian(rij_mag, vdw_params(:,typ), enrg, frc)
    case(3)
        call vdw_cosine(rij_mag, vdw_params(:,typ), enrg, frc)
    case(4)
        call vdw_lj_coul_debye(rij_mag, qi*qj, vdw_params(:,typ), enrg, frc)
    case(5)
        call vdw_lj_coul(rij_mag, qi*qj, vdw_params(:,typ), enrg, frc)
    case(6)
        call vdw_dpd(rij_mag, vdw_params(:,typ), enrg, frc)
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine vdw_lj_set(params, eps, sigma, rcut)
    !! Setter for 12-6 LJ (truncated & force-shifted) interaction.
    !!
    !! The potential `U` is given by:
    !!```
    !!   V = 4*eps*[(r/sigma)^12 - (r/sigma)^6]  
    !!   U = V - V(rcut) - (r - rcut)*dV/dr, r < rcut,
    !!       0, r >= rcut  
    !!```
    !! where `dV/dr` is evaluated at `r = rcut`.
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `eps`
    !! * params(2) = `sigma`
    !! * params(3) = `rcut`
    !!
    !! Internally stored parameters:
    !!
    !! * params(4) = `V(rcut)`
    !! * params(5) = `dV/dr(rcut)`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: eps
    real(rp), intent(in), optional :: sigma
    real(rp), intent(in), optional :: rcut
    real(rp) :: eps_, sigma_, rcut_
    real(rp) :: pot_rcut, pot_deriv_rcut
    real(rp) :: sir, sir2, sir6, sir12

    if (present(eps)) params(1) = eps
    if (present(sigma)) params(2) = sigma
    if (present(rcut)) params(3) = rcut

    eps_ = params(1); sigma_ = params(2); rcut_ = params(3)

    sir = sigma_/rcut_
    sir2 = sir*sir; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
    pot_rcut = 4*eps_*(sir12 - sir6)
    pot_deriv_rcut = -24*eps_*(2*sir12 - sir6)/rcut_

    params(4) = pot_rcut
    params(5) = pot_deriv_rcut

    end subroutine

!******************************************************************************

pure subroutine vdw_lj(r, params, enrg, frc)
    !! Evaluates the potential and its derivative for LJ interaction. See
    !! [[vdw_lj_set]].

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    real(rp) :: eps, sigma, rcut
    real(rp) :: pot_rcut, pot_deriv_rcut
    real(rp) :: sir, sir2, sir6, sir12

    enrg = 0.0_rp; frc = 0.0_rp

    eps = params(1); sigma = params(2); rcut = params(3)
    pot_rcut = params(4)
    pot_deriv_rcut = params(5)

    if ( r < rcut ) then
        sir = sigma/r
        sir2 = sir*sir; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
        enrg = 4*eps*(sir12 - sir6) - pot_rcut - (r - rcut)*pot_deriv_rcut
        frc = -24*eps*(2*sir12 - sir6)/r - pot_deriv_rcut
    end if

    end subroutine

!******************************************************************************

subroutine vdw_gaussian_set(params, A, B, rcut)
    !! Setter for gaussian interaction. The potential is truncated and
    !! force-shifted.
    !!
    !! The potential `U` is given by:
    !!```
    !!   V = A*exp(-B*r^2)  
    !!   U = V - V(rcut) - (r - rcut)*dV/dr, r < rcut  
    !!       0, r >= rcut,
    !!```
    !! where `dV/dr` is evaluated at `r = rcut`.
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `A`
    !! * params(2) = `B`
    !! * params(3) = `rcut`
    !!
    !! Internally stored parameters:
    !!
    !! * params(4) = `V(rcut)`
    !! * params(5) = `dV/dr(rcut)`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: A
    real(rp), intent(in), optional :: B
    real(rp), intent(in), optional :: rcut
    real(rp) :: pot_rcut, pot_deriv_rcut

    if (present(A)) params(1) = A
    if (present(B)) params(2) = B
    if (present(rcut)) params(3) = rcut

    pot_rcut = params(1)*exp(-params(2)*params(3)**2)
    pot_deriv_rcut = -2.0_rp*params(1)*params(2)*params(3) &
                     * exp(-params(2)*params(3)**2)

    params(4) = pot_rcut
    params(5) = pot_deriv_rcut

    end subroutine

!******************************************************************************

pure subroutine vdw_gaussian(r, params, enrg, frc)
    !! Calculates energy & its derivative for gaussian interaction. See
    !! [[vdw_gaussian_set]].

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    real(rp) :: A, B, rcut, pot_rcut, pot_deriv_rcut
    real(rp) :: exrs

    enrg = 0.0_rp; frc = 0.0_rp
    A = params(1); B = params(2); rcut = params(3)
    pot_rcut = params(4); pot_deriv_rcut = params(5)

    exrs = exp(-B*r*r)

    if ( r < rcut ) then
        enrg = A*exrs - pot_rcut - (r - rcut)*pot_deriv_rcut
        frc = -2*A*B*r*exrs - pot_deriv_rcut
    end if

    end subroutine

!******************************************************************************

subroutine vdw_cosine_set(params, A, rcut)
    !! Setter for cosine interaction.
    !!
    !! The potential `U` is given by:
    !!```
    !!   U = A*[1 + cos(pi*r/rcut)], r < rcut
    !!       0, r >= rcut
    !!```
    !! The potential as well as its derivative is zero at `r = rcut`.
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `A`
    !! * params(2) = `rcut`
    !!
    !! Internally stored parameters:
    !!
    !! * None

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: A
    real(rp), intent(in), optional :: rcut

    if (present(A)) params(1) = A
    if (present(rcut)) params(2) = rcut

    end subroutine

!******************************************************************************

pure subroutine vdw_cosine(r, params, enrg, frc)
    !! Evaluates the potential and its derivative for  cosine interaction.
    !! See [[vdw_cosine_set]].

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    real(rp) :: A, rcut, pf, pr

    enrg = 0.0_rp; frc = 0.0_rp
    A = params(1); rcut = params(2)

    pf = math_pi/rcut; pr = pf*r

    if ( r < rcut ) then
        enrg = A*( 1.0_rp + cos(pr) )
        frc = -A*pf*sin(pr)
    end if

    end subroutine

!******************************************************************************

subroutine vdw_lj_coul_debye_set(params, eps, sigma, rcut, rcut_coul, C, kappa)
    !! Setter for 12-6 LJ with screened Coulombic interaction.
    !!
    !! The potential `U` is given by:
    !!```
    !! V = 4*eps*[(r/sigma)^12 - (r/sigma)^6]
    !! W = C*qi*qj*exp(-kappa*r)/r
    !! if rcut_coul > 0:
    !!     U = V - V(rcut) + W - W(rcut_coul), r < rcut
    !!         W - W(rcut_coul), rcut <= r < rcut_coul
    !!         0, r >= rcut_coul
    !! if rcut_coul <= 0:
    !!     U = V - V(rcut) + W, r < rcut
    !!         W, r >= rcut
    !!```
    !!
    !!- The LJ potential `V` is cut & shifted at `r = rcut`.  
    !!- If `rcut_coul > 0`, the screened Coulombic potential `W` is cut & shifted
    !! at `r = rcut_coul`.  
    !!- If `rcut_coul <= 0`, no cutoff is applied on `W`.  
    !!- `rcut_coul` must be >= `rcut`.  
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `eps`
    !! * params(2) = `sigma`
    !! * params(3) = `rcut`
    !! * params(4) = `rcut_coul`
    !! * params(5) = `C`
    !! * params(6) = `kappa`
    !!
    !! Internally stored parameters:
    !!
    !! * params(7) = `V(rcut)`
    !! * params(8) = `C*exp(-kappa*rcut_coul)/rcut_coul`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: eps
    real(rp), intent(in), optional :: sigma
    real(rp), intent(in), optional :: rcut
    real(rp), intent(in), optional :: rcut_coul
    real(rp), intent(in), optional :: C
    real(rp), intent(in), optional :: kappa
    real(rp) :: eps_, sigma_, rcut_, rcut_coul_, C_, kappa_
    real(rp) :: pot_rcut, pot_rcut_coul
    real(rp) :: sir, sir2, sir6, sir12

    if (present(eps)) params(1) = eps
    if (present(sigma)) params(2) = sigma
    if (present(rcut)) params(3) = rcut
    if (present(rcut_coul)) params(4) = rcut_coul
    if (present(C)) params(5) = C
    if (present(kappa)) params(6) = kappa

    eps_ = params(1); sigma_ = params(2); rcut_ = params(3)
    rcut_coul_ = params(4); C_ = params(5); kappa_ = params(6)

    sir = sigma_/rcut_
    sir2 = sir*sir; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
    pot_rcut = 4*eps_*(sir12 - sir6)

    if (rcut_coul_ > 0.0_rp) then
        pot_rcut_coul = C_*exp(-kappa_*rcut_coul_)/rcut_coul_
    else
        pot_rcut_coul = 0.0_rp
    end if

    params(7) = pot_rcut
    params(8) = pot_rcut_coul

    end subroutine

!******************************************************************************

pure subroutine vdw_lj_coul_debye(r, qiqj, params, enrg, frc)
    !!Evaluates the potential and its derivative for screened Coulombic interaction
    !!combined with 12-6 LJ (cut & shifted). See [[vdw_lj_coul_debye_set]].

    real(rp), intent(in) :: r
    real(rp), intent(in) :: qiqj
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    real(rp) :: eps, sigma, rcut, rcut_coul, C, kappa
    real(rp) :: pot_rcut, pot_rcut_coul
    real(rp) :: sir, sir2, sir6, sir12, ekr

    enrg = 0.0_rp; frc = 0.0_rp

    eps = params(1); sigma = params(2); rcut = params(3)
    rcut_coul = params(4); C = params(5); kappa = params(6)
    pot_rcut = params(7); pot_rcut_coul = params(8)

    ekr = C*qiqj*exp(-kappa*r)/r
    if (rcut_coul > 0.0_rp) then
        !Coulombic cutoff at rcut_coul
        if ( r < rcut_coul ) then
            enrg = ekr - qiqj*pot_rcut_coul
            frc = -ekr*(1+kappa*r)/r
        end if
    else
        !No Coulombic cutoff
        enrg = ekr
        frc = -ekr*(1+kappa*r)/r
    end if

    if ( r < rcut ) then
        sir = sigma/r
        sir2 = sir*sir; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
        enrg = enrg + 4*eps*(sir12 - sir6) - pot_rcut
        frc = frc - 24*eps*(2*sir12 - sir6)/r
    end if

    end subroutine

!******************************************************************************

subroutine vdw_lj_coul_set(params, eps, sigma, rcut, rcut_coul, C)
    !! Setter for 12-6 LJ with Coulombic interaction.
    !!
    !! The potential `U` is given by:
    !!```
    !! V = 4*eps*[(r/sigma)^12 - (r/sigma)^6]
    !! W = C*qi*qj/r
    !! if rcut_coul > 0:
    !!     U = V - V(rcut) + W - W(rcut_coul), r < rcut
    !!         W - W(rcut_coul), rcut <= r < rcut_coul
    !!         0, r >= rcut_coul
    !! if rcut_coul <= 0:
    !!     U = V - V(rcut) + W, r < rcut
    !!         W, r >= rcut
    !!```
    !! - The LJ potential `V` is cut & shifted at `r = rcut`.
    !! - If `rcut_coul > 0`, the Coulombic potential `W` is cut & shifted at `r = rcut_coul`.
    !! - If `rcut_coul <= 0`, no cutoff is applied on `W`.
    !! - `rcut_coul` must be >= `rcut`.
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `eps`
    !! * params(2) = `sigma`
    !! * params(3) = `rcut`
    !! * params(4) = `rcut_coul`
    !! * params(5) = `C`
    !!
    !! Internally stored parameters:
    !!
    !! * params(6) = `V(rcut)`
    !! * params(7) = `C/rcut_coul`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: eps
    real(rp), intent(in), optional :: sigma
    real(rp), intent(in), optional :: rcut
    real(rp), intent(in), optional :: rcut_coul
    real(rp), intent(in), optional :: C
    real(rp) :: eps_, sigma_, rcut_, rcut_coul_, C_
    real(rp) :: pot_rcut, pot_rcut_coul
    real(rp) :: sir, sir2, sir6, sir12

    if (present(eps)) params(1) = eps
    if (present(sigma)) params(2) = sigma
    if (present(rcut)) params(3) = rcut
    if (present(rcut_coul)) params(4) = rcut_coul
    if (present(C)) params(5) = C

    eps_ = params(1); sigma_ = params(2); rcut_ = params(3)
    rcut_coul_ = params(4); C_ = params(5)

    sir = sigma_/rcut_
    sir2 = sir*sir; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
    pot_rcut = 4*eps_*(sir12 - sir6)

    if (rcut_coul_ > 0.0_rp) then
        pot_rcut_coul = C_/rcut_coul_
    else
        pot_rcut_coul = 0.0_rp
    end if

    params(6) = pot_rcut
    params(7) = pot_rcut_coul

    end subroutine

!******************************************************************************

pure subroutine vdw_lj_coul(r, qiqj, params, enrg, frc)
    !! Evaluates the potential and its derivative for Coulombic interaction
    !! combined with 12-6 LJ (cut & shifted). See [[vdw_lj_coul_set]].

    real(rp), intent(in) :: r
    real(rp), intent(in) :: qiqj
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    real(rp) :: eps, sigma, rcut, rcut_coul, C
    real(rp) :: pot_rcut, pot_rcut_coul
    real(rp) :: sir, sir2, sir6, sir12, ekr

    enrg = 0.0_rp; frc = 0.0_rp

    eps = params(1); sigma = params(2); rcut = params(3)
    rcut_coul = params(4); C = params(5)
    pot_rcut = params(6); pot_rcut_coul = params(7)

    ekr = C*qiqj/r
    if (rcut_coul > 0.0_rp) then
        !Coulombic cutoff at rcut_coul
        if ( r < rcut_coul ) then
            enrg = ekr - qiqj*pot_rcut_coul
            frc = -ekr/r
        end if
    else
        !No Coulombic cutoff
        enrg = ekr
        frc = -ekr/r
    end if

    if ( r < rcut ) then
        sir = sigma/r
        sir2 = sir*sir; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
        enrg = enrg + 4*eps*(sir12 - sir6) - pot_rcut
        frc = frc - 24*eps*(2*sir12 - sir6)/r
    end if

    end subroutine

!******************************************************************************

subroutine vdw_dpd_set(params, A, rcut)
    !! Setter for standard DPD interaction.
    !!
    !! The potential `U` is given by:
    !!```
    !!   U = (A/2)*rcut*(1 - (r/rcut))^2, r < rcut
    !        0, r >= rcut
    !!```
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `A`
    !! * params(2) = `rcut`
    !!
    !! Internally stored parameters:
    !!
    !! * None

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: A
    real(rp), intent(in), optional :: rcut

    if (present(A)) params(1) = A
    if (present(rcut)) params(2) = rcut

    end subroutine

!******************************************************************************

pure subroutine vdw_dpd(r, params, enrg, frc)
    !! Evaluates the potential and its derivative for standard DPD interaction.
    !! See [[vdw_dpd_set]].

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    real(rp) :: A, rcut, ir

    enrg = 0.0_rp; frc = 0.0_rp
    A = params(1); rcut = params(2)

    ir = 1.0_rp - (r/rcut)

    if ( r < rcut ) then
        enrg = 0.5_rp*A*rcut*ir*ir
        frc = -A*ir
    end if

    end subroutine

!******************************************************************************

subroutine vdw_expnrx_set(params,eps,alpha,Rm,exponent,sign,rcut)
    !! Setter for expn/rx interaction.
    !!
    !! The potential `U` is given by:
    !!```
    !!   V = eps/(alpha-6)*(sign*6*exp(alpha*(1-rij/Rm))-alpha*(Rm/rij)**exponent)
    !!   U = V - V(rcut) - (r-rcut)*dV/dr, r < rcut
    !        0, r >= rcut
    !!```
    !! where `dV/dr` is evaluated at `r = rcut`.
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `eps`
    !! * params(2) = `alpha`
    !! * params(3) = `Rm`
    !! * params(4) = `exponent`
    !! * params(5) = `sign` (optional:1,-1)
    !! * params(6) = `rcut`
    !!
    !! Internally stored parameters:
    !! * params(7) = `V(rcut)`
    !! * params(8) = `dV/dr(rcut)`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: eps
    real(rp), intent(in), optional :: alpha
    real(rp), intent(in), optional :: Rm
    integer, intent(in), optional :: exponent
    integer, intent(in), optional :: sign
    real(rp), intent(in), optional :: rcut
    real(rp) :: pot_rcut, pot_deriv_rcut
    real(rp) :: RmOverrcut

    if (present(eps)) params(1) = eps
    if (present(alpha)) params(2) = alpha
    if (present(Rm)) params(3) = Rm
    if (present(exponent)) params(4) = exponent
    if (present(sign)) params(5) = sign
    if (present(rcut)) params(6) = rcut

    RmOverrcut = params(3)/params(6);
    pot_rcut = params(1)/(params(2)-6)*(params(5)*6*exp(params(2) &
                *(1-1/RmOverrcut))-params(2)*RmOverrcut**params(4))
    pot_deriv_rcut = params(1)/(params(2)-6)*(-params(2)/params(3) &
                *params(5)*6*exp(params(2)*(1-1/RmOverrcut)) &
                +params(4)*params(2)*RmOverrcut**(params(4)+1))
    
    params(7) = pot_rcut
    params(8) = pot_deriv_rcut

    end subroutine

!******************************************************************************

pure subroutine vdw_expnrx(r, params, enrg, frc)
    !! Evaluates the potential and its derivative for expnrx interaction.
    !! See [[vdw_expnrx_set]].
    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    real(rp) :: eps, alpha, Rm, rcut
    real(rp) :: pot_rcut, pot_deriv_rcut
    integer :: sign, exponent
    real(rp) :: RmOverrcut
    
    enrg = 0.0_rp; frc = 0.0_rp
    
    eps = params(1); alpha = params(2); Rm = params(3)
    exponent = params(4); sign = params(5)
    rcut = params(6); pot_rcut = params(7); pot_deriv_rcut = params(8)
    
    if (r < rcut) then
        enrg = eps/(alpha-6)*(sign*6*exp(alpha*(1-r/Rm)) &
                -alpha*(Rm/r)**exponent) - pot_rcut - (r-rcut)*pot_deriv_rcut
        frc = eps/(alpha-6)*(-alpha/Rm*sign*6*exp(alpha*(1-r/Rm)) &
                +alpha*exponent*(Rm/r)**(exponent+1))
    end if

    end subroutine

!******************************************************************************

end module m_ia_vdw
