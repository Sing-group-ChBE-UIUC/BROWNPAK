module m_ia_tether
!! Tether potentials
!!
!! *Style 0: None
!! *Style 1: Rigid connector (not implemented)
!! *Style 2: Harmonic spring. See [[teth_harm_set]].

use m_precision
use m_constants_math
use m_globals

implicit none

private

public :: ia_tether_setup, ia_get_tether_force

contains

!******************************************************************************

subroutine ia_tether_setup()
    !! Sets up parameters for tether potentials 

    integer :: i
    integer :: sty

    !Set tether interactions
    do i = 1, num_tether_types
        sty = tether_styles(i)
        select case(sty)
        case(1)
            call teth_rigid_set(tether_params(:,i))
        case(2)
            call teth_harm_set(tether_params(:,i))
        case default
            continue
        end select
    end do

    end subroutine

!******************************************************************************

subroutine ia_get_tether_force(qmag, teth_typ, enrg, frc, ierr)
    !! Calculates the energy and its derivative due to a tether.

    real(rp), intent(in) :: qmag
        !! Distance between the tethered atom & the tether point
    integer, intent(in)   :: teth_typ
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    integer, intent(out)  :: ierr
    integer :: teth_styl

    ierr = 0
    enrg = 0.0_rp; frc = 0.0_rp

    teth_styl = tether_styles(teth_typ)

    select case (teth_styl)
    case (1)
        call teth_rigid(qmag, tether_params(:,teth_typ), enrg, frc, ierr) 
    case (2)
        call teth_harm(qmag, tether_params(:,teth_typ), enrg, frc) 
    case default
        continue
    end select

    end subroutine

!********************************************************************************

subroutine teth_rigid_set(params, r0, eps)
    !! Setter for rigid tether interaction.
    !!
    !! User-set parameters:
    !!
    !! * params(1) = `r0` (fixed separation distance)
    !! * params(2) = `eps` (allowed tolerance)

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: r0
    real(rp), intent(in), optional :: eps

    if (present(r0)) params(1) = r0
    if (present(eps)) params(2) = eps

    end subroutine

!********************************************************************************

subroutine teth_rigid(r, params, enrg, frc, ierr)
    !! Not implemented, needs constraint formalism. See [[teth_rigid_set]].
    !!
    !! Calculates energy for rigid tether interaction.

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc
    integer, intent(out) :: ierr

    ierr = 0; enrg = 0.0_rp; frc = 0.0_rp

    if ( abs(r-params(1)) < params(2) ) then
        enrg = 0.0_rp; frc = 0.0_rp
    else
        ierr = 1
        return
    end if

    end subroutine

!********************************************************************************

!params(1) = k, params(2) = r0 (equilibrium distance)
subroutine teth_harm_set(params, k, r0)
    !! Setter for harmonic tether interaction.
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

subroutine teth_harm(r, params, enrg, frc)
    !! Calculates energy and its derivative for harmonic tether interaction. See
    !! [[teth_harm_set]].

    real(rp), intent(in) :: r
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), intent(out) :: frc

    enrg = 0.5_rp*params(1)*(r-params(2))*(r-params(2))
    frc = params(1)*(r-params(2))

    end subroutine

!******************************************************************************

end module m_ia_tether
