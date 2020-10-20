module ia_angle_m
!! Angle potentials
!!
!! * Style 0: None (Only topology)
!! * Style 1: Cosine. See [[ang_cos_set]].

use constants_m
use atmcfg_m

implicit none

private

public :: ia_angle_setup, ia_get_angle_force

contains

!******************************************************************************

subroutine ia_angle_setup(num_angle_types, angle_styles, angle_params)
    !! Sets up parameters for angle potentials 

    integer, intent(in) :: num_angle_types
        !! Number of angle types
    integer, dimension(:), intent(in) :: angle_styles
        !! Styles for each type
    real(rp), dimension(:,:), intent(in out) :: angle_params
        !! Parameters for each type, depending on style
    integer :: i, sty

    !Set angular interactions
    do i = 1, num_angle_types
        sty = angle_styles(i)
        select case(sty)
        case(1)
            call ang_cos_set(angle_params(:,i))
        case default
            continue
        end select
    end do

    end subroutine

!******************************************************************************

subroutine ia_get_angle_force(q1, q2, sty, params, enrg, fim1, fi, fip1)
    !! Calculates the energy & force due to an angle.

    real(rp), dimension(3), intent(in) :: q1
    real(rp), dimension(3), intent(in) :: q2
    integer, intent(in) :: sty
    real(rp), dimension(:), intent(in) :: params
    real(rp), intent(out) :: enrg
    real(rp), dimension(3), intent(out) :: fim1
    real(rp), dimension(3), intent(out) :: fi
    real(rp), dimension(3), intent(out) :: fip1
    real(rp), dimension(3) :: q1hat, q2hat
    real(rp) :: kang
    real(rp) :: q1mag, q2mag, ctheta

    q1mag = norm2(q1); q2mag = norm2(q2)
    q1hat = q1/q1mag; q2hat = q2/q2mag
    ctheta = dot_product(q1hat, q2hat)
    !Floating point correction
    if (ctheta > 1.0_rp) ctheta = 1.0_rp 
    if (ctheta < -1.0_rp) ctheta = -1.0_rp

    select case(sty)
    case(1)
        kang = params(1)
        enrg = kang*(1.0_rp - ctheta)
        fim1 = kang*(-q2hat + ctheta*q1hat)/q1mag
        fip1 = kang*(q1hat - ctheta*q2hat)/q2mag
        fi = -(fim1 + fip1)
    case default
        continue
    end select

    end subroutine

!********************************************************************************

subroutine ang_cos_set(params, k)
    !! Setter for angular cosine interaction.
    !!
    !!```
    !!   U(theta) = k*(1 - cos theta),
    !!   where theta is the complementary angle between bonds i & (i+1).
    !!```
    !! User-set parameters:
    !!
    !! params(1) = `k`

    real(rp), dimension(:), intent(in out) :: params
    real(rp), intent(in), optional :: k

    if (present(k)) params(1) = k

    end subroutine

!******************************************************************************

end module ia_angle_m
