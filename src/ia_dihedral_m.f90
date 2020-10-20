module ia_dihedral_m
!! Dihedral potentials (none implemented)

use constants_m
use atmcfg_m

implicit none

private

public :: ia_dihedral_setup, ia_get_dihedral_force

contains

!******************************************************************************

subroutine ia_dihedral_setup(num_dihedral_types, dihedral_styles, dihedral_params)
    !! Sets up parameters for dihedral potentials 

    integer, intent(in) :: num_dihedral_types
        !! Number of dihedral types
    integer, dimension(:), intent(in) :: dihedral_styles
        !! Styles for each type
    real(rp), dimension(:,:), intent(in out) :: dihedral_params
        !! Parameters for each type, depending on style
    integer :: i, sty

    !Set dihedral interactions
    do i = 1, num_dihedral_types
        sty = dihedral_styles(i)
        select case(sty)
        case default
            continue
        end select
    end do

    end subroutine

!******************************************************************************

subroutine ia_get_dihedral_force(q1, q2, q3, sty, params, enrg, fi, fj, fk, fl)
    !! Calculates the force & energy due to a dihedral.

    real(rp), dimension(3), intent(in)  :: q1
    real(rp), dimension(3), intent(in)  :: q2
    real(rp), dimension(3), intent(in)  :: q3
    integer,                intent(in)  :: sty
    real(rp), dimension(:), intent(in)  :: params
    real(rp), dimension(3), intent(out) :: fi
    real(rp), dimension(3), intent(out) :: fj
    real(rp), dimension(3), intent(out) :: fk
    real(rp), dimension(3), intent(out) :: fl
    real(rp), intent(out) :: enrg

    enrg = 0.0_rp
    fi = 0.0_rp; fj = 0.0_rp; fk = 0.0_rp; fl = 0.0_rp

    select case(sty)
    case default
        continue
    end select

    end subroutine

!******************************************************************************

end module ia_dihedral_m
