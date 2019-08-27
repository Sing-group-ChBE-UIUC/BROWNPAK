module m_ia_dihedral
!! Dihedral potentials (none implemented)

use m_precision
use m_constants_math
use m_globals

implicit none

private

public :: ia_dihedral_setup, ia_get_dihedral_force

contains

!******************************************************************************

subroutine ia_dihedral_setup()
    !! Sets up parameters for dihedral potentials 

    integer :: i
    integer :: sty

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

subroutine ia_get_dihedral_force(q1, q2, q3, typ, enrg, fi, fj, fk, fl)
    !! Calculates the force & energy due to a dihedral.

    real(rp), dimension(3), intent(in)  :: q1
    real(rp), dimension(3), intent(in)  :: q2
    real(rp), dimension(3), intent(in)  :: q3
    real(rp), dimension(3), intent(out) :: fi
    real(rp), dimension(3), intent(out) :: fj
    real(rp), dimension(3), intent(out) :: fk
    real(rp), dimension(3), intent(out) :: fl
    integer, intent(in) :: typ
    real(rp), intent(out) :: enrg
    integer :: dhd_styl

    enrg = 0.0_rp
    fi = 0.0_rp; fj = 0.0_rp; fk = 0.0_rp; fl = 0.0_rp

    dhd_styl = dihedral_styles(typ)

    select case(dhd_styl)
    case default
        continue
    end select

    end subroutine

!******************************************************************************

end module m_ia_dihedral
