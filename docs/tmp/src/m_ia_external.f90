module m_ia_external
!! External potentials
!!
!! This module is meant to be a placeholder to any external fields that the user
!! wants to add. Accordingly replace/add to the existing routines. The
!! subroutines [[ia_external_setup]] and [[ia_add_external_forces]] must remain
!! for interfacing to the force calculation driver routine [[ia_calc_forces]].
!!
!! * Style 0: None
!! * Style 1: Pulling force along +ve x-axis
!! * Style 2: Hard planar wall

use m_precision
use m_constants_math
use m_globals

implicit none

private

public :: ia_external_setup, ia_add_external_forces

contains

!******************************************************************************

subroutine ia_external_setup()
    !! Sets up parameters for external potentials. Usually there is nothing to
    !! set for externals, but this acts as a placeholder for special cases.

    integer :: i
    integer :: sty

    !Set external interactions
    do i = 1, num_externals
        sty = external_styles(i)
        select case(sty)
        case default
            continue
        end select
    end do

    end subroutine

!******************************************************************************

subroutine ia_add_external_forces(ierr)
    !! Calculates the force and energy due to an external field and adds to
    !! `energy_external`, 'forces`, & `stress` in module `m_globals`.

    integer, intent(out)   :: ierr
    real(rp) :: enrg
    real(rp) :: frcx, v, sn
    integer  :: iext
    integer  :: styl, iatm, m

    ierr = 0
    do iext = 1, num_externals
        styl = external_styles(iext)

        select case(styl)
        case(1)
            ! Pulling force along +ve x-axis.
            iatm = int(external_params(1,iext))
            frcx = external_params(2,iext)
            enrg = -frcx*(coordinates(1,iatm) - coordinates(1,1))
            energy_external = energy_external + enrg
            forces(1,iatm) = forces(1,iatm) + frcx
        case(2)
            ! Rigid walls. Need to modify this (or another case) for repulsive walls.
            enrg = 0.0_rp
            m = int(external_params(1,iext)) 
            v = external_params(2,iext)
            sn = external_params(3,iext)
            if (sn > 0.0_rp) then
                if ( any(coordinates(m,:) < v) ) ierr = 1
            else
                if ( any(coordinates(m,:) > v) ) ierr = 1
            end if
        case default
            continue
        end select
    end do

    !Need to update stress

    end subroutine

!******************************************************************************

end module m_ia_external
