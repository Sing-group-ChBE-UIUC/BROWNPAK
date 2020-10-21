!********************************************************************************!
!                                                                                !
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2020 Sarit Dutta                                                 !
!                                                                                !
! Permission is hereby granted, free of charge, to any person obtaining a copy   !
! of this software and associated documentation files (the "Software"), to deal  !
! in the Software without restriction, including without limitation the rights   !
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      !
! copies of the Software, and to permit persons to whom the Software is          !
! furnished to do so, subject to the following conditions:                       !
!                                                                                !
! The above copyright notice and this permission notice shall be included in all !
! copies or substantial portions of the Software.                                !
!                                                                                !
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     !
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       !
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    !
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         !
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  !
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  !
! SOFTWARE.                                                                      !
!                                                                                !
!********************************************************************************!

module ia_external_m
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

use constants_m
use atmcfg_m

implicit none

private

public :: ia_external_setup, ia_add_external_forces

contains

!******************************************************************************

subroutine ia_external_setup(num_externals, external_styles, external_params)
    !! Sets up parameters for external potentials. Usually there is nothing to
    !! set for externals, but this acts as a placeholder for special cases.

    integer, intent(in) :: num_externals
        !! Number of external fields
    integer, dimension(:), intent(in) :: external_styles
        !! Styles for each field
    real(rp), dimension(:,:), intent(in out) :: external_params
        !! Parameters for each field, depending on style
    integer :: i, sty

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

subroutine ia_add_external_forces(num_externals, external_styles, &
        external_params, coordinates, energy_external, forces, stress, ierr)
    !! Calculates the force and energy due to an external field and adds to
    !! `energy_external`, 'forces`, & `stress`.

    integer, intent(in) :: num_externals
        !! Number of external fields
    integer, dimension(:), intent(in) :: external_styles
        !! Styles for each field
    real(rp), dimension(:,:), intent(in out) :: external_params
        !! Parameters for each field, depending on style
    real(rp), dimension(:,:), intent(in) :: coordinates
    real(rp), intent(out) :: energy_external
    real(rp), dimension(:,:), intent(in out) :: forces
    real(rp), dimension(3,3), intent(in out) :: stress
    integer, intent(out)   :: ierr
    real(rp) :: enrg
    real(rp) :: frcx, v, sn
    integer  :: iext
    integer  :: sty, iatm, m

    ierr = 0; energy_external = 0.0_rp

    do iext = 1, num_externals
        sty = external_styles(iext)

        select case(sty)
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

end module ia_external_m
