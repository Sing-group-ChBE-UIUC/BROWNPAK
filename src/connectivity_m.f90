!********************************************************************************!
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2020 Sarit Dutta <saridut@gmail.com>                             !
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
!********************************************************************************!

module connectivity_m
    !! Routines for building atom->bond, atom->angle, etc. tables and excluded
    !! atoms table.

use constants_m
use vector_m
use table_m

implicit none

private

public :: atbo_build, atat_build, atan_build, atdh_build, exat_build

contains

!********************************************************************************

subroutine atbo_build(num_atoms, num_bonds, bonds, atbo_tab)

    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_bonds
    integer, dimension(:,:), intent(in) :: bonds
    type(itable_t), intent(out) :: atbo_tab
        !! Atoms -> bonds table
    type(ivector_t), dimension(:), allocatable :: buf_map
    integer :: iatm, jatm, ibnd
    integer :: j

    !Build table as list of lists
    allocate(buf_map(num_atoms))
    do iatm = 1, num_atoms
        call ivector_init(buf_map(iatm))
    end do

    do ibnd = 1, num_bonds
        iatm = bonds(2,ibnd)
        jatm = bonds(3,ibnd)
        call buf_map(iatm)%append(ibnd)
        call buf_map(jatm)%append(ibnd)
    end do

    !Sort in ascending order
    do iatm = 1, num_atoms
        if (buf_map(iatm)%len > 1) call buf_map(iatm)%sort()
    end do

    !Copy over to table
    call itbl_init(atbo_tab, num_atoms)
    do iatm = 1, num_atoms
        do j = 1, buf_map(iatm)%len
            jatm = buf_map(iatm)%get_val(j)
            call atbo_tab%append(iatm, jatm)
        end do
    end do

    !Delete list of lists table
    do iatm = 1, num_atoms
        call buf_map(iatm)%delete()
    end do
    deallocate(buf_map)

    !Release additional memory
    call atbo_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine atan_build(num_atoms, num_angles, angles, atan_tab)

    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_angles
    integer, dimension(:,:), intent(in) :: angles
    type(itable_t), intent(out) :: atan_tab
        !! Atoms -> angles table
    type(ivector_t), dimension(:), allocatable :: buf_map
    integer :: iatm, jatm, katm, iang
    integer :: j

    !Build table as list of lists
    allocate(buf_map(num_atoms))
    do iatm = 1, num_atoms
        call ivector_init(buf_map(iatm))
    end do

    do iang = 1, num_angles
        iatm = angles(2,iang)
        jatm = angles(3,iang)
        katm = angles(4,iang)
        call buf_map(iatm)%append(iang)
        call buf_map(jatm)%append(iang)
        call buf_map(katm)%append(iang)
    end do

    !Sort in ascending order
    do iatm = 1, num_atoms
        if (buf_map(iatm)%len > 1) call buf_map(iatm)%sort()
    end do

    !Copy over to table
    call itbl_init(atan_tab, num_atoms)
    do iatm = 1, num_atoms
        do j = 1, buf_map(iatm)%len
            jatm = buf_map(iatm)%get_val(j)
            call atan_tab%append(iatm, jatm)
        end do
    end do

    !Delete list of lists table
    do iatm = 1, num_atoms
        call buf_map(iatm)%delete()
    end do
    deallocate(buf_map)

    !Release additional memory
    call atan_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine atdh_build(num_atoms, num_dihedrals, dihedrals, atdh_tab)

    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_dihedrals
    integer, dimension(:,:), intent(in) :: dihedrals
    type(itable_t), intent(out) :: atdh_tab
        !! Atoms -> dihedrals table
    type(ivector_t), dimension(:), allocatable :: buf_map
    integer :: iatm, jatm, katm, latm, idhd
    integer :: j

    !Build table as list of lists
    allocate(buf_map(num_atoms))
    do iatm = 1, num_atoms
        call ivector_init(buf_map(iatm))
    end do

    do idhd = 1, num_dihedrals
        iatm = dihedrals(2,idhd)
        jatm = dihedrals(3,idhd)
        katm = dihedrals(4,idhd)
        latm = dihedrals(5,idhd)
        call buf_map(iatm)%append(idhd)
        call buf_map(jatm)%append(idhd)
        call buf_map(katm)%append(idhd)
        call buf_map(latm)%append(idhd)
    end do

    !Sort in ascending order
    do iatm = 1, num_atoms
        if (buf_map(iatm)%len > 1) call buf_map(iatm)%sort()
    end do

    !Copy over to table
    call itbl_init(atdh_tab, num_atoms)
    do iatm = 1, num_atoms
        do j = 1, buf_map(iatm)%len
            jatm = buf_map(iatm)%get_val(j)
            call atdh_tab%append(iatm, jatm)
        end do
    end do

    !Delete list of lists table
    do iatm = 1, num_atoms
        call buf_map(iatm)%delete()
    end do
    deallocate(buf_map)

    !Release additional memory
    call atdh_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine atat_build(num_atoms, bonds, atbo_tab, atat_tab)

    integer, intent(in) :: num_atoms
    integer, dimension(:,:), intent(in) :: bonds
    type(itable_t), intent(in) :: atbo_tab
        !! Atoms -> bonds table
    type(itable_t), intent(out) :: atat_tab
        !! Atoms -> bonded atoms table (1-ring)
    integer, dimension(:), pointer :: inc_bonds => null()
    type(ivector_t) :: excl_atms
    integer :: iatm, jatm, jbnd, jbnd_atm1, jbnd_atm2
    integer :: j

    !Initialize table
    call itbl_init(atat_tab, num_atoms)
    !List of excluded atoms
    call ivector_init(excl_atms)

    do iatm = 1, num_atoms
        call excl_atms%clear()
        !Get incident bonds
        call atbo_tab%get_row(iatm, inc_bonds)
        !Add atoms from each of the incident bonds, excluding iatom
        do j = 1, size(inc_bonds)
            jbnd = inc_bonds(j)
            jbnd_atm1 = bonds(2,jbnd)
            jbnd_atm2 = bonds(3,jbnd)
            if (jbnd_atm1 /= iatm) call excl_atms%append(jbnd_atm1)
            if (jbnd_atm2 /= iatm) call excl_atms%append(jbnd_atm2)
        end do
        !Sort and remove duplicates. These are the atoms in the 1-ring.
        call excl_atms%unique()

        !Add to atat_tab. Note that iatm does not appear in the atat table
        !for iatm.
        do j = 1, excl_atms%len
            jatm = excl_atms%get_val(j)
            call atat_tab%append(iatm, jatm)
        end do
    end do

    !Release additional memory
    call atat_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

subroutine exat_build(num_atoms, excluded_atoms, atat_tab, exat_tab)
    !! Before a call to this first build atbo_tab.

    integer, intent(in) :: num_atoms
    integer, intent(in) :: excluded_atoms
    type(itable_t), intent(in) :: atat_tab
        !! Atoms -> bonded atoms table
    type(itable_t), intent(out) :: exat_tab
        !! Atoms -> excluded (from vdw calculation) atoms table
    integer, dimension(:), pointer :: nbr_atms => null()
    integer, dimension(:), pointer :: nbr2_atms => null()
    integer, dimension(:), pointer :: nbr3_atms => null()
    type(ivector_t) :: excl_atms
    integer :: iatm, jatm, j2atm, j3atm
    integer :: j, j2, j3

    !Initialize table
    call itbl_init(exat_tab, num_atoms)
    !List of excluded atoms
    call ivector_init(excl_atms)

    do iatm = 1, num_atoms
        call excl_atms%clear()
        if (excluded_atoms > 0) then
            !Get atoms in 1-ring neighborhood
            call atat_tab%get_row(iatm, nbr_atms)
            !Add atoms to excluded atoms list. First adding iatm itself.
            call excl_atms%append(iatm)
            do j = 1, size(nbr_atms)
                jatm = nbr_atms(j)
                call excl_atms%append(jatm)
                !Atoms for second ring neighbors
                if (excluded_atoms > 1) then
                    call atat_tab%get_row(jatm, nbr2_atms)
                    do j2 = 1, size(nbr2_atms)
                        j2atm = nbr2_atms(j2)
                        call excl_atms%append(j2atm)
                        !Atoms for third ring neighbors
                        if (excluded_atoms > 2) then
                            call atat_tab%get_row(j2atm, nbr3_atms)
                            do j3 = 1, size(nbr3_atms)
                                j3atm = nbr3_atms(j3)
                                call excl_atms%append(j3atm)
                            end do
                            !Finished adding all third ring neighbors
                        end if
                    end do
                    !Finished adding all second ring neighbors
                end if
                !Finished adding all first ring neighbors
            end do
        end if

        !Sort and remove duplicates
        call excl_atms%unique()

        !Add to exat_tab. Note that iatom appears in the excluded atoms list for
        !iatom.
        do j = 1, excl_atms%len
            jatm = excl_atms%get_val(j)
            call exat_tab%append(iatm, jatm)
        end do
    end do

    !Release additional memory
    call exat_tab%shrink_to_fit()

    end subroutine

!********************************************************************************

end module connectivity_m
