module m_verlet
    !! Routines for building Verlet neighbor table.

use m_precision
use m_vector
use m_table
use m_globals

implicit none

private

public :: verlet_init, verlet_build, verlet_delete
public :: verlet_tab

type(itable_t) :: verlet_tab
    !! Verlet table
real(rp), dimension(:,:), allocatable :: coordinates_save
    !!  (3, *num_atoms*) array
real(rp), dimension(:,:), allocatable :: coordinates_dr
    !!  (3, *num_atoms*) array

real(rp) :: rskin_sq = 0.0_rp
real(rp) :: tskin_sq = 0.0_rp

contains

!******************************************************************************

subroutine verlet_init(rskin, tskin)

    real(rp), intent(in) :: rskin
    real(rp), intent(in) :: tskin

    rskin_sq = rskin**2
    tskin_sq = tskin**2

    call itbl_init(verlet_tab, num_atoms-1)

    allocate( coordinates_save(3,num_atoms) )
    allocate( coordinates_dr(3,num_atoms)   )

    !Initializing to zero
    coordinates_save = 0.0_rp
    coordinates_dr = 0.0_rp

    end subroutine

!******************************************************************************

subroutine verlet_delete()

    if (allocated(coordinates_save)) deallocate(coordinates_save)
    if (allocated(coordinates_dr)  ) deallocate(coordinates_dr)
    call verlet_tab%delete()

    end subroutine

!******************************************************************************

subroutine verlet_build()

    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rj
    real(rp), dimension(3) :: rij
    real(rp) :: dr_sq_max
    real(rp) :: rij_sq
    logical, save :: first_call = .true.
    integer :: i, j

    !On first call no check for rebuilding
    if (.not. first_call) then
        !Check whether rebuilding the list is necessary
        coordinates_dr = coordinates(:,1:num_atoms) - coordinates_save
        if (imcon /= 0) then
            do i = 1, num_atoms
                call simbox%get_image(coordinates_dr(:,i))
            end do
        end if
        dr_sq_max = maxval(sum(coordinates_dr**2, dim=1))
        if ( 4*dr_sq_max < tskin_sq ) then
            return
        end if
    end if

    first_call = .false.

    !Clear table
    call verlet_tab%clear()

    !Loop over all pairs to build list
    do i = 1, num_atoms-1
        ri = coordinates(:,i)
        do j = i+1, num_atoms
            rj = coordinates(:,j)
            rij = rj - ri
            if (imcon /= 0) call simbox%get_image(rij)
            rij_sq = sum(rij**2)
            if (rij_sq < rskin_sq) call verlet_tab%append(i,j)
        end do
    end do

    !Release additional memory
    call verlet_tab%shrink_to_fit()

    !Back up positions
    coordinates_save = coordinates(:,1:num_atoms)

    end subroutine

!********************************************************************************

end module m_verlet
