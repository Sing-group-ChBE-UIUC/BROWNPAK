module m_cell_list

use m_precision
use m_vector
use m_globals

implicit none

private

public :: cl_init, cl_delete, cl_print, cl_build, &
    cl_get_num_cells, cl_get_contents, cl_get_nbr_cells 

real(rp), dimension(3)   :: cell_size = 0.0_rp
    !! Cell size along *x*, *y*, & *z*.
integer,  dimension(3)   :: nc_max = 0
    !! Maximum number of cells along *x*, *y*, & *z*.
integer,  dimension(3)   :: nc = 0
    !! Number of cells along *x*, *y*, & *z*.
integer  :: nct_max = 0
    !! Maximum total number of cells.
integer  :: nct = 0
    !! Total of cells.
integer, dimension(:), allocatable, target :: cells
    !! *(na_max,)* array. Listing atoms in each cell.
integer, dimension(:), allocatable :: cells_pos
    !! *(0:nct_max,)* index array. Note: 0-based indexing.
type(ivector_t) :: cell_nbrs
    !! Lists neighbor cells for each cell.
integer, dimension(:), allocatable :: cell_nbrs_pos
    !! *(0:nct_max,)* index array. Note: 0-based indexing.
integer, dimension(:), allocatable :: host_cells
    !! *(na_max,)* array. *host_cells(i)* stores the linear index of the cell
    !! containing atom *i*. *na_max* is the total number of atoms under consideration.
integer, dimension(:), allocatable :: cell_pop
    !! *(0:nct_max-1,)* array storing population of each cell. Note: 0-based indexing.
integer, dimension(3,13), parameter :: d = reshape( [           &
                &                            1, 0, 0,           &
                &    1, 1, 0,   -1, 1, 0,    0, 1, 0,           &
                &    0, 0, 1,   -1, 0, 1,    1, 0, 1,           &
                &   -1,-1, 1,    0,-1, 1,    1,-1, 1,           &
                &   -1, 1, 1,    0, 1, 1,    1, 1, 1 ], [ 3, 13] )

contains

!******************************************************************************

subroutine cl_init()
    !! Allocates memory for cells.

    integer :: na_max
        !! Maximum number of atoms to be handled
    
    na_max = num_atoms

    nc_max(1) = floor( simbox%basis(1,1)/rcutoff )
    nc_max(2) = floor( simbox%basis(2,2)/rcutoff )
    nc_max(3) = floor( simbox%basis(3,3)/rcutoff )

    nct_max = product(nc_max)

    cell_size(1) = simbox%basis(1,1)/nc_max(1)
    cell_size(2) = simbox%basis(2,2)/nc_max(2)
    cell_size(3) = simbox%basis(3,3)/nc_max(3)

    nc = nc_max; nct = nct_max

    allocate(cells(na_max))
    allocate(cells_pos(0:nct_max))

    call ivector_init(cell_nbrs, 13*nct_max)
    allocate(cell_nbrs_pos(0:nct_max))

    allocate(host_cells(na_max))
    allocate(cell_pop(0:nct_max-1))

    !Find cell neighbors for non-deforming box
    if (.not. simbox%is_deforming) call cl_build_cell_nbrs()

    end subroutine

!******************************************************************************

subroutine cl_build_cell_nbrs()

    integer :: j
    integer :: ncx, ncy, ncz
    integer :: icx, icy, icz, ic
    integer :: jcx, jcy, jcz, jc

    ncx = nc(1); ncy = nc(2); ncz = nc(3)

    cell_nbrs_pos = 1

    do icz = 0, ncz-1
        do icy = 0, ncy-1
            do icx = 0, ncx-1
                ic = icz*ncx*ncy + icy*ncx + icx
                cell_nbrs_pos(ic+1) = cell_nbrs_pos(ic)
                do j = 1, 13
                    jcx = modulo( icx + d(1,j), ncx )
                    jcy = modulo( icy + d(2,j), ncy )
                    jcz = modulo( icz + d(3,j), ncz )
                    jc = jcz*ncx*ncy + jcy*ncx + jcx
                    call cell_nbrs%append(jc)
                    cell_nbrs_pos(ic+1) = cell_nbrs_pos(ic+1) + 1
                end do
            end do
        end do
    end do

    end subroutine

!******************************************************************************

subroutine cl_delete()
    !! Deallocates memory allocated in `cl_init`.

    if (allocated(cells))         deallocate(cells)
    if (allocated(cells_pos))     deallocate(cells_pos)

    call cell_nbrs%delete()
    if (allocated(cell_nbrs_pos)) deallocate(cell_nbrs_pos)

    if (allocated(host_cells))    deallocate(host_cells)
    if (allocated(cell_pop))      deallocate(cell_pop)

    cell_size = 0.0_rp
    nc_max = 0; nc = 0
    nct_max = 0; nct = 0

    end subroutine

!******************************************************************************

subroutine cl_build()
    !! Sorts atoms into cells for calculating short-range interations

    real(rp), dimension(3) :: ri
    integer :: i, j, iatm 
    integer :: na, ncx, ncy, ncz
    integer :: icx, icy, icz, ic

    cells = 0
    host_cells = 0
    cell_pop = 0
    ncx = nc(1); ncy = nc(2); ncz = nc(3)

    !Loop over particles and put into cells
    do iatm = 1, num_atoms
        ri = coordinates(:,iatm)
        icx = int( ri(1)/cell_size(1) )
        icy = int( ri(2)/cell_size(2) )
        icz = int( ri(3)/cell_size(3) )

        !Take care of roundoff error
        if ( any( [icx,icy,icz] > (nc-1) ) ) then
            print*, 'i: ', iatm, 'r: ', ri
            print*, 'f: ', forces(:,iatm)
            print*, icx, icy, icz
            print*, 'NTS ', nts
            stop
        end if

        if ( any( [icx,icy,icz] < 0 ) ) then
            print*, 'i: ', iatm, 'r: ', ri
            print*, 'f: ', forces(:,iatm)
            print*, icx, icy, icz
            print*, 'NTS ', nts
            stop
        end if

        !If atoms are exactly on the box edge
        if ( icx > (ncx-1) ) icx = ncx - 1
        if ( icy > (ncy-1) ) icy = ncy - 1
        if ( icz > (ncz-1) ) icz = ncz - 1

        ic = icz*ncx*ncy + icy*ncx + icx
        host_cells(iatm) = ic
        cell_pop(ic) = cell_pop(ic) + 1
    end do

    !Loop over all cells to set up pointers
    cells_pos(0) = 1
    do ic = 0, (nct-1)
        cells_pos(ic+1) = cells_pos(ic) + cell_pop(ic)
    end do

    !Loop over all atoms
    do iatm = 1, num_atoms
        ic = host_cells(iatm)
        j = cells_pos(ic)
        cells(j) = iatm
        cells_pos(ic) = cells_pos(ic) + 1
    end do

    !Loop over all cells to set up pointers
    cells_pos(0) = 1
    do ic = 0, (nct-1)
        cells_pos(ic+1) = cells_pos(ic) + cell_pop(ic)
    end do

    !call cl_print()
    !stop
    end subroutine

!******************************************************************************

function cl_get_num_cells() result (res)
    !! Returns the total number of cells

    integer :: res

    res = nct

    end function

!******************************************************************************

subroutine cl_get_contents(ic, res)
    !! Returns a pointer to the entries of cell with linear index *ic*.

    integer, intent(in) :: ic
    integer, dimension(:), pointer, intent(out) :: res
    integer :: ibeg, iend

    res => null()
    ibeg = cells_pos(ic); iend = cells_pos(ic+1) - 1
    res => cells(ibeg:iend)

    end subroutine

!******************************************************************************

subroutine cl_get_nbr_cells(ic, res)
    !! Returns a pointer to the neighbor cells of cell with linear index *ic*.

    integer, intent(in) :: ic
    integer, dimension(:), pointer, intent(out) :: res
    integer :: ibeg, iend

    res => null()
    ibeg = cell_nbrs_pos(ic); iend = cell_nbrs_pos(ic+1) - 1
    call cell_nbrs%get_data(res, ibeg, iend)

    end subroutine

!********************************************************************************

subroutine cl_print()
    !! Prints a cell list

    integer, dimension(:), pointer :: aic => null()
    integer, dimension(:), pointer :: nbrc => null()
    integer :: na, ncx, ncy, ncz
    integer :: icx, icy, icz, ic

    ncx = nc(1); ncy = nc(2); ncz = nc(3)

    write(*,'("ncx: ", i0, " ncy: ", i0, " ncz: ", i0)') ncx, ncy, ncz
    write(*,'("lcx: ", g0.6, " lcy: ", g0.6, " lcz: ", g0.6)') cell_size
    write(*, *) 'CELL CONTENTS'
    do icz = 0, ncz-1
        do icy = 0, ncy-1
            do icx = 0, ncx-1
                ic = icz*ncx*ncy + icy*ncx + icx
                write(*,'("(",i0,",",i0,",",i0,")[",i0,"] ")', advance='no') icx, icy, icz, ic
                call cl_get_contents(ic, aic)
                write(*, *) size(aic), aic
            end do
        end do
    end do

    write(*, *)
    write(*, *) 'NBR CELLS'
    do icz = 0, ncz-1
        do icy = 0, ncy-1
            do icx = 0, ncx-1
                ic = icz*ncx*ncy + icy*ncx + icx
                write(*,'("(",i0,",",i0,",",i0,")[",i0,"] ")', advance='no') icx, icy, icz, ic
                call cl_get_nbr_cells(ic, nbrc)
                write(*, *) size(nbrc), nbrc
            end do
        end do
    end do

    end subroutine

!********************************************************************************

end module m_cell_list
