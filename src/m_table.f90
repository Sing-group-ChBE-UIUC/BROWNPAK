module m_table
    !! Implements a table with contiguously stored rows.

use m_precision
use m_vector

implicit none

private
public :: itable_t
public :: itbl_init

type itable_t
    integer :: num_rows = 0 !Must be > 0
    type(ivector_t) :: buffer
    integer, dimension(:), allocatable :: row_indx

    contains
        procedure :: delete   => itbl_delete
        procedure :: clear    => itbl_clear
        procedure :: append   => itbl_append
        procedure :: set_val  => itbl_set_val
        procedure :: is_in    => itbl_is_in
        procedure :: get_val  => itbl_get_val
        procedure :: get_row  => itbl_get_row
        procedure :: shrink_to_fit  => itbl_shrink_to_fit
        procedure :: print    => itbl_print
end type itable_t

contains

!******************************************************************************

subroutine itbl_init(this, num_rows, ierr)
    !! Creates an empty *itable_t* with *num_rows* rows and all rows having
    !! zero elements.

    type(itable_t), intent(in out) :: this
    integer, intent(in)             :: num_rows
        !! Must be > 0
    integer, intent(out), optional :: ierr
    integer :: istat

    if (present(ierr)) ierr = 0

    !Set initial size guess to be num_rows
    call ivector_init(this%buffer, num_rows)

    allocate(this%row_indx(num_rows+1), stat=istat)
    if (istat /= 0) then
        write(*, *) 'error: allocation failure for this%row_indx of size', num_rows
        if (present(ierr)) ierr = 1
        return
    end if

    !Set the number of rows
    this%num_rows = num_rows

    !Empty table: All rows have zero elements
    this%row_indx = 1

    end subroutine

!******************************************************************************

subroutine itbl_delete(this)
    !! Deletes an *itable_t*. No access is allowed to this object after this call.

    class(itable_t), intent(in out) :: this

    call this%buffer%delete()
    if (allocated(this%row_indx)) deallocate(this%row_indx)
    this%num_rows = 0

    end subroutine

!******************************************************************************

subroutine itbl_clear(this)
    !! Clears all rows. Does not deallocate memory.

    class(itable_t), intent(in out) :: this

    call this%buffer%clear()
    this%row_indx = 1

    end subroutine

!******************************************************************************

subroutine itbl_append(this, irow, val)
    !! Appends an element `val` to row *irow*.

    class(itable_t), intent(in out) :: this
    integer, intent(in) :: irow
    integer, intent(in) :: val

    call this%buffer%append(val)

    this%row_indx(irow+1:) = this%row_indx(irow+1:) + 1

    end subroutine

!******************************************************************************

subroutine itbl_set_val(this, irow, j, val)
    !! Sets the value of the *j*th element of row *irow*.

    class(itable_t), intent(in out) :: this
    integer, intent(in) :: irow
    integer, intent(in) :: j
    integer, intent(in) :: val
    integer :: k, n

    n = this%row_indx(irow+1) - this%row_indx(irow)
    if (j > n) then
        write(*,*) 'error: out-of-bounds index ', j
        write(*,*) 'for row length', n
        stop
    else
        k = this%row_indx(irow) + j - 1
        call this%buffer%set_val(k, val)
    end if

    end subroutine

!******************************************************************************

function itbl_is_in(this, irow, val) result(res)
    !! Returns .true. if *val* is in row *irow*, .false. otherwise.

    class(itable_t), intent(in) :: this
    integer, intent(in) :: irow
    integer, intent(in) :: val
    logical :: res
    integer :: i

    res = .false.

    do i = this%row_indx(irow), this%row_indx(irow+1)-1
        if (this%buffer%get_val(i) == val) then
            res = .true.
            exit
        end if
    end do
    
    end function

!******************************************************************************

function itbl_get_val(this, irow, j) result(res)
    !! Returns the *j*th element of row *irow*

    class(itable_t), intent(in) :: this
    integer, intent(in) :: irow
    integer, intent(in) :: j
    integer :: res
    integer :: k, n

    n = this%row_indx(irow+1) - this%row_indx(irow)
    if (j > n) then
        write(*,*) 'error: out-of-bounds index ', j
        write(*,*) 'for row length', n
        stop
    else
        k = this%row_indx(irow) + j - 1
        res = this%buffer%get_val(k)
    end if
    
    end function

!******************************************************************************

subroutine itbl_get_row(this, irow, res)
    !! Returns a pointer to the row data of *irow*. No bounds checking is performed.

    class(itable_t), target, intent(in) :: this
    integer, intent(in) :: irow
    integer, dimension(:), pointer, intent(out) :: res
    integer :: ibeg, iend

    res => null()
    ibeg = this%row_indx(irow)
    iend = this%row_indx(irow+1) - 1

    call this%buffer%get_data(res, ibeg, iend)

    end subroutine

!******************************************************************************

subroutine itbl_shrink_to_fit(this, ierr)
    !! Releases additional memory to fit underlying data.

    class(itable_t), intent(in out) :: this
    integer, intent(out), optional  :: ierr

    if (present(ierr)) then
        call this%buffer%shrink_to_fit(ierr)
    else
        call this%buffer%shrink_to_fit()
    end if

    end subroutine

!******************************************************************************

subroutine itbl_print(this)
    !! Prints an *itable_t*.

    class(itable_t), intent(in) :: this
    integer :: i, j

    do i = 1, this%num_rows
        write(*,'(i0,": ",2x)', advance='no') i
        do j = this%row_indx(i), this%row_indx(i+1)-1
            write(*, '(i0,2x)', advance='no') this%buffer%get_val(j)
        end do
        write(*,*)
    end do
    
    end subroutine

!******************************************************************************

end module m_table
