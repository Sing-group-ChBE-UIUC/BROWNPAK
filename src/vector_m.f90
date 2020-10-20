module vector_m
    !! Implements a vector as an extendable array.

use constants_m
use qsort_m

implicit none

type ivector_t
    integer :: len_init = 8 !Must be > 0
    integer :: len = 0
    integer :: len_max = 16 !Must be > len_init
    integer, dimension(:), allocatable :: buffer

    contains
        procedure :: delete => ivector_delete
        procedure :: clear => ivector_clear
        procedure :: get_val => ivector_get_val
        procedure :: set_val => ivector_set_val
        procedure :: get_data => ivector_get_data
        procedure :: append => ivector_append
        procedure :: pop => ivector_pop
        procedure :: resize => ivector_resize
        procedure :: shrink_to_fit => ivector_shrink_to_fit
        procedure :: sort => ivector_sort
        procedure :: unique => ivector_unique
        procedure :: print => ivector_print
end type ivector_t

type dvector_t
    integer :: len_init = 8 !Must be > 0
    integer :: len = 0
    integer :: len_max = 16 !Must be > len_init
    real(rp), dimension(:), allocatable :: buffer

    contains
        procedure :: delete => dvector_delete
        procedure :: clear => dvector_clear
        procedure :: get_val => dvector_get_val
        procedure :: set_val => dvector_set_val
        procedure :: get_data => dvector_get_data
        procedure :: append => dvector_append
        procedure :: pop => dvector_pop
        procedure :: resize => dvector_resize
        procedure :: shrink_to_fit => dvector_shrink_to_fit
        procedure :: sort => dvector_sort
        procedure :: unique => dvector_unique
        procedure :: print => dvector_print
end type dvector_t

interface ivector_init
    module procedure ivector_init
    module procedure ivector_from_array
end interface

interface dvector_init
    module procedure dvector_init
    module procedure dvector_from_array
end interface

interface assignment(=)
    module procedure ivector_assign
    module procedure dvector_assign
    module procedure i_dvector_assign
end interface


contains

!******************************************************************************

!> Creates an empty *ivector* with all elements equal to zero.
subroutine ivector_init(this, len_init, ierr)

    type(ivector_t), intent(in out) :: this
    integer, intent(in), optional :: len_init
        !! Must be > 0
    integer, intent(out), optional :: ierr
    integer :: istat

    if (present(ierr)) ierr = 0

    if (present(len_init)) then
        if (len_init > 0) this%len_init = len_init
        this%len_max = 2*len_init
    end if

    allocate(this%buffer(this%len_max), stat=istat)
    if (istat /= 0) then
        write(*, *) 'error: allocation failure for this%buffer of size', this%len_max
        if (present(ierr)) ierr = 1
        return
    end if

    this%buffer = 0

    end subroutine

!******************************************************************************

!> Creates an empty *dvector* with all elements equal to zero.
subroutine dvector_init(this, len_init, ierr)

    type(dvector_t), intent(in out) :: this
    integer, intent(in), optional :: len_init
        !! Must be > 0
    integer, intent(out), optional :: ierr
    integer :: istat

    if (present(ierr)) ierr = 0

    if (present(len_init)) then
        if (len_init > 0) this%len_init = len_init
        this%len_max = 2*len_init
    end if

    allocate(this%buffer(this%len_max), stat=istat)
    if (istat /= 0) then
        write(*, *) 'error: allocation failure for this%buffer of size', this%len_max
        if (present(ierr)) ierr = 1
        return
    end if

    this%buffer = 0.0_rp

    end subroutine

!******************************************************************************

!> Creates an *ivector* with all elements from an array
subroutine ivector_from_array(this, x)

    type(ivector_t), intent(in out) :: this
    integer, dimension(:), intent(in) :: x
    integer, dimension(:), allocatable :: temp
    integer :: n

    call ivector_init(this)
    n = size(x,1)

    if (this%len_max < n) then
        allocate(temp(2*n))
        temp = 0
        temp(1:n) = x
        call move_alloc(temp, this%buffer)
        this%len = n
        this%len_max = 2*n
    else
        this%buffer(1:n) = x
        this%len = n
    end if

    end subroutine

!******************************************************************************

!> Creates a *dvector* with all elements from an array
subroutine dvector_from_array(this, x)

    type(dvector_t), intent(in out) :: this
    real(rp), dimension(:), intent(in) :: x
    real(rp), dimension(:), allocatable :: temp
    integer :: n

    call dvector_init(this)
    n = size(x,1)

    if (this%len_max < n) then
        allocate(temp(2*n))
        temp = 0
        temp(1:n) = x
        call move_alloc(temp, this%buffer)
        this%len = n
        this%len_max = 2*n
    else
        this%buffer(1:n) = x
        this%len = n
    end if

    end subroutine

!******************************************************************************

!> Deletes an *ivector*. No access is allowed to this object after this call.
subroutine ivector_delete(this)

    class(ivector_t), intent(in out) :: this

    if (allocated(this%buffer)) deallocate(this%buffer)

    this%len = 0

    end subroutine

!******************************************************************************

!> Deletes a *dvector*. No access is allowed to this object after this call.
subroutine dvector_delete(this)

    class(dvector_t), intent(in out) :: this

    if (allocated(this%buffer)) deallocate(this%buffer)

    this%len = 0

    end subroutine

!******************************************************************************

!> Clears an *ivector*. Access allowed after a call to clear.
subroutine ivector_clear(this)

    class(ivector_t), intent(in out) :: this

    this%len = 0
    this%buffer = 0

    end subroutine

!******************************************************************************

!> Clears a *dvector*. Access allowed after a call to clear.
subroutine dvector_clear(this)

    class(dvector_t), intent(in out) :: this

    this%len = 0
    this%buffer = 0.0_rp

    end subroutine

!******************************************************************************

!> Copies the contents of *ivector* `other` to *ivector* `this`
!!
subroutine ivector_assign(this, other)

    class(ivector_t), intent(in out) :: this
    class(ivector_t), intent(in) :: other
    integer :: n

    n = other%len

    if (this%len >= n) then
        !Copy contents if buffer is large enough
        this%buffer(1:n) = other%buffer
        this%len = n
    else
        !Need to enlarge buffer. Reallocate and copy
        deallocate(this%buffer)
        allocate(this%buffer, source=other%buffer(1:n))
        this%len = n
        this%len_max = n
    end if

    end subroutine

!******************************************************************************

!> Copies the contents of *dvector* `other` to *dvector* `this`
!!
subroutine dvector_assign(this, other)

    class(dvector_t), intent(in out) :: this
    class(dvector_t), intent(in) :: other
    integer :: n

    n = other%len

    if (this%len >= n) then
        !Copy contents if buffer is large enough
        this%buffer(1:n) = other%buffer(1:n)
        this%len = n
    else
        !Need to enlarge buffer. Reallocate and copy
        deallocate(this%buffer)
        allocate(this%buffer, source=other%buffer(1:n))
        this%len = n
        this%len_max = n
    end if

    end subroutine

!******************************************************************************

subroutine i_dvector_assign(this, other)
    !!  Copies the contents of *ivector* `other` to *dvector* `this`. Integers
    !!  are cast to reals.
    !""

    class(dvector_t), intent(in out) :: this
    class(ivector_t), intent(in) :: other
    integer :: n

    n = other%len

    if (this%len >= n) then
        !Copy contents if buffer is large enough
        this%buffer(1:n) = dble(other%buffer(1:n))
        this%len = n
    else
        !Need to enlarge buffer. Reallocate and copy
        deallocate(this%buffer)
        allocate(this%buffer(n))
        this%buffer = dble(other%buffer(1:n))
        this%len = n
        this%len_max = n
    end if

    end subroutine

!******************************************************************************

!> Returns the length of an *ivector*
function ivector_get_len(this) result(res)

    class(ivector_t), intent(in) :: this
    integer :: res

    res = this%len

    end function

!******************************************************************************

!> Returns the length of a *dvector*
function dvector_get_len(this) result(res)

    class(dvector_t), intent(in) :: this
    integer :: res

    res = this%len

    end function

!******************************************************************************

!> Returns the ith element of an *ivector*.
function ivector_get_val(this, i) result(res)

    class(ivector_t), intent(in) :: this
    integer, intent(in) :: i
    integer :: res

    if (i > this%len) then
        write(*,*) 'error: out-of-bounds index ', i
        write(*,*) 'for ubound', this%len
        stop
    else
        res = this%buffer(i)
    end if

    end function

!******************************************************************************

!> Returns the ith element of a *dvector*.
function dvector_get_val(this, i) result(res)

    class(dvector_t), intent(in) :: this
    integer, intent(in) :: i
    real(rp) :: res

    if (i > this%len) then
        write(*,*) 'error: out-of-bounds index ', i
        write(*,*) 'for ubound', this%len
        stop
    else
        res = this%buffer(i)
    end if

    end function

!******************************************************************************

!> Sets the value of the ith element of an *ivector*. No bounds check is performed.
subroutine ivector_set_val(this, i, val)

    class(ivector_t), intent(in out) :: this
    integer, intent(in) :: i
    integer, intent(in) :: val

    if (i > this%len) then
        write(*,*) 'error: out-of-bounds index ', i
        write(*,*) 'for ubound', this%len
        stop
    else
        this%buffer(i) = val
    end if

    end subroutine

!******************************************************************************

!> Sets the value of the ith element of a *dvector*. No bounds check is performed.
subroutine dvector_set_val(this, i, val)

    class(dvector_t), intent(in out) :: this
    integer, intent(in) :: i
    real(rp), intent(in) :: val

    if (i > this%len) then
        write(*,*) 'error: out-of-bounds index ', i
        write(*,*) 'for ubound', this%len
        stop
    else
        this%buffer(i) = val
    end if

    end subroutine

!******************************************************************************

!> Adds an element to the end of an *ivector*. Reallocation will take place if required.
subroutine ivector_append(this, val)

    class(ivector_t), intent(in out) :: this
    integer, intent(in) :: val
    integer, dimension(:), allocatable :: temp

    !Double len_max if current length is equal to len_max
    if (this%len == this%len_max) then

        allocate(temp(2*this%len_max))
        temp = 0
        temp(1:this%len) = this%buffer(1:this%len)

        call move_alloc(temp, this%buffer)
        this%len_max = 2*this%len_max

    end if

    this%buffer(this%len+1) = val
    this%len = this%len + 1

    end subroutine

!******************************************************************************

!> Adds an element to the end of a *dvector*. Reallocation will take place if required.
subroutine dvector_append(this, val)

    class(dvector_t), intent(in out) :: this
    real(rp), intent(in) :: val
    real(rp), dimension(:), allocatable :: temp

    !Double len_max if current length is equal to len_max
    if (this%len == this%len_max) then

        allocate(temp(2*this%len_max))
        temp = 0
        temp(1:this%len) = this%buffer(1:this%len)

        call move_alloc(temp, this%buffer)
        this%len_max = 2*this%len_max

    end if

    this%buffer(this%len+1) = val
    this%len = this%len + 1

    end subroutine

!******************************************************************************

!> Removes the last element and returns it. Calling this method on an empty list
! will generate an error.
function ivector_pop(this) result(val)

    class(ivector_t), intent(in out) :: this
    integer :: val

    val = this%buffer(this%len)
    this%len = this%len - 1

    end function

!******************************************************************************

!> Removes the last element and returns it. Calling this method on an empty list
! will generate an error.
function dvector_pop(this) result(val)

    class(dvector_t), intent(in out) :: this
    real(rp) :: val

    val = this%buffer(this%len)
    this%len = this%len - 1

    end function

!******************************************************************************

!> Returns a pointer to the underlying data of an *ivector*
! No bounds checking is performed
subroutine ivector_get_data(this, res, ibeg, iend)

    class(ivector_t), target, intent(in) :: this
    integer, dimension(:), pointer, intent(out) :: res
    integer, intent(in), optional :: ibeg
    integer, intent(in), optional :: iend
    integer :: ibeg_
    integer :: iend_

    res => null()
    ibeg_ = 1
    iend_ = this%len

    if (present(ibeg)) then
        ibeg_ = ibeg
    end if

    if (present(iend)) then
        iend_ = iend
    end if

    res => this%buffer(ibeg_:iend_)

    end subroutine

!******************************************************************************

!> Returns a pointer to the underlying data of a *dvector*
! No bounds checking is performed
subroutine dvector_get_data(this, res, ibeg, iend)

    class(dvector_t), target, intent(in) :: this
    real(rp), dimension(:), pointer, intent(out) :: res
    integer, intent(in), optional :: ibeg
    integer, intent(in), optional :: iend
    integer :: ibeg_
    integer :: iend_

    res => null()
    ibeg_ = 1
    iend_ = this%len

    if (present(ibeg)) then
        ibeg_ = ibeg
    end if

    if (present(iend)) then
        iend_ = iend
    end if

    res => this%buffer(ibeg_:iend_)

    end subroutine

!******************************************************************************

!> Resizes a vector to a given size. Existing data is truncated if desired size
!> is smaller than current size. Otherwise, the empty spaces are filled with
!> zero.
subroutine ivector_resize(this, new_size, ierr)

    class(ivector_t), intent(in out) :: this
    integer, intent(in) :: new_size
    integer, intent(out), optional :: ierr
    integer, dimension(:), allocatable :: temp
    integer :: istat

    if (this%len_max >= new_size) then
        this%buffer(this%len+1:new_size) = 0
        this%len = new_size
    else
       allocate(temp(new_size), stat=istat) 
       if (istat /= 0) then
            write(*,*) 'error: allocation failure for `temp` from `this%buffer`.'
            if (present(ierr)) ierr = 1
            return
        end if
        temp = 0; temp(1:this%len_max) = this%buffer
        call move_alloc(temp, this%buffer)
        this%len_max = new_size; this%len = new_size
    end if

    if (present(ierr)) ierr = 0

    end subroutine

!******************************************************************************

!> Resizes a vector to a given size. Existing data is truncated if desired size
!> is smaller than current size. Otherwise, the empty spaces are filled with
!> zero.
subroutine dvector_resize(this, new_size, ierr)

    class(dvector_t), intent(in out) :: this
    integer, intent(in) :: new_size
    integer, intent(out), optional :: ierr
    real(rp), dimension(:), allocatable :: temp
    integer :: istat

    ierr = 0
    if (this%len_max > new_size) then
        this%buffer(this%len+1:new_size) = 0.0_rp
        this%len = new_size
    else
       allocate(temp(new_size), stat=istat) 
       if (istat /= 0) then
            write(*,*) 'error: allocation failure for `temp` from `this%buffer`.'
            if (present(ierr)) ierr = 1
            return
        end if
        temp = 0.0_rp; temp(1:this%len_max) = this%buffer
        call move_alloc(temp, this%buffer)
        this%len_max = new_size; this%len = new_size
    end if

    if (present(ierr)) ierr = 0

    end subroutine

!******************************************************************************

!> Releases additional memory to fit underlying data of an *ivector* to a size
!> this%len_init.
subroutine ivector_shrink_to_fit(this, ierr)

    class(ivector_t), intent(in out) :: this
    integer, intent(out), optional :: ierr
    integer, dimension(:), allocatable :: temp
    integer :: istat

    if (present(ierr)) ierr = 0

    if (this%len > this%len_init) then
        allocate(temp, source=this%buffer(1:this%len), stat=istat)
    else
        allocate(temp(this%len_init), stat=istat)
        temp = 0
        temp(1:this%len) = this%buffer(1:this%len)
    end if

    if (istat /= 0) then
        write(*,*) 'error: allocation failure for `temp` from `this%buffer`.'
        if (present(ierr)) ierr = 1
        return
    end if

    call move_alloc(temp, this%buffer)

    this%len_max = size(this%buffer)

    end subroutine

!******************************************************************************

!> Releases additional memory to fit underlying data of a *dvector* to a size
!> this%len_init.
subroutine dvector_shrink_to_fit(this, ierr)

    class(dvector_t), intent(in out) :: this
    integer, intent(out), optional :: ierr
    real(rp), dimension(:), allocatable :: temp
    integer :: istat

    if (present(ierr)) ierr = 0

    if (this%len > this%len_init) then
        allocate(temp, source=this%buffer(1:this%len), stat=istat)
    else
        allocate(temp(this%len_init), stat=istat)
        temp = 0.0_rp
        temp(1:this%len) = this%buffer(1:this%len)
    end if

    if (istat /= 0) then
        write(*,*) 'error: allocation failure for `temp` from `this%buffer`.'
        if (present(ierr)) ierr = 1
        return
    end if

    call move_alloc(temp, this%buffer)

    this%len_max = size(this%buffer)

    end subroutine

!******************************************************************************

!>  Sorts an *ivector* in ascending order. If `order` is provided, it will contain
!!  the sorted indices.  The size of `order` must be at least `this%len`.
subroutine ivector_sort(this, order)

    class(ivector_t), intent(in out) :: this
    integer, dimension(:), intent(out), optional :: order
    integer :: n

    n = this%len
    if (n > 1) then
        if (present(order)) then
            call iqsort(this%buffer(1:n), order)
        else
            call iqsort(this%buffer(1:n))
        end if
    end if

    end subroutine

!******************************************************************************

!>  Sorts a *dvector* in ascending order. If `order` is provided, it will contain
!!  the sorted indices.  The size of `order` must be at least `this%len`.
subroutine dvector_sort(this, order)

    class(dvector_t), intent(in out) :: this
    integer, dimension(:), intent(out), optional :: order
    integer :: n

    n = this%len
    if (n > 1) then
        if (present(order)) then
            call dqsort(this%buffer(1:n), order)
        else
            call dqsort(this%buffer(1:n))
        end if
    end if

    end subroutine

!******************************************************************************

!> Sorts and removes all duplicate entries of an *ivector*. Note that the internal
!! buffer size remains unchanged. To reduce the buffer size, call
!! `ivector_shrink_to_fit`.
subroutine ivector_unique(this)

    class(ivector_t), intent(in out) :: this
    integer, dimension(this%len) :: temp
    integer :: n
    integer :: i

    if (this%len < 2) return !Already sorted

    temp = this%buffer(1:this%len)

    call iqsort(temp)

    !The first element is unique by default
    this%buffer(1) = temp(1)
    n = 1
    !Successively copy non-duplicate elements
    do i = 2, this%len
        if (temp(i) /= this%buffer(n)) then
            this%buffer(n+1) = temp(i)
            n = n + 1
        end if
    end do

    !Set the length to number of unique elements
    this%len = n

    end subroutine

!******************************************************************************

!> Sorts and removes all duplicate entries of a *dvector*. Note that the internal
!! buffer size remains unchanged. To reduce the buffer size, call
!! `dvector_shrink_to_fit`.
subroutine dvector_unique(this)

    class(dvector_t), intent(in out) :: this
    real(rp), dimension(this%len) :: temp
    integer :: n
    integer :: i

    if (this%len < 2) return !Already sorted

    temp = this%buffer

    call dqsort(temp)

    !The first element is unique by default
    this%buffer(1) = temp(1)
    n = 1
    !Successively copy non-duplicate elements
    do i = 2, this%len
        if (temp(i) /= this%buffer(n)) then
            this%buffer(n+1) = temp(i)
            n = n + 1
        end if
    end do

    !Set the length to number of unique elements
    this%len = n

    end subroutine

!******************************************************************************

!> Prints an *ivector*
subroutine ivector_print(this)

    class(ivector_t), intent(in) :: this
    integer :: i

    do i = 1, this%len
        write(*,'(i0,1x,i0)') i, this%buffer(i)
    end do
    
    end subroutine

!******************************************************************************

!> Prints a *dvector*
subroutine dvector_print(this)

    class(dvector_t), intent(in) :: this
    integer :: i

    do i = 1, this%len
        write(*,'(i0,1x,f0.15)') i, this%buffer(i)
    end do
    
    end subroutine

!******************************************************************************

end module vector_m
