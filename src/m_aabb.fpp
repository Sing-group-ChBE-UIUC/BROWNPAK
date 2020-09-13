#:include 'asserts.fypp'

module m_aabb

use m_precision
use m_strings

implicit none

type aabb_t
    real(rp), dimension(3) :: lbnd
    real(rp), dimension(3) :: ubnd
    real(rp), dimension(3) :: center
    real(rp) :: srfarea
    contains
        procedure :: init => aabb_init
        procedure :: print => aabb_print
        procedure :: clear => aabb_clear
        procedure :: get_extent => aabb_get_extent
        procedure :: update => aabb_update
        procedure :: fatten => aabb_fatten
        procedure :: includes => aabb_includes
        procedure :: overlaps => aabb_overlaps
        procedure, private :: calc_center => aabb_calc_center
        procedure, private :: calc_srfarea => aabb_calc_srfarea
end type aabb_t

interface operator(+)
    module procedure :: aabb_union
end interface

contains

!*******************************************************************************

subroutine aabb_init(this, lbnd, ubnd)
    !! Creates an *aabb_t* instance from lower and upper bounds.

    class(aabb_t), intent(out) :: this
    !! An *aabb_t* instance.
    real(rp), dimension(3), intent(in) :: lbnd
    !! Lower bound
    real(rp), dimension(3), intent(in) :: ubnd
    !! Upper bound

    @:ASSERT(all(lbnd <= ubnd))

    this%lbnd = lbnd; this%ubnd = ubnd
    call this%calc_center()
    call this%calc_srfarea()

    end subroutine

!*******************************************************************************

subroutine aabb_print(this, frmt, str)
    !! Prints an AABB.

    class(aabb_t), intent(in) :: this
        !! An *aabb_t* instance.
    character(len=*), intent(in), optional :: frmt
        !! Fortran-style format string for a real number. Default: *(g0.6)*.
    character(len=:), allocatable, intent(out), optional :: str
        !! If present, the output is printed to this string instead of STDOUT.
    character(len=:), allocatable :: frmt_, buf

    frmt_ = '(g0.6)'
    if (present(frmt)) frmt_ = trim(adjustl(frmt))

    buf = 'lbnd: ['//str_from_num(this%lbnd(1), frmt_)// ', ' &
        &//str_from_num(this%lbnd(2), frmt_) // ', ' &
        &//str_from_num(this%lbnd(3), frmt_) // ']'

    buf = buf // ', ubnd: ['//str_from_num(this%ubnd(1), frmt_)// ', ' &
        &//str_from_num(this%ubnd(2), frmt_) // ', ' &
        &//str_from_num(this%ubnd(3), frmt_) // ']'

    buf = buf // ', center: ['//str_from_num(this%center(1), frmt_)// ', ' &
        &//str_from_num(this%center(2), frmt_) // ', '&
        &//str_from_num(this%center(3), frmt_) // ']'

    buf = buf // ', srfarea: ' // str_from_num(this%srfarea, frmt_)

    if (present(str)) then
        str = buf
    else
        write(*,*) buf
    end if

    end subroutine

!*******************************************************************************

subroutine aabb_clear(this)
    !! Clears all attributes of an AABB and sets them to zero.

    class(aabb_t), intent(in out) :: this
        !! An *aabb_t* instance.

    this%lbnd = 0.0_rp; this%ubnd = 0.0_rp
    this%center = 0.0_rp; this%srfarea = 0.0_rp

    end subroutine

!*******************************************************************************

subroutine aabb_get_extent(this, extent)
    !! Calculates the extent of an *aabb*. The *extent* of an AABB is defined as
    !! the difference between its upper and lower bounds.

    class(aabb_t), intent(in) :: this
        !! An *aabb_t* instance.
    real(rp), dimension(3), intent(out) :: extent
        !! Extent of an AABB.

    extent = this%ubnd - this%lbnd

    end subroutine

!*******************************************************************************

subroutine aabb_update(this, lbnd, ubnd)
    !! Updates an AABB with new bounds.

    class(aabb_t), intent(in out) :: this
        !! An *aabb_t* instance.
    real(rp), dimension(3), intent(in), optional :: lbnd
        !! Lower bound
    real(rp), dimension(3), intent(in), optional :: ubnd
        !! Upper bound

    !Nothing to do if no bounds are present
    if ( (.not. present(lbnd)) .and. (.not. present(ubnd)) ) return

    !Update bounds based on input
    if (present(lbnd)) this%lbnd = lbnd 
    if (present(ubnd)) this%ubnd = ubnd 

    @:ASSERT( all(this%lbnd <= this%ubnd) )

    call this%calc_center()
    call this%calc_srfarea()

    end subroutine

!*******************************************************************************

subroutine aabb_fatten(this, frac)
    !! Fattens an AABB by a fraction of its base extent.

    class(aabb_t), intent(in out) :: this
        !! An *aabb_t* instance.
    real(rp), intent(in) :: frac
        !! Fraction of AABB base extent.
    real(rp), dimension(3) :: extent

    extent = this%ubnd - this%lbnd
    !New bounds
    this%lbnd = this%lbnd - frac*extent; this%ubnd = this%ubnd + frac*extent

    @:ASSERT( all(this%lbnd <= this%ubnd) )

    call this%calc_center()
    call this%calc_srfarea()

    end subroutine

!*******************************************************************************

function aabb_includes(this, other) result(res)
    !! Returns *true* if *this* includes *other*, *false* otherwise. Inclusion
    !! is considered in a strict sense.

    class(aabb_t), intent(in) :: this
    !! An *aabb_t* instance.
    type(aabb_t), intent(in) :: other
    !! An *aabb_t* instance.
    logical :: res

    res = ( all(this%lbnd < other%lbnd) .and. all(this%ubnd > other%ubnd) )

    end function

!*******************************************************************************

function aabb_overlaps(this, other) result(res)
    !! Returns *true* if *this* overlaps *other*, *false* otherwise. Touching
    !! does not count as an overlap.

    class(aabb_t), intent(in) :: this
    !! An *aabb_t* instance.
    type(aabb_t), intent(in) :: other
    !! An *aabb_t* instance.
    logical :: res

    res = ( all(this%ubnd > other%lbnd) .and. all(other%ubnd > this%lbnd) )

    end function

!*******************************************************************************

function aabb_union(x, y) result(z)
    !! Combines AABBs *x* and *y* to return a new AABB *z*.

    type(aabb_t), intent(in) :: x
    !! An *aabb_t* instance.
    type(aabb_t), intent(in) :: y
    !! An *aabb_t* instance.
    type(aabb_t) :: z

    z%lbnd = merge( x%lbnd, y%lbnd, (x%lbnd < y%lbnd) )
    z%ubnd = merge( x%ubnd, y%ubnd, (x%ubnd > y%ubnd) )

    call z%calc_center()
    call z%calc_srfarea()

    end function

!*******************************************************************************

subroutine aabb_calc_center(this)
    !! Calculates the center of an *aabb*.

    class(aabb_t), intent(in out) :: this
    !! An *aabb_t* instance.

    this%center = 0.5_rp*(this%ubnd + this%lbnd)

    end subroutine

!*******************************************************************************

subroutine aabb_calc_srfarea(this)
    !! Calculates the surface area of an *aabb*.

    class(aabb_t), intent(in out) :: this
    !! An *aabb_t* instance.
    real(rp) :: dx, dy, dz

    dx = this%ubnd(1) - this%lbnd(1)
    dy = this%ubnd(2) - this%lbnd(2)
    dz = this%ubnd(3) - this%lbnd(3)

    this%srfarea = 2.0_rp*(dx*dy + dx*dz + dy*dz)

    end subroutine

!*******************************************************************************

end module m_aabb
