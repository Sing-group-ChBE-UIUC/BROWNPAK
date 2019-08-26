module m_simbox
!! Implements a simulation box with appropriate boundary conditions.

use m_precision
use m_constants_math
use m_ran_num

implicit none

type smbx_t
    real(rp), dimension(3,3) :: basis
    real(rp), dimension(3,3) :: dl_basis
    real(rp) :: volume
    logical :: is_deforming
    logical :: is_aligned
        !! Whether the basis vectors are aligned with the laboratory frame

    contains
        procedure :: set_basis => smbx_set_basis
        procedure :: freeze => smbx_freeze
        procedure :: unfreeze => smbx_unfreeze
        procedure :: get_image => smbx_get_image
        procedure :: wrap_all => smbx_wrap_all
        procedure :: to_center => smbx_to_center
        procedure :: get_rnd_points => smbx_get_rnd_points
end type smbx_t


contains

!********************************************************************************

subroutine smbx_init(this)
    !!Creates an instance of *smbx_t*. Can also be called to reset.

    type(smbx_t), intent(in out) :: this

    !Initialize to an identity matrix
    this%basis = 0.0_rp
    this%basis(1,1) = 1.0_rp
    this%basis(2,2) = 1.0_rp
    this%basis(3,3) = 1.0_rp
    this%volume = 1.0_rp

    this%dl_basis = this%basis
    this%is_deforming = .false.
    this%is_aligned = .true.

    end subroutine

!********************************************************************************

subroutine smbx_set_basis(this, bv)
    !!Sets all three basis vectors.

    class(smbx_t), intent(in out) :: this
    real(rp), dimension(3,3), intent(in) :: bv
    real(rp), dimension(3) :: a, b, c

    this%basis = bv

    a = this%basis(:,1); b = this%basis(:,2); c = this%basis(:,3)

    this%volume = a(1)*( b(2)*c(3)-c(2)*b(3) ) - a(2)*( b(1)*c(3)-c(1)*b(3) ) &
            + a(3)*( b(1)*c(2)-c(1)*b(2) )

    end subroutine

!********************************************************************************

subroutine smbx_freeze(this)
    !!Specifies *this* as non-deforming.

    class(smbx_t), intent(in out) :: this

    this%is_deforming = .false.

    end subroutine

!********************************************************************************

subroutine smbx_unfreeze(this)
    !!Specifies *this* as deforming.

    class(smbx_t), intent(in out) :: this

    this%is_deforming = .true.

    end subroutine

!********************************************************************************

subroutine smbx_get_image(this, r)
    !!Returns the image of *r* under PBC.
    !https://scicomp.stackexchange.com/questions/20165/periodic-boundary-conditions-for-triclinic-box

    class(smbx_t), intent(in) :: this
    real(rp), dimension(3), intent(in out) :: r
    real(rp), dimension(3) :: rf
    
    if (this%is_aligned) then
        r(1) = r(1) - this%basis(1,1)*nint( r(1)/this%basis(1,1) )
        r(2) = r(2) - this%basis(2,2)*nint( r(2)/this%basis(2,2) )
        r(3) = r(3) - this%basis(3,3)*nint( r(3)/this%basis(3,3) )
    else
        rf = matmul(this%dl_basis, r)
        rf = rf - nint(rf)
        r = matmul(this%basis, rf)
    end if

    end subroutine

!********************************************************************************

subroutine smbx_wrap_all(this, coords)
    !!Wraps atom positions w.r.t. periodic boundary conditions.
    !https://scicomp.stackexchange.com/questions/20165/periodic-boundary-conditions-for-triclinic-box

    class(smbx_t), intent(in) :: this
    real(rp), dimension(:,:), intent(in out) :: coords
    real(rp), dimension(3) :: rf
    real(rp), dimension(3) :: diag
    integer :: n, i
    
    n = size(coords,2)

    if (this%is_aligned) then
        diag = [this%basis(1,1), this%basis(2,2), this%basis(3,3)]
        do i = 1, n
            coords(:,i) = coords(:,i) - diag*floor( coords(:,i)/diag )
        end do
    else
        do i = 1, n
            rf = matmul(this%dl_basis, coords(:,i))
            rf = rf - floor(rf)
            coords(:,i) = matmul(this%basis, rf)
        end do
    end if

    end subroutine

!********************************************************************************

subroutine smbx_to_center(this, coords, com)
    !!Adjusts atom positions such that the c.o.m. of the atoms is at the center 
    !! of the box. Assumes all atoms to have the same mass and aligned axis.
    !! Optionally returns the original c.o.m.

    class(smbx_t), intent(in) :: this
    real(rp), dimension(:,:), intent(in out) :: coords
    real(rp), dimension(3), intent(out), optional :: com
    real(rp), dimension(3) :: half_diag, com_
    integer :: n, i
    
    half_diag = 0.5_rp*[this%basis(1,1), this%basis(2,2), this%basis(3,3)]
    n = size(coords,2)
    com_ = sum(coords, 2)/n

    do i = 1, n
        coords(:,i) = coords(:,i) - com_ + half_diag
    end do

    if (present(com)) com = com_ - half_diag

    end subroutine

!********************************************************************************

subroutine smbx_get_rnd_points(this, coords)
    !!Returns uniformly distributed points within the box.

    class(smbx_t), intent(in) :: this
    real(rp), dimension(:,:), intent(out) :: coords
    real(rp), dimension(3) :: r
    integer :: n, i
    
    n = size(coords,2)
    do i = 1, n
        call get_rv_uniform(0.0_rp, 1.0_rp, r)
        coords(:,i) = matmul(this%basis, r)
    end do

    end subroutine

!********************************************************************************

end module m_simbox
