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

module utils_math_m
    !!Various (mostly linear algebra) functions, particularly for use with small
    !!matrices.

use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
use constants_m

implicit none

interface allclose
    !! Checks if two arrays are elementwise close within tolerance
    module procedure allclose_rank1
    module procedure allclose_rank2
    module procedure allclose_rank3
end interface 

interface swap
    !! Swaps two arrays
    module procedure swap_integer
    module procedure swap_real
    module procedure swap_complex
end interface swap

contains

!************************************************************************

elemental subroutine rad2deg(rad, deg)

    real(rp), intent(in) :: rad
    real(rp), intent(out) :: deg

    deg = (180.0_rp/math_pi)*rad

    end subroutine

!************************************************************************

elemental subroutine deg2rad(deg, rad)

    real(rp), intent(in) :: deg
    real(rp), intent(out) :: rad

    rad = (math_pi/180.0_rp)*deg

    end subroutine

!************************************************************************

subroutine cross(a, b, c)
    !!  Calculates the cross product between two 3-element vectors
    real(rp), dimension(3), intent(in)  :: a
    real(rp), dimension(3), intent(in)  :: b
    real(rp), dimension(3), intent(out) :: c
        !! Cross product of `a` and `b`; **c** = **a** x **b**
    
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

    end subroutine

!******************************************************************************

subroutine cross_mat(a, mat)
    !! Calculates the cross product matrix of a 3-element vector. The cross
    !! product matrix **A** of **a** is defined as **a** x **b** = **A** . **b**,
    !! where **b** is another 3-element vector.
    real(rp), dimension(3), intent(in)    :: a
    real(rp), dimension(3,3), intent(out) :: mat
        !! Cross product matrix of `a`
    
    mat = 0.0_rp
    mat(2,1) =  a(3)
    mat(3,1) = -a(2)
    mat(1,2) = -a(3)
    mat(3,2) =  a(1)
    mat(1,3) =  a(2)
    mat(2,3) = -a(1)

    end subroutine

!******************************************************************************

subroutine outer(a, b, c)
    !! Calculates the outer product of two vectors, \(c_{ij} = a_i  b_j\).
    real(rp), dimension(:), intent(in)    :: a
        !! (m,) array
    real(rp), dimension(:), intent(in)    :: b
        !! (n,) array
    real(rp), dimension(:,:), intent(out) :: c
        !! (m,n) array; Outer product
    integer :: m, n
    integer :: i, j

    m = size(a)
    n = size(b)

    do j = 1, n
        do i = 1, m
            c(i,j) = a(i)*b(j)
        end do
    end do

    end subroutine

!******************************************************************************

function scalar_triple_product (a, b, c) result(res)
    !! Returns the scalar triple product **a**.(**b** x **c**)
    real(rp), dimension(3), intent(in)  :: a
    real(rp), dimension(3), intent(in)  :: b
    real(rp), dimension(3), intent(in)  :: c
    real(rp)  :: res

    res = a(1)*( b(2)*c(3)-c(2)*b(3) ) - a(2)*( b(1)*c(3)-c(1)*b(3) ) &
            + a(3)*( b(1)*c(2)-c(1)*b(2) )

    end function

!******************************************************************************

subroutine vector_triple_product (a, b, c, d)
    !! Returns the vector triple product **d** = **a** x (**b** x **c**)
    real(rp), dimension(3), intent(in)  :: a
    real(rp), dimension(3), intent(in)  :: b
    real(rp), dimension(3), intent(in)  :: c
    real(rp), dimension(3), intent(out) :: d
     !! Vector triple product

    d = b*( a(1)*c(1) + a(2)*c(2) + a(3)*c(3) )  &
        - c*( a(1)*b(1) + a(2)*b(2) + a(3)*b(3) )

    end subroutine

!******************************************************************************

function det(A) result(res)
    !!  Returns the determinant of an (N x N) matrix, where N = 2, 3, or 4.
    !!
    !!  Original routine by [David Simpson](http://www.davidgsimpson.com/software.html)
    !!  @note For a general NxN matrix do an LU decomp
    !""
    real(rp), dimension(:,:), intent(in) :: A
        !! (N,N) array, where N = 2, 3, or 4.
    real(rp) :: res
    integer :: nrows

    nrows = size(A, 1)

    if (nrows==2) then
        res = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    else if (nrows==3) then
        res = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    else if (nrows==4) then
        res = A(1,1) * ( A(2,2) * (A(3,3)*A(4,4) - A(3,4)*A(4,3))    &
                       + A(2,3) * (A(3,4)*A(4,2) - A(3,2)*A(4,4))    &
                       + A(2,4) * (A(3,2)*A(4,3) - A(3,3)*A(4,2)) )  &
            - A(1,2) * ( A(2,1) * (A(3,3)*A(4,4) - A(3,4)*A(4,3))    &
                       + A(2,3) * (A(3,4)*A(4,1) - A(3,1)*A(4,4))    &       
                       + A(2,4) * (A(3,1)*A(4,3) - A(3,3)*A(4,1)) )  &
            + A(1,3) * ( A(2,1) * (A(3,2)*A(4,4) - A(3,4)*A(4,2))    &
                       + A(2,2) * (A(3,4)*A(4,1) - A(3,1)*A(4,4))    &
                       + A(2,4) * (A(3,1)*A(4,2) - A(3,2)*A(4,1)) )  &
            - A(1,4) * ( A(2,1) * (A(3,2)*A(4,3) - A(3,3)*A(4,2))    &       
                       + A(2,2) * (A(3,3)*A(4,1) - A(3,1)*A(4,3))    &
                       + A(2,3) * (A(3,1)*A(4,2) - A(3,2)*A(4,1)) )

     else
         stop 'matrix size must be 2, 3, or 4'
     end if

    end function

!******************************************************************************

function trace(mat) result (res)
    !! Returns the trace of a square matrix

    real(rp), intent(in) :: mat(:,:)
        !! (N,N) array
    real(rp) :: res
    integer :: nrows
    integer :: i

    nrows = size(mat, 1)
    res = 0.0_rp

    do i = 1, nrows
        res = res + mat(i,i)
    end do

    end function

!******************************************************************************
 
elemental logical function isclose(a, b, rel_tol, abs_tol)
    !!  Checks if two floating point numbers of type double are close within
    !!  tolerance.
    !!
    !!  Based on python implementation at
    !!  <https://github.com/PythonCHB/close_pep/blob/master/is_close.py>.
    !!  The *method='weak'* option is used here.
    !""

    real(rp), intent(in)    :: a
    real(rp), intent(in)    :: b
    real(rp), intent(in), optional    :: rel_tol
        !! Relative tolerance, `rel_tol` >= 0, default 1e-10
    real(rp), intent(in), optional    :: abs_tol
        !! Absolute tolerance, `abs_tol` >= 0, default 0.0
    real(rp) :: rel_tol_
    real(rp) :: abs_tol_
    real(rp) :: diff

    rel_tol_ = 1e-10_rp
    abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol
    
    if (a == b) then  ! short-circuit exact equality
        isclose = .true.
    end if
    
    if ((.not. ieee_is_finite(a)) .or. (.not. ieee_is_finite(b))) then
         ! Includes the case of two infinities of opposite sign, or
         ! one infinity and one finite number. Two infinities of opposite sign
         ! would otherwise have an infinite relative tolerance.
         isclose = .false.
    end if
    
    diff = abs(b - a)
    isclose = ( ((diff <= abs(rel_tol_*b)) .or. (diff <= abs(rel_tol_*a))) &
            .or. (diff <= abs_tol_) )

    end function

!******************************************************************************

logical function allclose_rank1(a, b, rel_tol, abs_tol)
    !!  Checks if two rank-1 floating point arrays of type double are close within
    !!  tolerance.
    !""
    real(rp), dimension(:), intent(in) :: a
        !! (m,) array
    real(rp), dimension(:), intent(in) :: b
        !! (m,) array
    real (rp), intent(in), optional :: rel_tol
        !! Relative tolerance; default 1e-10
    real (rp), intent(in), optional :: abs_tol
        !! Absolute tolerance; default 0.0
    real (rp) :: rel_tol_ = 1e-10_rp
    real (rp) :: abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol

    allclose_rank1 = all(isclose(a, b, rel_tol_, abs_tol_))

    end function


logical function allclose_rank2(a, b, rel_tol, abs_tol)
    !!  Checks if two rank-2 floating point arrays of type double are close within
    !!  tolerance.
    !""
    real(rp), dimension(:,:), intent(in) :: a
        !! (m,n) array
    real(rp), dimension(:,:), intent(in) :: b
        !! (m,n) array
    real (rp), intent(in), optional :: rel_tol
        !! Relative tolerance; default 1e-10
    real (rp), intent(in), optional :: abs_tol
        !! Absolute tolerance; default 0.0
    real (rp) :: rel_tol_ = 1e-10_rp
    real (rp) :: abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol

    allclose_rank2 = all(isclose(a, b, rel_tol_, abs_tol_))

    end function


logical function allclose_rank3(a, b, rel_tol, abs_tol)
    !!  Checks if two rank-3 floating point arrays of type double are close within
    !!  tolerance.
    !""
    real(rp), dimension(:,:,:), intent(in) :: a
        !! (m,n,p) array
    real(rp), dimension(:,:,:), intent(in) :: b
        !! (m,n,p) array
    real (rp), intent(in), optional :: rel_tol
        !! Relative tolerance; default 1e-10
    real (rp), intent(in), optional :: abs_tol
        !! Absolute tolerance; default 0.0
    real (rp) :: rel_tol_ = 1e-10_rp
    real (rp) :: abs_tol_ = 0.0_rp

    if (present(rel_tol)) rel_tol_ = rel_tol
    if (present(abs_tol)) abs_tol_ = abs_tol

    allclose_rank3 = all(isclose(a, b, rel_tol_, abs_tol_))

    end function

!******************************************************************************

elemental subroutine swap_integer(a, b)
    integer, intent(in out) :: a
    integer, intent(in out) :: b
    integer :: temp

    temp = a
    a = b
    b = temp
    end subroutine

elemental subroutine swap_real(a, b)
    real(rp), intent(in out) :: a
    real(rp), intent(in out) :: b
    real(rp) :: temp

    temp = a
    a = b
    b = temp
    end subroutine

elemental subroutine swap_complex(a, b)
    complex(rp), intent(in out) :: a
    complex(rp), intent(in out) :: b
    complex(rp) :: temp

    temp = a
    a = b
    b = temp
    end subroutine

!******************************************************************************

subroutine unitize(a)
    !! Normalizes a vector in-place. If the magnitude of the vector is nearly
    !! zero, no normalization takes place and the vector is returned as is with
    !! a warning message.
    real(rp), dimension(:), intent(in out)  :: a
        !! (m,) array
    real(rp) :: norm

    norm = norm2(a)

    if (isclose(norm, 0.0_rp, rel_tol=1.e-15_rp, abs_tol=0.0_rp)) then
        write(*,*) '[unitize] norm close to zero'
    else
        a = a/norm2(a)
    end if

    end subroutine

!******************************************************************************

subroutine linspace(start, finish, num, val, step)
    !! Generates evenly spaced numbers over a specified interval. Both end
    !! points are included. If `start` < `finish`, the returned step size (if
    !! `step` is present) will be negative.
    !""
    real(rp), intent(in) :: start
        !! Starting point
    real(rp), intent(in) :: finish
        !! Ending point, `finish` /= `start`
    integer, intent(in) :: num
        !! Number of values to generate, `num >= 2`
    real(rp), dimension(:), intent(out) :: val
        !! (`num`,) array; Generated values
    real(rp), intent(out), optional :: step
        !! Step size
    real(rp) :: step_
    integer :: i

    if (num < 2) then
        write(*,*) '[linspace] `num` must be >= 2'
        stop
    end if

    if (num < size(val)) then
        write(*,*) '[linspace] `size(val)` must be >= `num`'
        stop
    end if

    if (start == finish) then
        write(*,*) '[linspace] `start` must not be equal to `finish`'
        stop
    end if

    val = 0.0_rp

    step_ = (finish-start)/(num-1)

    if (present(step)) step = step_

    val(1) = start
    val(num) = finish

    do i = 2, (num-1)
        val(i) = start + (i-1)*step_
    end do

    end subroutine

!******************************************************************************

subroutine logspace(start, finish, num, val, base)
    !!  Generates numbers spaced evenly on a log scale.
    !!  In linear space, the sequence starts at `base ** start`
    !!  (`base` to the power of `start`) and ends with `base ** stop`
    !""
    real(rp), intent(in) :: start
        !! Starting point, `base ** start` is the starting value
    real(rp), intent(in) :: finish
        !! Ending point, `finish` /= `start`,  `base ** start`
        !! is the ending value
    integer, intent(in) :: num
        !! Number of values to generate, `num >= 2`
    real(rp), dimension(:), intent(out) :: val
        !! (`num`,) array; Generated values
    real(rp), intent(in), optional :: base
        !! Base of the logspace, default 10
    real(rp) :: base_
    integer :: i

    if (num < 2) then
        write(*,*) 'logspace `num` must be >= 2'
        stop
    end if

    if (num < size(val)) then
        write(*,*) 'logspace `size(val)` must be >= `num`'
        stop
    end if

    if (start == finish) then
        write(*,*) 'logspace `start` must not be equal to `finish`'
        stop
    end if

    val = 0.0_rp
    base_ = 10.0_rp
    if (present(base)) base_ = base

    call linspace(start, finish, num, val)

    do i = 1, num
        val(i) = base_**val(i)
    end do

    end subroutine

!******************************************************************************

subroutine identity(mat_eye)
    !! Creates an identity matrix of size n x n.
    real(rp), dimension(:,:), intent(out) :: mat_eye
        !! (n,n) array
    integer :: nrows
    integer :: i

    nrows = size(mat_eye, 1)

    mat_eye = 0.0_rp

    do i = 1, nrows
        mat_eye(i,i) = 1.0_rp
    end do
    
    end subroutine

!******************************************************************************

subroutine get_diagonal(mat, d)
    !! Returns the diagonal elements of a square matrix.
    real(rp), dimension(:,:), intent(in) :: mat
        !! (n,n) array
    real(rp), dimension(:), intent(out) :: d
        !! (n,) array; contains the entries of the main diagonal
    integer :: n
    integer :: i
    
    n = size(mat, 1)

    do i = 1, n
        d(i) = mat(i,i)
    end do

    end subroutine

!******************************************************************************

subroutine add_transpose(mat)
    !! Adds a square matrix and its transpose in place: \(A_{ij} = A_{ij } + A_{ji}\)
    real(rp), dimension(:,:), intent(in out) :: mat
        !! (n,n) array
    integer :: n
    integer :: i
    integer :: j

    n = size(mat,1)

    !Dealing with the upper triangular part (including diagonal)
    do j = 1, n
        do i = 1, j
            mat(i,j) = mat(i,j) + mat(j,i)
        end do
    end do

    !Dealing with the strictly lower triangular part (i.e. excluding diagonal)
    do j = 1, (n-1)
        do i = (j+1), n
            mat(i,j) = mat(j,i)
        end do
    end do

    end subroutine

!******************************************************************************

subroutine subtract_transpose(mat)
    !! Calculates the difference of a square matrix and its transpose in place:
    !! \(A_{ij} = A_{ij } - A_{ji}\)
    real(rp), dimension(:,:), intent(in out) :: mat
    integer :: n
    integer :: i
    integer :: j

    n = size(mat,1)

    !Dealing with the strictly upper triangular part (i.e. excluding diagonal)
    do j = 2, n
        do i = 1, (j-1)
            mat(i,j) = mat(i,j) - mat(j,i)
        end do
    end do

    !Dealing with the strictly lower triangular part (i.e. excluding diagonal)
    do j = 1, (n-1)
        do i = (j+1), n
            mat(i,j) = -mat(j,i)
        end do
    end do

    !Putting the diagonal to zero
    do i = 1, n
        mat(i,j) = 0.0_rp
    end do

    end subroutine

!******************************************************************************

subroutine multiply_transpose(A, B)
    !! Multiplies a matrix with its transpose:
    !! \(\mathbf{\mathrm{B}} = \mathbf{\mathrm{A}} \cdot \mathbf{\mathrm{A}}^T\)
    !""
    real(rp), dimension(:,:), intent(in) :: A
        !! (m,n) array
    real(rp), dimension(:,:), intent(out) :: B
        !! (m,m) array
        !A is m x n, A^T is n x m, so B is m x m

    B = matmul(A, transpose(A))

    end subroutine

!******************************************************************************

function get_quad_form(A, x) result(res)
    !!  Calculates the quadratic form *x^T A x*, where A is an *n x n* matrix and *x* is a
    !!  vector of length *n*
    !""
    real(rp), dimension(:,:), intent(in) :: A
        !! (n,n) array
    real(rp), dimension(:), intent(in) :: x
        !! (n,) array
    real(rp) :: res
    integer :: n
    integer :: i, j

    n = size(A,2)
    res = 0.0_rp

    do j = 1, n
        do i = 1, n
            res = res + A(i,j)*x(j)*x(i)
        end do
    end do

    end function

!******************************************************************************

subroutine orth(a)
    !! Orthogonalizes a set of vectors in-place using Gram-Schmidt orthonormalization
    !! 
    !! *Reference:* Golub and Van Loan, Matrix Computations, 3rd edition, Section 5.2.8,
    !! Algorithm 5.2.5, p. 231.

    real(rp), dimension(:,:), intent(in out) :: a
        !! (m,n) array, where m <= n. The first m columns of the matrix are
        !!  overwritten with the orthogonal basis vectors.
    integer :: m
    integer :: i, j

    m = size(a, 1)

    do i = 1, m
        a(i,:) = a(i,:)/norm2(a(i,:)) 
        do j = (i+1), m
            a(j,:) = a(j,:) - dot_product(a(j,:),a(i,:))*a(i,:)
        end do
    end do

    end subroutine

!******************************************************************************

subroutine invert_mat33(a, inv_a)
    !! Inverts a 3x3 matrix.
    !!
    !!*Reference:* https://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf

    real(rp), dimension(3,3),  intent(in) :: a
    real(rp), dimension(3,3), intent(out) :: inv_a
    real(rp) :: det_a

    inv_a(1,1) =   a(2,2)*a(3,3) - a(2,3)*a(3,2)
    inv_a(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))
    inv_a(3,1) =   a(2,1)*a(3,2) - a(2,2)*a(3,1)

    inv_a(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))
    inv_a(2,2) =   a(1,1)*a(3,3) - a(1,3)*a(3,1)
    inv_a(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))

    inv_a(1,3) =   a(1,2)*a(2,3) - a(1,3)*a(2,2)
    inv_a(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))
    inv_a(3,3) =   a(1,1)*a(2,2) - a(1,2)*a(2,1)

    !Determinant
    det_a = a(1,1)*inv_a(1,1) + a(1,2)*inv_a(2,1) + a(1,3)*inv_a(3,1)
    inv_a = inv_a/det_a
    
    end subroutine

!******************************************************************************

subroutine eigval_33rsym(a, ev)
    !! Calculates the eigenvalues of a 3 x 3 real symmetric matrix. The
    !! eigenvalues calculated are in decreasing order. Only the diagonal and
    !! lower triangular part of the matrix is accessed.
    !!
    !! *Reference:* https://en.wikipedia.org/wiki/Eigenvalue_algorithm#cite_note-Smith-12
    !!
    !! See also David Eberly's notes and implementation at 
    !! https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf

    real(rp), dimension(3,3), intent(in) :: a
    real(rp), dimension(3), intent(out) :: ev
    real(rp), dimension(3,3) :: eye
    real(rp), dimension(3,3) :: b
    real(rp) :: p, q, r
    real(rp) :: p1, p2
    real(rp) :: phi

    !Sum of the elements in the lower triangular part
    p1 = a(2,1)**2 + a(3,1)**2 + a(3,2)**2 

    if ( isclose(p1, 0.0_rp) ) then
        !Matrix a is diagonal
        ev(1) = a(1,1); ev(2) = a(2,2); ev(3) = a(3,3)
        return
    else
        q = (a(1,1) + a(2,2) + a(3,3))/3.0_rp
        p2 = (a(1,1) - q)**2 + (a(2,2) - q)**2 + (a(3,3) - q)**2 + 2*p1
        p = sqrt(p2/6.0_rp)

        call identity(eye)

        b = (a - q*eye)/p
        r = 0.5_rp*det(b)

        !r must be within [-1, 1] in exact computation; but need to handle
        !slightly out of range in computation.
        if (r <= -1.0_rp) then
            phi = math_pi/3.0_rp
        else if ( r >= 1.0_rp) then
            phi = 0.0_rp
        else
            phi = acos(r)/3.0_rp
        end if

        !Eigen values are ordered as ev(3) <= ev(2) <= ev(1) 
        ev(1) = q + 2*p*cos(phi)
        ev(3) = q + 2*p*cos(phi + (2*math_pi/3.0_rp))
        ev(2) = 3*q - ev(1) - ev(3)
    end if

    end subroutine

!******************************************************************************

subroutine dsyevc3(a, w)
    !!author: Joachim Kopp
    !!date: 2006
    !!
    !! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
    !! analytical algorithm.
    !! Only the diagonal and upper triangular parts of A are accessed. The access
    !! is read-only.
    !!
    !! Copyright (C) 2006  Joachim Kopp
    !  https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    ! .. Arguments ..
      REAL(RP), DIMENSION(3,3), INTENT(IN) :: A
      REAL(RP), DIMENSION(3), INTENT(OUT) ::  W(3)

     !.. Local Variables ..
      REAL(RP) ::  M, C1, C0
      REAL(RP) ::  DE, DD, EE, FF
      REAL(RP) ::  P, SQRTP, Q, C, S, PHI
  
     !Determine coefficients of characteristic poynomial. We write
     !      | A   D   F  |
     ! A =  | D*  B   E  |
     !      | F*  E*  C  |

      DE    = A(1,2) * A(2,3)
      DD    = A(1,2)**2
      EE    = A(2,3)**2
      FF    = A(1,3)**2
      M     = A(1,1) + A(2,2) + A(3,3)
      C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) &
               - (DD + EE + FF)
      C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) &
               - 2.0_RP * A(1,3)*DE

      P     = M**2 - 3.0_RP * C1
      Q     = M*(P - (3.0_RP/2.0_RP)*C1) - (27.0_RP/2.0_RP)*C0
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0_RP * ( 0.25_RP * C1**2 * (P - C1) &
                + C0 * (Q + (27.0_RP/4.0_RP)*C0) )
      PHI   = (1.0_RP/3.0_RP) * ATAN2(SQRT(ABS(PHI)), Q)

      C     = SQRTP * COS(PHI)
      S     = (1.0_RP/MATH_SQRT3) * SQRTP * SIN(PHI)

      W(2) = (1.0_RP/3.0_RP) * (M - C)
      W(3) = W(2) + S
      W(1) = W(2) + C
      W(2) = W(2) - S

      END SUBROUTINE

!******************************************************************************

end module utils_math_m
