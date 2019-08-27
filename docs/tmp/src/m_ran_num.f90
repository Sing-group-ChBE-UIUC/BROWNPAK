module m_ran_num
    !!  Provides random number generation procedures, mostly calling
    !! routine from Intel MKL VSL.

use m_precision
use mkl_vsl_type
use mkl_vsl

implicit none

private

public ::             & 
    init_stream,      &
    delete_stream,    &
    load_stream,      &
    save_stream,      &
    save_seed,        &
    get_iuniform,     &
    get_rv_iuniform,  &
    get_uniform,      &
    get_rv_uniform,   &
    get_rv_gaussian,  &
    ransphere

integer(ip), save :: seed
type(VSL_STREAM_STATE), save :: stream

contains

!*******************************************************************

subroutine init_stream(fn)
    !! Initializes a BRNG stream.

    !! The stream is initialized with a seed from the file `fn`. If
    !! `fn` is an empty string, the seed is obtained from `/dev/urandom`.

    character(len=*) :: fn
        !! Name of the file containing the RNG seed.
    integer(ip) :: fu
    integer(ip) :: stat

    if (len_trim(fn)==0) then
        open(newunit=fu, file='/dev/urandom', action='read', &
            form='unformatted', access='stream', status='old')
       read(fu) seed
       close(fu)
       seed = abs(seed)
    else 
       open(newunit=fu, file=fn, action='read', status='old')
       read(fu, *) seed
       close(fu)
    end if

    stat = vslNewStream(stream, VSL_BRNG_MT19937, seed)

    if (stat /= VSL_STATUS_OK) then
        write(*, *) 'error: vslNewStream'
    end if

    end subroutine

!*******************************************************************

subroutine delete_stream()
    !! Deletes a BRNG stream.

    integer(ip) :: stat

    stat = vslDeleteStream(stream)

    if (stat /= VSL_STATUS_OK) then
        write(*, *) 'error: vslDeleteStream'
    end if

    end subroutine

!*******************************************************************

subroutine load_stream(fn)
    !! Loads a BRNG stream from file `fn`.

    character(len=*), intent(in) :: fn
    integer(ip) :: stat

    stat = vslLoadStreamF(stream, fn)

    if (stat /= VSL_STATUS_OK) then
        write(*, *) 'error: vslLoadStream'
    end if

    end subroutine

!*******************************************************************

subroutine save_seed(fn)
    !! Saves the RNG seed to file `fn`.

    character(len=*), intent(in) :: fn
    integer(ip) :: fu

    open(newunit=fu, file=fn, action='write', status='replace')
    write(fu, *) seed
    close(fu)

    end subroutine

!*******************************************************************

subroutine save_stream(fn)
    !! Saves a BRNG stream to file `fn`.

    character(len=*), intent(in) :: fn
    integer(ip) :: stat

    stat = vslSaveStreamF(stream, fn)

    if (stat /= VSL_STATUS_OK) then
        write(*, *) 'error: vslSaveStream'
    end if

    end subroutine

!*******************************************************************

function get_uniform(lb, ub) result(res)
    !! Returns a random number from a uniform distribution.

    real(rp), intent(in) :: lb !! Lower bound
    real(rp), intent(in) :: ub !! Upper bound
    real(rp) :: res
    real(rp), dimension(1) :: rv
    integer(ip) :: stat

    stat = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, &
                        1, rv, lb, ub)
    res = rv(1)

    end function

!*******************************************************************

subroutine get_rv_uniform(lb, ub, rv, block_size)
    !! Returns a random vector from a uniform distribution.

    !! If a `block_size > 0` is provided, it fills `rv` in blocks of
    !! size `block_size`.

    real(rp), intent(in) :: lb !! Lower bound
    real(rp), intent(in) :: ub !! Upper bound
    real(rp), dimension(:), intent(out) :: rv
    integer(ip), intent(in), optional :: block_size
    integer(ip) :: stat
    integer(ip) :: num_blocks
    integer(ip) :: size_rv
    integer(ip) :: ini
    integer(ip) :: fin
    integer(ip) :: n
    
    if (.not. present(block_size)) then
        stat = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, &
                        size(rv), rv, lb, ub)
    else
        size_rv = size(rv)
        num_blocks = ceiling(size_rv/real(block_size))
        ini = 0
        fin = 0

        do
            if (fin >= size_rv) then
                exit
            end if

            if (size_rv <= block_size) then
                n = size_rv
            else if ((size_rv-fin) < block_size) then
                n = size_rv - fin
            else
                n = block_size
            end if

            ini = fin + 1
            fin = fin + n

            stat = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, &
                n, rv(ini:fin), lb, ub)
        end do
    end if

    end subroutine

!*******************************************************************

function get_iuniform(lb, ub) result(res)
    !! Returns a random integer from a uniform distribution.

    integer(ip), intent(in) :: lb
    integer(ip), intent(in) :: ub
    integer(ip) :: res
    integer(ip), dimension(1) :: rv
    integer(ip) :: stat

    stat = virnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, &
                        1, rv, lb, ub)
    res = rv(1)

    end function

!*******************************************************************

subroutine get_rv_iuniform(lb, ub, rv, block_size)
    !! Returns a random vector of integers from a uniform distribution.

    integer(ip), intent(in) :: lb
    integer(ip), intent(in) :: ub
    integer(ip), dimension(:), intent(out) :: rv
    integer(ip), intent(in), optional :: block_size
    integer(ip) :: stat
    integer(ip) :: num_blocks
    integer(ip) :: size_rv
    integer(ip) :: ini
    integer(ip) :: fin
    integer(ip) :: n
    
    if (.not. present(block_size)) then
        stat = virnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, &
                        size(rv), rv, lb, ub)
    else
        size_rv = size(rv)
        num_blocks = ceiling(size_rv/real(block_size))
        ini = 0
        fin = 0

        do
            if (fin >= size_rv) then
                exit
            end if

            if (size_rv <= block_size) then
                n = size_rv
            else if ((size_rv-fin) < block_size) then
                n = size_rv - fin
            else
                n = block_size
            end if

            ini = fin + 1
            fin = fin + n

            stat = virnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, &
                            n, rv(ini:fin), lb, ub)
        end do
    end if

    end subroutine

!*******************************************************************

subroutine get_rv_gaussian(mean, std_dev, rv, block_size)
    !! Generates a random vector of integers from a gaussian distribution.

    real(rp), intent(in) :: mean
    real(rp), intent(in) :: std_dev
    real(rp), dimension(:), intent(out) :: rv
    integer(ip), intent(in), optional :: block_size
    integer(ip) :: stat
    integer(ip) :: num_blocks
    integer(ip) :: size_rv
    integer(ip) :: ini
    integer(ip) :: fin
    integer(ip) :: n

    if (.not. present(block_size)) then
        stat = vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, &
                        size(rv), rv, mean, std_dev)
    else
        size_rv = size(rv)
        num_blocks = ceiling(size_rv/real(block_size))
        ini = 0
        fin = 0

        do
            if (fin >= size_rv) then
                exit
            end if

            if (size_rv <= block_size) then
                n = size_rv
            else if ((size_rv-fin) < block_size) then
                n = size_rv - fin
            else
                n = block_size
            end if

            ini = fin + 1
            fin = fin + n

            stat = vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, &
                        n, rv(ini:fin), mean, std_dev)
        end do
    end if

    end subroutine

!*******************************************************************

subroutine ransphere(r)
    !! Generates a random vector from the surface of a unit sphere.
    !!
    !! Algorithm from Allen & Tildesley (ed 1) p. 349.

    real(rp), dimension(3), intent(out) :: r(3)
    real(rp), dimension(2) :: zeta
    real(rp) :: zetasq
    real(rp) :: rt
    
    r = 0.0_rp
    zetasq = 2.0_rp ! Any value greater than 1
    
    do 
        call get_rv_uniform(-1.0_rp, 1.0_rp, zeta)
        zetasq = zeta(1)*zeta(1) + zeta(2)*zeta(2)
        rt = sqrt(1.0_rp - zetasq)
        r(1) = 2.0_rp*zeta(1)*rt
        r(2) = 2.0_rp*zeta(2)*rt
        r(3) = 1.0_rp - 2.0_rp*zetasq
    
        if (zetasq <= 1.0_rp) exit
    end do

    end subroutine

!*******************************************************************

end module m_ran_num
