module m_precision

use, intrinsic ::  iso_fortran_env, only: int32, int64, real64

implicit none

integer, parameter :: ip = int32
    !! Default integer precision
integer, parameter :: ip_long = int64
    !! Default long integer precision

integer, parameter :: rp = real64
    !! Default real precision

integer, parameter :: sizeof_char = 1
    !! Size of a char in bytes
integer, parameter :: sizeof_int  = 4
    !! Size of a default int in bytes
integer, parameter :: sizeof_long_int  = 8
    !! Size of a default long int in bytes
integer, parameter :: sizeof_real = 8
    !! Size of a default real in bytes

end module m_precision
