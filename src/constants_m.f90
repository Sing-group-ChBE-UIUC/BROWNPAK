module constants_m
    !! Various math constants.

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

real (rp), parameter :: math_third = 0.333333333333333_rp

!>pi
real (rp), parameter :: math_pi = 3.1415926535897931_rp

!>pi divided by two
real (rp), parameter :: math_pi_2 = 1.5707963267948966_rp

!>pi divided by four
real (rp), parameter :: math_pi_4 = 0.78539816339744828_rp

!>reciprocal of pi
real (rp), parameter :: math_1_pi = 0.31830988618379069_rp

!>two times reciprocal of pi
real (rp), parameter :: math_2_pi = 0.63661977236758138_rp

!>two times the reciprocal of the square root of pi. 
real (rp), parameter :: math_2_sqrtpi = 1.1283791670955126_rp

!>square root of two
real (rp), parameter :: math_sqrt2 = 1.4142135623730951_rp

!>cube root of two
real (rp), parameter :: math_cbrt2 = 1.2599210498948732_rp

!>sixth root of two
real (rp), parameter :: math_sxrt2 = 1.122462048309373_rp

!>reciprocal of the square root of two
real (rp), parameter :: math_sqrt1_2 = 0.70710678118654746_rp

!>square root of three
real (rp), parameter :: math_sqrt3 = 1.7320508075688772_rp

!>square root of M_E
real (rp), parameter :: math_sqrt_e = 1.6487212707001282_rp

!>square root of pi
real (rp), parameter :: math_sqrt_pi = 1.7724538509055159_rp

!>The base of natural logarithms
real (rp), parameter :: math_e = 2.7182818284590451_rp

!>The logarithm of M_E to base two
real (rp), parameter :: math_log2e = 1.4426950408889634_rp

!>The logarithm of M_E to base 10
real (rp), parameter :: math_log10e = 0.43429448190325182_rp

!>The natural logarithm of two
real (rp), parameter :: math_ln2 = 0.69314718055994529_rp

!>The natural logarithm of 10
real (rp), parameter :: math_ln10 = 2.3025850929940459_rp

end module constants_m
