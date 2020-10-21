!********************************************************************************!
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2020 Sarit Dutta <saridut@gmail.com>                             !
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
!********************************************************************************!

program main

use constants_m
use strings_m
use logger_m
use control_m
use setup_m

implicit none

!*******************************************************************************

character(len=64) :: cla
character(len=:), allocatable :: key
character(len=:), allocatable :: val
character(len=:), allocatable :: job_tag
character(len=:), allocatable :: fn_control
type(ctrlpar_t) :: cpar
integer :: ierr
integer :: icla
integer :: ncla ! number of command line arguments, without the command name

!Two command line arguments may be provided -- (i) fn_control=val and 
!(ii) job_tag=val
!
!Usage: ./brownpak fn_control=<str> job_tag=<int>

fn_control = 'control.txt'

ncla = command_argument_count()
do icla = 1, ncla
    call get_command_argument(icla, value=cla, status=ierr)
    if (ierr > 0) then
        write(*, '(a,1x,i0)') "read failure for command argument", icla
        stop
    else if (ierr == -1) then
        write(*, '(a,1x,i0)') "command argument truncated", icla
        stop
    end if

    call str_get_keyval(cla, key, val, '=')
    if (key == 'fn_control') then
        fn_control = val
    else if (key == 'job_tag') then
        job_tag = '.'//val
    end if
end do

call cpar%read(fn_control)

call logger%init('brownpak.log'//job_tag, .true.)

call run(cpar, job_tag)

call logger%finish()

!*******************************************************************************

end program
