program main

use m_precision
use m_strings
use m_logger
use m_globals, only: job_tag
use m_control_io, only: read_control
use m_setup, only: setup, finish, run

implicit none

!*******************************************************************************

character(len=64) :: cla
character(len=:), allocatable :: key
character(len=:), allocatable :: val
character(len=:), allocatable :: fn_control
character(len=128) :: msg
integer :: ierr
integer :: icla
integer :: ncla ! number of command line arguments, without the command name
integer(int64) :: ccbeg, ccend, crate

!Two command line arguments may be provided -- (i) fn_control=val and 
!(ii) job_tag=val

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

call read_control(fn_control)

call logger_init(logger, 'brownpak.log'//job_tag, .true.)

call setup()

call system_clock(ccbeg, crate)

call run()

call system_clock(ccend, crate)

write(msg,'(a,es12.5)') 'execution time(s) = ', (ccend-ccbeg)/real(crate,rp)

call logger%log_msg(msg)

call logger%finish()

call finish()

!*******************************************************************************

end program
