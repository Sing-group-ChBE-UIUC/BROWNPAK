program main

use m_precision
use m_strings
use m_logger, only: logger_init, logger => master_logger
use m_globals, only: sim_style, job_tag
use m_control_io, only: read_control
use m_setup, only: setup, finish
use m_relax, only: rlx_run
use m_bd_solver, only: bds_run
use m_mpcd, only: mpcd_run

implicit none

!*******************************************************************************

character(len=64) :: cla
character(len=:), allocatable :: key
character(len=:), allocatable :: val
character(len=:), allocatable :: fn_control
character(len=512) :: msg
real(rp):: t_start, t_end
integer :: ierr
integer :: icla
integer :: ncla ! number of command line arguments, without the command name

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

!Initialize logger
call logger_init('brownpak.log'//trim(adjustl(job_tag)), stderr_threshold=100, &
        stdout_threshold=100, logfile_threshold=0)

call read_control(fn_control)

call setup()

call cpu_time(t_start)

select case(sim_style)
case(0)
    call rlx_run()
case(1)
    call bds_run()
case(2)
    call mpcd_run()
case default
    continue
end select

call cpu_time(t_end)

write(msg,'(g0.6)') (t_end-t_start)
call logger%info('Total execution time (s)', trim(adjustl(msg)))

call finish()

call logger%destroy()

!*******************************************************************************

end program
