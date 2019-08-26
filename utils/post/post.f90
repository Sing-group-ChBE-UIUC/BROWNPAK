program post
!TODO: Postprocess depending on stuff to be calculated

use m_precision
use m_constants_math
use m_strings
use m_trajectory
use m_globals
use m_config_io

implicit none

character(len=256) :: cla_buf
character(len=:), allocatable :: frm_dir
character(len=:), allocatable :: fn_ldf
integer :: ifrm_beg, ifrm_end, ifrm_stp
integer :: i, iframe, oframe, ierr, nofrms

call get_command_argument(1, cla_buf)
fn_cfg = trim(adjustl(cla_buf))

call get_command_argument(2, cla_buf)
fn_traj = trim(adjustl(cla_buf))

call get_command_argument(3, cla_buf)
frm_dir = trim(adjustl(cla_buf))

call get_command_argument(4, cla_buf)
ifrm_beg = str_to_i(cla_buf)

call get_command_argument(5, cla_buf)
ifrm_end = str_to_i(cla_buf)

call get_command_argument(6, cla_buf)
ifrm_stp = str_to_i(cla_buf)

!Read cfg file
call read_config(fn_cfg)

!Open file for reading trajectory
call traj%open(fn_traj, 'r')
write(*,*) 'numframes ', traj%num_frames

if (ifrm_end == -1) ifrm_end = traj%num_frames

!Number of output frames
nofrms = (ifrm_end - (ifrm_beg-1))/ifrm_stp
if (nofrms > 1000) then
    write(*,*) 'More than 1000 frames will be written'
    write(*,'(a)', advance='no') 'Press "y" if ok: '
    read(*,*) cla_buf
    !Will fall through if not 'y' 
    if ( trim(adjustl(cla_buf)) /= 'y' ) ifrm_end = ifrm_beg - 1
end if

do iframe = ifrm_beg, ifrm_end, ifrm_stp
    oframe = (iframe - ifrm_beg + 1)/ifrm_stp
    write(*,*) oframe, iframe
    call traj%read(iframe, nts, ierr, coordinates)
    fn_ldf = frm_dir//'/frame-'//str_from_num(oframe)//'.txt'
    call write_ldf(fn_ldf, 'Frame-'//str_from_num(iframe))
end do

!Close trajectory
call traj%close()

!*******************************************************************************

end program
