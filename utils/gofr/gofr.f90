program main

use m_precision
use m_constants_math
use m_trajectory
use m_globals
use m_config_io
use mkl

implicit none

character(len=256) :: cla_buf
character(len=:), allocatable :: fn_sf, fn_sf_bb, fn_sf_sc
real(rp), dimension(:), allocatable :: r_bin !Binned |r|
real(rp), dimension(:), allocatable :: gr_bin !Binned pdf
real(rp), dimension(:), allocatable :: normalization
real(rp), dimension(3) :: ri, rj, rij
real(rp) :: rij_mag, rmin, rmax, bin_size
real(rp) :: rho
integer, dimension(:), allocatable :: count_bin
integer, dimension(:), allocatable :: sc_pos
integer :: natoms, numframes
integer :: num_bins
integer :: teq
integer :: istep, iatm, jatm, ibin, iframe, ierr, ios
integer :: fu_sf, fu_sf_bb, fu_sf_sc
integer :: ibr, ibeg, iend, ia_beg, ia_end
integer :: i, j, k
integer(ip_long) :: mode
real(rp) :: lbox
integer :: mflag

call get_command_argument(1, cla_buf)
fn_cfg = trim(adjustl(cla_buf))

call get_command_argument(2, cla_buf)
fn_traj = trim(adjustl(cla_buf))

call get_command_argument(3, cla_buf)
fn_sf = trim(adjustl(cla_buf))

call get_command_argument(4, cla_buf)
num_bins = int(str_to_d(cla_buf))

call get_command_argument(5, cla_buf)
rmax = str_to_d(cla_buf)

!call get_command_argument(6, cla_buf) 
!natoms = int(str_to_d(cla_buf))

call get_command_argument(6, cla_buf)
teq = int(str_to_d(cla_buf))

call get_command_argument(7, cla_buf)
lbox = str_to_d(cla_buf)

bin_size = rmax/num_bins

!write(*,'(a,t15,g0.8)') 'BOXL', lbox
!write(*,'(a,t15,g0.8)') 'RMAX', rmax
!write(*,'(a,t15,g0.8)') 'NATOMS', natoms
!write(*,'(a,t15,g0.8)') 'BIN_SIZE', bin_size
!write(*,'(a,t15,g0.8)') 'NUM_BINS', num_bins
!write(*,'(a,t15,g0.8)') 'NUM_TRAJ', num_traj

allocate(r_bin(1:num_bins))
allocate(gr_bin(1:num_bins))
allocate(normalization(1:num_bins))

numframes = 0
gr_bin = 0.0
rmin = 0.0

!Read final (can be initial as well) configuration file
call read_config(fn_cfg)

natoms = num_atoms
rho = natoms/lbox**3


do i = 1, num_bins
    r_bin(i) = (i-1)*(rmax-rmin)/num_bins
end do

do j = 1, num_bins-1
    normalization(j) = ((r_bin(j+1))**3-(r_bin(j))**3)*(4.0/3.0*math_pi*rho)*(natoms-1)/2
end do

normalization(num_bins) = (rmax**3-(rmax-bin_size)**3)*(4.0/3.0*math_pi*rho)*(natoms-1)/2

!write(*,*) 'rbin', r_bin
!write(*,*) 'grbin', gr_bin
!write(*,*) 'normalization', normalization

mflag = 1
mode =  vmlSetMode(VML_LA)


!Open file for writing sf data & write header line
open(newunit=fu_sf, file=fn_sf, action='write', status='replace')
write(fu_sf, '(a16,2x,a16)') 'r', 'gr'

!Open file for reading trajectory
call traj%open(fn_traj, 'r', ierr)
write(*,*) 'numframes ', traj%num_frames
do iframe = 1, traj%num_frames
    call traj%read(iframe, nts, ierr, mflag, coordinates)
    !write(*,*) 'nts ', nts
    if (nts < teq) cycle
    write(*,*) 'iframe ', iframe
    numframes = numframes + 1
    do iatm = 1, num_atoms-1
        ri = coordinates(:,iatm)
        do jatm = iatm+1, num_atoms
            rj = coordinates(:,jatm)
            rij = rj - ri
            rij(1) = rij(1) - nint(rij(1)/lbox)*lbox
            rij(2) = rij(2) - nint(rij(2)/lbox)*lbox
            rij(3) = rij(3) - nint(rij(3)/lbox)*lbox            
            rij_mag = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))
            if (rij_mag > rmax) cycle
            ibin = floor(rij_mag/bin_size) + 1
            !r_bin(ibin) = rij(ibin) + rij_mag
            gr_bin(ibin) = gr_bin(ibin) + 1
        end do
    end do
end do

!r_bin = r_bin/gr_bin
gr_bin = gr_bin/normalization/numframes

!write(*,*) 'gr_bin', gr_bin
do i = 1, num_bins
    write(fu_sf, '(es16.7,2x,es16.7)') r_bin(i), gr_bin(i)
end do

!Close files
close(fu_sf)

!*******************************************************************************

end program
