program main

use m_precision
use m_constants_math
use m_trajectory
use m_globals
use m_config_io
use m_table
use m_qsort
use mkl

implicit none

type(itable_t) :: neighbor_table
type(itable_t) :: cluster_table
character(len=256) :: cla_buf
character(len=:), allocatable :: fn_sf, fn_probagg, fn_rg, fn_label
integer, dimension(:), allocatable :: label
real(rp), dimension(:,:), allocatable :: new_coordinates
real(rp), dimension(:,:), allocatable :: com_clusters
real(rp), dimension(:,:), allocatable :: molbuf
real(rp), dimension(:), allocatable :: r_bin
real(rp), dimension(:), allocatable :: gr_bin !Binned rdf
real(rp), dimension(:), allocatable :: normalization
! real(rp), dimension(:), allocatable :: prob_nn
real(rp), dimension(:), allocatable :: prob_aggnum
real(rp), dimension(:,:), allocatable :: dist_table
real(rp), dimension(:), allocatable :: kth_dist_table
real(rp), dimension(:), allocatable :: dist_buffer
integer, dimension(:), allocatable :: num_neighbor
real(rp), dimension(3) :: ri, rj, rij, com
real(rp), dimension(3) :: comi, comj, comij
real(rp) :: Rgsq
real(rp) :: rij_mag, rmax, rmin, rcrit
integer :: numframes
integer :: nchains, nbeads
integer :: teq
integer :: istep, iatm, jatm, ibin, iframe, ierr, ios
integer :: num_solvobead, num_totalbead
integer :: fu_sf, fu_probagg, fu_rg, fu_label
integer :: ibr, ibeg, iend, ia_beg, ia_end
integer :: i, j, k, c, idx
integer :: ichain, jchain
integer :: idx_k, minpts
integer :: num_nn
integer :: num_bins
integer :: cluster_count, ineighbor, neighbor_idx, nofn
integer :: agg_num, total_cluster_count
integer(ip_long) :: mode
real(rp) :: lbox
real(rp) :: bin_size
real(rp) :: rho
integer :: mflag
integer, dimension(:), pointer :: res



call get_command_argument(1, cla_buf)
fn_cfg = trim(adjustl(cla_buf))

call get_command_argument(2, cla_buf)
fn_traj = trim(adjustl(cla_buf))

call get_command_argument(3, cla_buf)
fn_sf = trim(adjustl(cla_buf))

call get_command_argument(4, cla_buf)
fn_probagg = trim(adjustl(cla_buf))

call get_command_argument(5, cla_buf)
fn_rg = trim(adjustl(cla_buf))

call get_command_argument(6, cla_buf)
teq = int(str_to_d(cla_buf))

call get_command_argument(7, cla_buf)
num_solvobead = int(str_to_d(cla_buf))

call get_command_argument(8, cla_buf)
num_totalbead = int(str_to_d(cla_buf))

call get_command_argument(9, cla_buf)
lbox = str_to_d(cla_buf)

call get_command_argument(10, cla_buf)
num_bins = str_to_d(cla_buf)

fn_label = 'label_all'

!Read final (can be initial as well) configuration file
call read_config(fn_cfg)
nchains = num_atoms/15
nbeads = nchains*num_solvobead
rmax = lbox/2
bin_size = rmax/num_bins

! allocate(r_bin(1:num_bins))
! allocate(gr_bin(1:num_bins))
! allocate(normalization(1:num_bins))

! allocate(prob_nn(1:100))

! allocate arrays for building distance table
allocate(dist_table(nbeads,nbeads))
! allocate arrays for building neighbor table 
allocate(num_neighbor(1:nbeads))
allocate(kth_dist_table(1:nchains))
allocate(dist_buffer(1:nchains))
! allocate arrays for recording cluster info
allocate(label(1:nchains))
allocate(prob_aggnum(1:nchains))
! allocate arrays for calculating RDF of clusters
allocate(r_bin(1:num_bins))
allocate(gr_bin(1:num_bins))
allocate(normalization(1:num_bins))

! initialize variables
numframes = 0
prob_aggnum = 0.0
total_cluster_count = 0
gr_bin = 0.0
r_bin = 0.0
normalization = 0.0
rmin = 0.0
rcrit = 4.2
idx_k = 4
minpts = 4

write(*,'(a,t15,g0.8)') 'BOXL', lbox
write(*,'(a,t15,g0.8)') 'RMAX', rmax
write(*,'(a,t15,g0.8)') 'NCHAINS', nchains
write(*,'(a,t15,g0.8)') 'BIN_SIZE', bin_size
write(*,'(a,t15,g0.8)') 'NUM_BINS', num_bins

allocate(new_coordinates(3,nbeads))
allocate(com_clusters(3,nchains))
allocate( molbuf(3,num_solvobead) )

mflag = 1

!Open file for writing sf data & write header line
open(newunit=fu_sf, file=fn_sf, action='write', status='replace')

!Open file for Rg writing
open(newunit=fu_rg, file=fn_rg, action='write', status='replace')

!Open file for reading trajectory
call traj%open(fn_traj, 'r', ierr)
write(*,*) 'numframes ', traj%num_frames
do iframe = 1, traj%num_frames
    call traj%read(iframe, nts, ierr, mflag, coordinates)
    !write(*,*) 'nts ', nts
    if (nts < teq) cycle
    write(*,*) 'iframe ', iframe
    numframes = numframes + 1
    k = 0
    do ichain = 1, nchains
        do iatm = 1, num_solvobead
            k = k + 1
            new_coordinates(:,k) = coordinates(:,num_totalbead*(ichain-1)+iatm)
        end do
    end do

    kth_dist_table = 0.0
    dist_buffer = 0.0
    dist_table = rmax
    label = -1
    num_neighbor = 0    
    ! build distance table for all atoms
    do i = 1, nbeads - 1 
        ri = new_coordinates(:,i)
        do j = i + 1, nbeads
            rj = new_coordinates(:,j)
            rij = rj - ri
            rij(1) = rij(1) - nint(rij(1)/lbox)*lbox
            rij(2) = rij(2) - nint(rij(2)/lbox)*lbox
            rij(3) = rij(3) - nint(rij(3)/lbox)*lbox            
            rij_mag = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))
            dist_table(j,i) = rij_mag
            dist_table(i,j) = rij_mag
        end do
    end do

    ! ! Build the neighbor table. Get the k_th nearest neighbor for each 
    ! ! "com of chain"(atom), and store the distance between this neighbor and 
    ! ! the atom under consideration
    ! call itbl_init(neighbor_table, nchains)
    ! do ichain = 1, nchains
    !     ! Find neighbors              
    !     do jchain = 1, nchains
    !         if (jchain == ichain) cycle
    !         if (dist_table(ichain,jchain) < rcrit) then 
    !             num_neighbor(ichain) = num_neighbor(ichain) + 1
    !             call neighbor_table%append(ichain,jchain)
    !         end if
    !     end do
    !     ! sort the distance table to find the k_th neighbor of each atom
    !     dist_buffer = dist_table(:,ichain)
    !     !write(*,*) dist_buffer
    !     call dqsort(dist_buffer)
    !     kth_dist_table(ichain) = dist_buffer(idx_k)
   
    ! end do
    ! ! write(*,*) kth_dist_table     

    ! create number of neighbor list
    call itbl_init(neighbor_table, nbeads)
    call itbl_init(cluster_table, nchains)
    do i = 1, nbeads
        ! Find neighbors              
        do j = 1, nbeads
            if (j == i) cycle
            if (dist_table(j,i) < rcrit) then 
                num_neighbor(i) = num_neighbor(i) + 1
                call neighbor_table%append(i,j)
            end if
        end do
    end do
    ! write(*,*) num_neighbor
    ! call neighbor_table%print

    ! Cluster counter
    cluster_count = 0
    do i = 1, nbeads
        ! Previously processed in inner loop
        if (label(i) /= -1) cycle
        ! Density check 
        if (num_neighbor(i) < minPts) then
            ! Label all beads in this chain as Noise
            ichain = (i/num_solvobead)*num_solvobead 
            do j = ichain+1, ichain+8
                label(j) = 0
            end do
            cycle
        end if
        ! next cluster label 
        cluster_count = cluster_count + 1   
        ! Label initial point
        label(i) = cluster_count                        
        ! Neighbors to expand 
        do k = 1, num_neighbor(i)
            neighbor_idx = neighbor_table%get_val(i,k)
            call cluster_table%append(cluster_count,neighbor_idx)
        end do
        idx = 1
        do while( .TRUE. )
            call cluster_table%get_row(cluster_count,res)
            if (size(res) < idx) exit
            neighbor_idx = cluster_table%get_val(cluster_count,idx)
            ! ! case if this neighbor is previously labeled as noise
            ! if (label(neighbor_idx) == 0) then
            !     ! change noise to border point
            !     label(neighbor_idx) = cluster_count     
            !     ! Find neighbors of this neighbor    
            !     if (num_neighbor(neighbor_idx) >= minpts) then
            !         do j = 1, num_neighbor(neighbor_idx)
            !             nofn = neighbor_table%get_val(neighbor_idx,j)
            !             if (cluster_table%is_in(cluster_count,nofn) .OR. nofn == i) cycle   
            !             call cluster_table%append(cluster_count,nofn)
            !         end do
            !     end if
            ! end if 
            ! case if this neighbor is already in the cluster
            idx = idx + 1
            if (neighbor_idx == i) cycle
            ! case if this neighbor is unlabeled
            if (label(neighbor_idx) == -1) then
                ! Label neighbor 
                label(neighbor_idx) = cluster_count              
                ! Find neighbors of this neighbor
                if (num_neighbor(neighbor_idx) >= minpts) then
                    do j = 1, num_neighbor(neighbor_idx)
                        nofn = neighbor_table%get_val(neighbor_idx,j)
                        if (cluster_table%is_in(cluster_count,nofn) .OR. nofn == ichain) cycle   
                        call cluster_table%append(cluster_count,nofn)
                    end do
                end if
            end if
        end do
    end do
    ! call cluster_table%print

    ! loop over all clusters, calculate statistics (aggregation number; Rg; etc.)
    do c = 1, cluster_count
        ! calculate aggregation number
        call cluster_table%get_row(c,res)
        agg_num = size(res)
        prob_aggnum(agg_num) = prob_aggnum(agg_num) + 1
        ! calculate COM of cluster based on COM's of chains
        ichain = cluster_table%get_val(c,1)
        com = new_coordinates(:,ichain)
        do i = 2, agg_num
            ichain = cluster_table%get_val(c,i)
            com = com + new_coordinates(:,ichain)
        end do
        com = com/agg_num
        ! record COM of this cluster in an array
        com_clusters(:,c) = com
        ! calculate Rg
        Rgsq = 0.0
        do i = 1, agg_num
            ichain = cluster_table%get_val(c,i)
            do iatm = 1, num_solvobead
                molbuf(:,iatm)= coordinates(:,num_totalbead*(ichain-1)+iatm)
                ri = molbuf(:,iatm) - com
                Rgsq = Rgsq + ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
            end do
        end do
        Rgsq = Rgsq/agg_num/num_solvobead
        write(fu_rg, '(es16.7)') Rgsq
    end do
    total_cluster_count = total_cluster_count + cluster_count
    
    ! loop over all clusters' COM, calculate RDF of clusters
    do c = 1, cluster_count - 1
        ri = com_clusters(:,c)
        do k = c + 1, cluster_count
            rj = com_clusters(:,k)
            rij = rj - ri
            rij(1) = rij(1) - nint(rij(1)/lbox)*lbox
            rij(2) = rij(2) - nint(rij(2)/lbox)*lbox
            rij(3) = rij(3) - nint(rij(3)/lbox)*lbox         
            rij_mag = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))
            if (rij_mag > rmax) cycle
            ibin = floor(rij_mag/bin_size) + 1
            gr_bin(ibin) = gr_bin(ibin) + 1
        end do
    end do

    open(newunit=fu_label, file=fn_label, action='write', status='replace')
    do i = 1, nchains
        do iatm = 1, num_solvobead
            k = (i-1)*num_solvobead+iatm
            rij = coordinates(:,num_totalbead*(i-1)+iatm)
            write(fu_label, '(i0,2x,i0,2x,*(g0.14,2x))') k, label(i)+1, rij
        end do
    end do
    close(fu_label)

    call neighbor_table%delete
    call cluster_table%delete

    ! probability function of neighbors, no longer used
    ! do ichain = 1, nchains
    !     ri = com_coordinates(:,ichain)
    !     num_nn = 0
    !     do jchain = 1, nchains
    !         if (jchain == ichain) cycle
    !         rj = com_coordinates(:,jchain)
    !         rij = rj - ri
    !         rij(1) = rij(1) - nint(rij(1)/lbox)*lbox
    !         rij(2) = rij(2) - nint(rij(2)/lbox)*lbox
    !         rij(3) = rij(3) - nint(rij(3)/lbox)*lbox            
    !         rij_mag = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))
    !         if (rij_mag < rcrit) then 
    !             num_nn = num_nn + 1
    !         end if
    !     end do
    !     prob_nn(num_nn) = prob_nn(num_nn) + 1
    ! end do
    ! do i = 1, nchains
    !     !write(*,*) com_coordinates(:,i)
    !     write(fu_sf, '(i0,2x,i0,2x,*(g0.14,2x))') i, mflag, com_coordinates(:,i)
    ! end do
    ! write(*,*) prob_nn
end do

!gr_bin = gr_bin/normalization/numframes
!prob_nn = prob_nn/nchains/numframes
prob_aggnum = prob_aggnum/total_cluster_count

do i = 1, num_bins
    r_bin(i) = (i-1)*(rmax-rmin)/num_bins
end do

rho = total_cluster_count/numframes/lbox**3

do j = 1, num_bins-1
    normalization(j) = ((r_bin(j+1))**3-(r_bin(j))**3)*(4.0/3.0*math_pi*rho)*(total_cluster_count/numframes-1)/2
end do
normalization(num_bins) = (rmax**3-(rmax-bin_size)**3)*(4.0/3.0*math_pi*rho)*(total_cluster_count/numframes-1)/2

gr_bin = gr_bin/normalization/numframes
! write(*,*) 'gr_bin', gr_bin

write(fu_sf, '(es16.7)') rho
do i = 1, num_bins
    write(fu_sf, '(es16.7,2x,es16.7)') r_bin(i), gr_bin(i)
end do

! do i = 1, nchains
!     ! write(fu_sf, '(I4,2x,es16.7)') i, kth_dist_table(i)
!     write(fu_sf, '(i0,2x,i0,2x,*(g0.14,2x))') i, label(i)+1, com_coordinates(:,i)
! end do

!Close files
close(fu_sf)

close(fu_rg)

open(newunit=fu_probagg, file=fn_probagg, action='write', status='replace')
do i = 1, 100
    write(fu_probagg, '(I4,2x,es16.7)') i, prob_aggnum(i)
end do

close(fu_probagg)

!*******************************************************************************

end program
