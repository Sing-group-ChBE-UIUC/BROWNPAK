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

module control_m
    !! Routines for reading and witing control file.

use constants_m
use strings_m

implicit none

integer, parameter :: mxrdln = 1024
    !! Maximum length of character string for input line buffer.

type ctrlpar_t
    real(rp) :: rcutoff = 0.0_rp
    real(rp) :: tskin = 0.0_rp
    character(len=5) :: mth_ptgen = ''
    logical  :: lelst = .false.
    logical  :: lhdia = .false.
    logical  :: lvdw = .false.
    integer  :: excluded_atoms = 0
    character(len=4) :: mob_fctr = ''
    integer  :: lanc_mxitr = 0
    real(rp) :: lanc_tol = 0.0_rp
    integer  :: se_nlmxitr = 0
    integer  :: se_kdmax = 0
    character(len=4) :: bdintg = ''
    real(rp), dimension(2) :: se_tol = 0.0_rp
    integer, dimension(2) :: stats_binsize = 0

    real(rp) :: tim_stp = 0.0_rp
    integer(ip_long) :: nts_sim = 0
    integer          :: nts_mobsam = 0
    integer(ip_long) :: nts_dump = 0
    integer(ip_long) :: nts_samp = 0
    integer(ip_long) :: nts_log = 0


    logical :: lrevive = .false.
        !! {T, F}. Whether the simulation is restarted.
    logical :: read_seed = .false.
        !! {T, F}. Whether to initialize the random number generator by reading
        !! a seed from a file. If `read_seed` == T, the seed will be read from
        !! a file 'random_seed.txt'
    logical :: write_seed = .false.
        !! {T, F}. Whether to write the random number generator seed. If
        !!  `write_seed` == T the seed will be written to a file named
        !!  'random_seed.txt'
    logical :: write_traj = .false.
        !! Should the trajectory be written to file? {T, F}
    character(len=:), allocatable :: fn_cfg
        !! Name of the file containing the initial configuration
    character(len=:), allocatable :: fn_revive
        !! Name of the revive file
    character(len=:), allocatable :: fn_stats
        !! Name of the statistics file
    character(len=:), allocatable :: fn_traj
        !! Name of the trajectory file
    contains
        procedure :: read => control_read
        procedure :: write => control_write
end type ctrlpar_t

contains

!******************************************************************************

subroutine control_read(this, fn)
    !! Reads simulation control parameters from file

    class(ctrlpar_t), intent(out) :: this
        !! A *ctrlpar_t* instance.
    character(len=*), intent(in) :: fn
        !! Name of parameters file.
    character(len=:), allocatable :: key
    character(len=:), allocatable :: val
    character(len=mxrdln) :: line
    character(len=1) :: cstr = '#' !Comment string
    integer :: fu
    integer :: ios

    open (newunit=fu, file = fn, action = 'read', status = 'old')
    
    do 
        call readline(fu, line, cstr, ios)
        if (ios /= 0) return

        call str_get_keyval(line, key, val)

        if (key=='rcutoff')        read(val,*) this%rcutoff
        if (key=='tskin')          read(val,*) this%tskin
        if (key=='mth_ptgen')      read(val,*) this%mth_ptgen
        if (key=='lelst')          read(val,*) this%lelst
        if (key=='lhdia')          read(val,*) this%lhdia
        if (key=='lvdw')           read(val,*) this%lvdw
        if (key=='excluded_atoms') read(val,*) this%excluded_atoms
        if (key=='mob_fctr')       read(val,*) this%mob_fctr
        if (key=='lanc_mxitr')     read(val,*) this%lanc_mxitr
        if (key=='lanc_tol')       read(val,*) this%lanc_tol
        if (key=='bdintg')         read(val,*) this%bdintg
        if (key=='se_nlmxitr')     read(val,*) this%se_nlmxitr
        if (key=='se_kdmax')       read(val,*) this%se_kdmax
        if (key=='se_tol')         read(val,*) this%se_tol
        if (key=='stats_binsize') read(val,*) this%stats_binsize

        if (key=='tim_stp')    read(val,*) this%tim_stp
        if (key=='nts_sim')    this%nts_sim    = int(str_to_d(val), ip_long)
        if (key=='nts_mobsam') this%nts_mobsam = int(str_to_d(val))
        if (key=='nts_dump')   this%nts_dump   = int(str_to_d(val), ip_long)
        if (key=='nts_samp')   this%nts_samp   = int(str_to_d(val), ip_long)
        if (key=='nts_log')    this%nts_log    = int(str_to_d(val), ip_long)

        if (key == 'lrevive')    read(val,*) this%lrevive
        if (key == 'read_seed')  read(val,*) this%read_seed
        if (key == 'write_seed') read(val,*) this%write_seed
        if (key == 'write_traj') read(val,*) this%write_traj

        if (key == 'fn_cfg')    this%fn_cfg    = val
        if (key == 'fn_revive') this%fn_revive = val
        if (key == 'fn_stats')  this%fn_stats  = val
        if (key == 'fn_traj')   this%fn_traj   = val
    end do

    close (fu)

    end subroutine

!******************************************************************************

subroutine control_write(this, fn)
    !! Write simulation parameters to file

    class(ctrlpar_t), intent(out) :: this
        !! A *ctrlpar_t* instance.
    character(len=*), intent(in) :: fn
        !! Output file name
    integer :: fu

    open(newunit=fu, file=fn, action='write', status='replace')

    write(fu, '(a,t20,g0.6)') 'rcutoff', this%rcutoff
    write(fu, '(a,t20,g0.6)') 'tskin', this%tskin
    write(fu, '(a,t20,a)'   ) 'mth_ptgen', trim(adjustl(this%mth_ptgen))
    write(fu, '(a,t20,l1)'  ) 'lelst', this%lelst
    write(fu, '(a,t20,l1)'  ) 'lhdia', this%lhdia
    write(fu, '(a,t20,l1)'  ) 'lvdw', this%lvdw
    write(fu, '(a,t20,i0)'  ) 'excluded_atoms', this%excluded_atoms
    write(fu, '(a,t20,a)'   ) 'mob_fctr', this%mob_fctr
    write(fu, '(a,t20,i0)'  ) 'lanc_mxitr', this%lanc_mxitr
    write(fu, '(a,t20,a)'   ) 'lanc_tol', str_from_num(this%lanc_tol,'(es24.15)')
    write(fu, '(a,t20,a)'   ) 'bdintg', this%bdintg
    write(fu, '(a,t20,i0)'  ) 'se_nlmxitr', this%se_nlmxitr
    write(fu, '(a,t20,i0)'  ) 'se_kdmax', this%se_kdmax
    write(fu, '(a,t20,es12.5,1x,es12.5)') 'se_tol', this%se_tol
    write(fu, '(a,t20,i0,1x,i0)') 'stats_binsize', this%stats_binsize

    write(fu, *)
    write(fu, '(a,t20,a)')  'tim_stp',    str_from_num(this%tim_stp,'(es24.15)')
    write(fu, '(a,t20,i0)') 'nts_sim',    this%nts_sim
    write(fu, '(a,t20,i0)') 'nts_mobsam', this%nts_mobsam
    write(fu, '(a,t20,i0)') 'nts_dump',   this%nts_dump
    write(fu, '(a,t20,i0)') 'nts_samp',   this%nts_samp
    write(fu, '(a,t20,i0)') 'nts_log',    this%nts_log

    write(fu, *)
    write(fu, '(a,t20,l1)') 'lrevive',    this%lrevive
    write(fu, '(a,t20,l1)') 'read_seed',  this%read_seed
    write(fu, '(a,t20,l1)') 'write_seed', this%write_seed
    write(fu, '(a,t20,l1)') 'write_traj', this%write_traj

    write(fu, *)
    if (allocated(this%fn_cfg)) then
        write(fu,'(a,t20,a)') 'fn_cfg', this%fn_cfg
    else
        write(fu,'(a)') 'fn_cfg'
    end if

    if (allocated(this%fn_revive)) then
        write(fu,'(a,t20,a)') 'fn_revive', this%fn_revive
    else
        write(fu,'(a)') 'fn_revive'
    end if

    if (allocated(this%fn_stats)) then
        write(fu,'(a,t20,a)') 'fn_stats', this%fn_stats
    else
        write(fu,'(a)') 'fn_stats'
    end if

    if (allocated(this%fn_traj)) then
        write(fu, '(a,t20,a)') 'fn_traj',   this%fn_traj
    else
        write(fu, '(a)') 'fn_traj'
    end if

    close(fu)
    
    end subroutine

!******************************************************************************

end module control_m
