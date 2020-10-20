module setup_m
    !! Routines for doing allocation, etc. in preparation for simulation run.

use constants_m
use strings_m
use random_m
use logger_m
use simbox_m
use atmcfg_m
use config_io_m
use control_m
use stats_m
use interaction_m
use trajectory_m
use bd_solver_m

implicit none

contains

!*******************************************************************************

subroutine run(cpar, job_tag)

    type(ctrlpar_t),  intent(in) :: cpar
    character(len=*), intent(in) :: job_tag
        !! This string will be used as is, i.e. if there are leading or trailing
        !! white spaces, they will appear in any I/O file names.
    character(len=:), allocatable :: fn
    character(len=128) :: msg
    integer(ip_long) :: ccbeg, ccend, crate
    integer(ip_long) :: nts
    type(smbx_t)     :: simbox
    type(atmcfg_t)   :: atc
    integer :: ierr
    logical :: lexists

    ierr =0

    !Read configuration from appropriate files.
    !Note: Reading from file includes initialization of the simulation box and
    !atom configuration.
    if ( cpar%lrevive ) then
        !Restarting simulation: Read revive file
        fn = cpar%fn_revive//job_tag
        call read_dump(nts, simbox, atc, fn)
        call logger%log_msg('read revive file `'//fn//'`')
    else 
        !New simulation: Read initial configuration file
        fn = cpar%fn_cfg//job_tag
        call read_config(simbox, atc, fn)
        call logger%log_msg('read config file `'//fn//'`')
        nts = 0
    end if

    !Allocate memory for forces. Forces are not saved in revive file or
    !in config file as they can be calculated from position data.
    allocate( atc%forces(3,atc%num_atoms) )
    atc%forces = 0.0_rp

    !Initialize stats collection
    call stats_init(cpar, simbox, atc, job_tag)
    call logger%log_msg('initialized stats collection')

    !Set up interactions
    call ia_setup(cpar, simbox, atc)
    call logger%log_msg('set up interaction params')

    !Set up trajectory if necessary. There can be only one trajectory per
    !simulation.
    if (cpar%write_traj) then
        fn = cpar%fn_traj//job_tag
        if (cpar%lrevive) then
            !Append/write to existing/new trajectory file
            !Check if trajectory file exists
            inquire(file=fn, exist=lexists)
            if (lexists) then
                !Open existing file for appending
                call traj%init(fn, 'rw', ierr)
                if (ierr /=0 ) then
                    call logger%log_msg('failed to open trajectory file ' &
                        //'`'//fn//'`')
                    call finish(atc); return
                else
                    call logger%log_msg('trajectory file opened')
                end if
            else
                !Create new trajectory file
                call traj%init(fn, atc%num_atoms, [integer:: 1, 0, 0, 0])
                call logger%log_msg('trajectory file created')
            end if
        else
            !Create new trajectory file
            call traj%init(fn, atc%num_atoms, [integer:: 1, 0, 0, 0])
            call logger%log_msg('trajectory file created')
        end if
    end if

    !Initialize random number generator
    if (cpar%read_seed) then
        call init_stream('random_seed.txt'//job_tag)
    else 
        call init_stream('')
    end if
    if (cpar%write_seed) call save_seed('random_seed.txt'//job_tag)
    call logger%log_msg('rng initialized')

    !Set up solver.
    call bds_init(cpar, atc%num_atoms, job_tag, ierr)
    if (ierr /= 0) then
        call logger%log_msg('bd solver initialization failed')
        call finish(atc); return
    end if
    call logger%log_msg('bd solver initialized')

    call system_clock(ccbeg, crate)

    call bds_run(nts, simbox, atc, ierr)

    call system_clock(ccend, crate)

    write(msg,'(a,es12.5)') 'execution time(s) = ', (ccend-ccbeg)/real(crate,rp)

    call logger%log_msg(msg)

    call finish(atc)

    end subroutine

!*******************************************************************************

subroutine finish(atc)

    type(atmcfg_t), intent(in out) :: atc

    call bds_finish()
    call traj%delete()
    call ia_finish()
    call stats_finish()
    call atmcfg_delete(atc)

    end subroutine

!*******************************************************************************

end module setup_m
