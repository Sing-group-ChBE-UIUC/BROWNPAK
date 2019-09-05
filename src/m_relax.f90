module m_relax
!! Performs structure relaxation (energy minimization) using steepest descent.
!!
!! Atom positions evolve following the gradient direction, in steps of size
!! `tim_stp`. If `tim_stp` is too large, bond constraints are likely to be
!! violated.
!!
!! If MPCD atoms are present, they do not take part in structure relaxation.

use m_precision
use m_ran_num
use m_strings
use m_globals
use m_interaction, only: ia_calc_forces
use m_stats_io
use m_config_io
use m_logger, only: logger => master_logger

implicit none

private

public :: rlx_run

contains

!******************************************************************************

subroutine rlx_run()
    !! Driver for relaxation

    character(len=:), allocatable :: fn_ldf, fn_rlx
    real(rp), parameter :: abs_tol = 1.0E-3_rp
    real(rp), parameter :: rel_tol = 1.0E-8_rp
    real(rp) :: enrg0, enrg_dif
    integer :: ierr
    logical :: zf, has_converged

    fn_rlx = str_strip(fn_cfg, '.cfg', 'r')//'-rlx.cfg'
    fn_ldf = str_strip(fn_cfg, '.cfg', 'r')//'-rlx.txt'

    ierr = 0

    !For isolated untethered molecule, ensure c.o.m. is at the center of the
    !box.
    if ( (num_tethers == 0) .and. (imcon == 0) ) then
        call simbox%to_center(coordinates)
    end if

    !Log & dump starting configuration
    call logger%info('rlx_run', ' nts: '//str_from_num(nts))
    call write_dump(fn_revive//trim(adjustl(job_tag)))

    !Relaxation run
    has_converged = .false.

    !Assume a very high initial energy
    enrg0 = huge(1.0_rp)

    do while (nts < nts_sim)
        if (has_converged) exit
        call rlx_integrate(ierr, zf)
        nts = nts + 1

        !For isolated untethered molecule bring c.o.m. to the center of the box
        if ( (num_tethers == 0) .and. (imcon == 0) ) then
            call simbox%to_center(coordinates)
        end if

        !Wrap around periodic boundaries
        if (imcon /= 0) call simbox%wrap_all(coordinates(:,1:num_atoms))

        !Logging
        if (mod(nts,nts_log) == 0) then
            call logger%info('rlx_run', ' nts: '//str_from_num(nts))
        end if

        !Dump revive file
        if (mod(nts,nts_dump)==0) then
            call write_dump(fn_revive//trim(adjustl(job_tag)))
        end if

        !Relaxation stats
        if (mod(nts,nts_samp)==0) call stats_write()

        !Convergence check
        !I am not using relative tolerance here.
        if (zf) has_converged = .true.
        enrg_dif = abs(energy_tot - enrg0)
        if ( enrg_dif <= (rel_tol*energy_tot+abs_tol) ) has_converged = .true.
        enrg0 = energy_tot
    end do

    !Log end of run
    if (ierr /= 0) then
        call logger%info('finish_rlx_run', 'Relaxation completed with violations')
    else
        call logger%info('finish_rlx_run', 'Relaxation completed')
    end if

    !Dump the final configuration
    call write_dump(fn_revive//trim(adjustl(job_tag)))

    !Export final configuration to config & ldf files
    call write_config(fn_rlx//trim(adjustl(job_tag)), ' ')
    call write_ldf(fn_ldf//trim(adjustl(job_tag)), ' ')

    end subroutine

!******************************************************************************

subroutine rlx_integrate(ierr, zf)
    !! Performs one step of relaxation

    integer, intent(out) :: ierr
    logical, intent(out) :: zf
    real(rp) :: scale_factor
    integer :: i

    zf = .false.

    !Calculate forces
    call ia_calc_forces(ierr)

    !If all forces are zero, nothing to do
    if ( all(forces == 0.0_rp) ) then
        zf = .true.
        return
    end if

    scale_factor = 1.0_rp/maxval(abs(forces))
    forces = scale_factor*forces

    scale_factor = 1.0_rp/sum(forces**2)
    forces = scale_factor*forces

    !Update position
    do i = 1, num_atoms
        coordinates(:,i) = coordinates(:,i) + tim_stp*forces(:,i)
    end do

    end subroutine

!********************************************************************************

end module m_relax
