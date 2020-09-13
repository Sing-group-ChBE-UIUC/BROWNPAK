module m_bd_solver
    !! Routines implementing Brownian Dynamics (BD) solver.

use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
use mkl_blas
use mkl_rci

use m_precision
use m_constants_math
use m_logger
use m_globals
use m_interaction
use m_brown
use m_stats_io
use m_config_io

implicit none

private
public :: bds_init, bds_run, bds_finish

integer  :: cntr_mobsam
real(rp) :: sqrt_two_dt
integer  :: nlmxitr = 150
integer  :: lmxitr = 150
real(rp) :: tol_res = 1.0e-6_rp
real(rp) :: tol_upd = 1.0e-8_rp
real(rp) :: atol = 0.0_rp
real(rp) :: rtol = 1.0e-6_rp

real(rp), dimension(:  ), allocatable :: drift
real(rp), dimension(:,:), allocatable :: diffusion
real(rp), dimension(:,:), allocatable :: mob

real(rp), dimension(:), allocatable :: crd0
real(rp), dimension(:), allocatable :: crdk
real(rp), dimension(:), allocatable :: fvalk
real(rp), dimension(:), allocatable :: fval
real(rp), dimension(:), allocatable :: rhs
real(rp), dimension(:), allocatable :: sol
real(rp), dimension(:), allocatable :: fgmres_tmp
real(rp), dimension(:), allocatable :: crdn
real(rp), dimension(:), allocatable :: h2
real(rp), dimension(:), allocatable :: h3
real(rp), dimension(:), pointer :: pvcrd
real(rp), dimension(:), pointer :: pvfrc

contains
 
!*******************************************************************************
 
subroutine bds_init(ierr)
    !! Initializes the BD solver.

    integer, intent(out) :: ierr
    integer :: ltmp
    logical :: lexists
    type(c_ptr) :: loc

    ierr = 0
    sqrt_two_dt = sqrt(2.0_rp*tim_stp)

    !Assign pointers
    loc = c_loc(coordinates)
    call c_f_pointer(loc, pvcrd, [3*num_atoms])
    loc = c_loc(forces)
    call c_f_pointer(loc, pvfrc, [3*num_atoms])

    !Allocate memory for BD drift-diffusion SDE
    if (bdintg=='EM') then
        allocate(drift (3*num_atoms))
        if (lhdia) then
            allocate(mob(3*num_atoms,3*num_atoms))
            allocate(diffusion(3*num_atoms,nts_mobsam))
            if (mob_fctr == 'CHOL') then
                call brn_init(3*num_atoms, nts_mobsam, mob_fctr)
            else if (mob_fctr == 'LANC') then
                call brn_init(3*num_atoms, nts_mobsam, mob_fctr, &
                            lanc_mxitr, lanc_tol)
            end if
            cntr_mobsam = 0
        else
            allocate(diffusion(3*num_atoms,1))
        end if
    end if

    if (bdintg == 'SE') then
        allocate(crd0(3*num_atoms))
        allocate(crdk(3*num_atoms))
        allocate(fvalk(3*num_atoms))
        allocate(fval(3*num_atoms))
        allocate(rhs(3*num_atoms))
        allocate(sol(3*num_atoms))
        ltmp = (2*lmxitr+1)*3*num_atoms + lmxitr*(lmxitr+9)/2 + 1
        allocate(fgmres_tmp(ltmp))
        allocate(drift (3*num_atoms))
        if (lhdia) then
            allocate(mob(3*num_atoms,3*num_atoms))
            allocate(diffusion(3*num_atoms,nts_mobsam))
            if (mob_fctr == 'CHOL') then
                call brn_init(3*num_atoms, nts_mobsam, mob_fctr)
            else if (mob_fctr == 'LANC') then
                call brn_init(3*num_atoms, nts_mobsam, mob_fctr, &
                            lanc_mxitr, lanc_tol)
            end if
            cntr_mobsam = 0
        else
            allocate(diffusion(3*num_atoms,1))
        end if
    end if

    if (bdintg == 'RK') then
        !Works only with `mob_fctr = CHOL`
        allocate(crdn(3*num_atoms))
        allocate(h2(3*num_atoms))
        allocate(h3(3*num_atoms))
        allocate(drift(3*num_atoms))
        if (lhdia) then
            allocate(mob(3*num_atoms,3*num_atoms))
            allocate(diffusion(3*num_atoms,2*nts_mobsam))
            call brn_init(3*num_atoms, 2*nts_mobsam, 'CHOL')
            cntr_mobsam = 0
        else
            allocate(diffusion(3*num_atoms,2))
        end if
    end if

    !Opening trajectory file
    if (write_traj) then
        if ( lrevive ) then
            !Append (write) to existing (new) trajectory file
            !Check if trajectory file exists
            inquire(file=fn_traj//trim(adjustl(job_tag)), exist=lexists)
            if (lexists) then
                !Open existing file for appending
                call traj%open(fn_traj//trim(adjustl(job_tag)), 'rw', ierr)
                if (ierr /= 0) return
            else
                !Create new trajectory file
                call traj%create(fn_traj//trim(adjustl(job_tag)), num_atoms, &
                    0, traj_frmcmp)
            end if
        else
            !Create new trajectory file
            call traj%create(fn_traj//trim(adjustl(job_tag)), num_atoms, &
                0, traj_frmcmp)
        end if
    end if

    end subroutine

!******************************************************************************

subroutine bds_finish()
    !! Clears up memory allocated in [[bds_init]].

    call brn_finish()
    call traj%close()

    if (allocated(mob)) deallocate(mob)
    if (allocated(drift)) deallocate(drift)
    if (allocated(diffusion)) deallocate(diffusion)

    if (allocated(crd0)) deallocate(crd0)
    if (allocated(crdk)) deallocate(crdk)
    if (allocated(fvalk)) deallocate(fvalk)
    if (allocated(fval)) deallocate(fval)
    if (allocated(rhs)) deallocate(rhs)
    if (allocated(sol)) deallocate(sol)
    if (allocated(fgmres_tmp)) deallocate(fgmres_tmp)

    if (allocated(crdn)) deallocate(crdn)
    if (allocated(h2)) deallocate(h2)
    if (allocated(h3)) deallocate(h3)

    end subroutine

!******************************************************************************

subroutine bds_run()
    !! Driver for BD integrator.
    !!
    !! Repeatedly calls [[bds_integrate_fd]] or [[bds_integrate_hi]] to update
    !! atom positions.

    real(rp), dimension(3) :: com
    real(rp), dimension(3,3) :: delta
    integer :: ierr
    integer :: counter

    ierr = 0

    !For isolated untethered molecule, ensure c.o.m. is at the center of the
    !box.
    if ( (num_tethers == 0) .and. (imcon == 0) ) then
        call simbox%to_center(coordinates)
    end if

    !Log & dump starting configuration
    call logger%log_msg('nts = '//str_from_num(nts))
    call write_dump(fn_revive//trim(adjustl(job_tag)))

    counter = 0
    do while (nts < nts_sim)
        if (bdintg=='EM') then
            call integrate_em(ierr)
        else if (bdintg=='SE') then
            call integrate_se(ierr)
        else if (bdintg=='RK') then
            call integrate_rk(ierr)
        end if
        if (ierr /= 0) return

        !Accumulate stress
        counter = counter + 1
        delta = stress - stress_accu
        stress_accu = stress_accu + (delta/counter)
        nts = nts + 1

        !Apply boundary conditions: For isolated untethered molecule,
        !update molecule_com & bring c.o.m. to the center of the box.
        if ( (num_tethers == 0) .and. (imcon == 0) ) then
            call simbox%to_center(coordinates, com)
            molecule_com = molecule_com + com
        end if

        !For PBC, wrap atom positions
        if (imcon /= 0) call simbox%wrap_all(coordinates)

        !Logging
        if (mod(nts,nts_log) == 0) then
            call logger%log_msg('nts = '//str_from_num(nts))
        end if

        !Dump revive file
        if (mod(nts,nts_dump)==0) then
            call write_dump(fn_revive//trim(adjustl(job_tag)))
        end if

        if (mod(nts,nts_samp)==0) then
            !Write stats
            call stats_write()
            !Write traj
            if (write_traj) call traj%append_frame(nts, coordinates, &
                velocities, forces, charge)
            !Empty stress accumulator
            stress_accu = 0.0_rp; counter = 0
        end if
    end do

    !Dump the final configuration
    call write_dump(fn_revive//trim(adjustl(job_tag)))

    end subroutine

!******************************************************************************

subroutine integrate_em(ierr)
    !! Performs one step of BD integration using explicit Euler-Maruyama scheme.

    integer, intent(out) :: ierr

    ierr = 0
    !Calculate the diffusion term
    call calc_diffusion(ierr)
    if (ierr /= 0) return

    !Calculate the drift term
    call calc_drift(ierr)
    if (ierr /= 0) return

    !Update positions
    if (lhdia) then
        pvcrd = pvcrd + drift + diffusion(:,cntr_mobsam)
        if (cntr_mobsam == nts_mobsam) cntr_mobsam = 0
    else
        pvcrd = pvcrd + drift + diffusion(:,1)
    end if

    end subroutine

!******************************************************************************

subroutine integrate_se(ierr)
    !! Performs one step of BD integration using semi-implicit Euler scheme.

    integer, intent(out) :: ierr
    real(rp) :: err_res, err_upd
    real(rp) :: nrm2_fval0
    integer :: k, itercount
    logical :: lconv
    character(len=256) :: msg

    ierr = 0; lconv = .false.

    !Calculate the diffusion term
    call calc_diffusion(ierr)
    if (ierr /= 0) return

    !Evaluate the nonlinear function arising out of discretizing the SDE in
    !semi-implicit form
    crd0 = pvcrd !Coordinates at the end of last time step
    call se_calc_fval(ierr)
    if (ierr /= 0) return

    !Save the norm of the function value
    nrm2_fval0 = norm2(fval)

    !Update coordinates via NR iterations
    do k = 1, nlmxitr
        !Save the coordinates and the function value
        crdk = pvcrd; fvalk = fval
        !Solve the linear system
        call se_linsolve(lconv, itercount, ierr)
        if (lconv) then
            !write(msg,'(a,i0,a)')  &
            !    'fgmres iterations at convergence = ', itercount
            !call logger%log_msg(msg)
        else
            write(msg,'(a)') 'fgmres failed to converge'
            call logger%log_msg(msg)
        end if
        if (ierr /= 0) return
        !Update the coordinates
        pvcrd = crdk + 0.25_rp*sol
        !Evaluate the function
        call se_calc_fval(ierr)
        if (ierr /= 0) return
        !Calculate errors
        err_upd = norm2(sol)/norm2(pvcrd)
        err_res = norm2(fval)/nrm2_fval0
        !write(*,*) 'err_upd = ', err_upd, ' err_res = ', err_res
        lconv = .false.
        if ((err_upd<=tol_upd) .or. (err_res<=tol_res)) then
            lconv = .true.
            !write(msg,'(a,i0)') 'newton iterations at convergence = ', k
            !call logger%log_msg(msg)
            exit
        end if
    end do

    if (.not. lconv) then
        write(msg,'(a,es12.5)') 'newton iterations not converged. err_res = ', &
                                err_res
        call logger%log_msg(msg)
        ierr = 1
    end if

    if (lhdia) then
        if (cntr_mobsam == nts_mobsam) cntr_mobsam = 0
    end if

    end subroutine

!******************************************************************************

subroutine se_calc_fval(ierr)

    integer, intent(out) :: ierr

    ierr = 0; fval = 0.0_rp

    call calc_drift(ierr)
    if (ierr /= 0) return

    if (lhdia) then
        fval = pvcrd - drift - crd0 - diffusion(:,cntr_mobsam)
    else
        fval = pvcrd - drift - crd0 - diffusion(:,1)
    end if

    end subroutine

!******************************************************************************

subroutine se_linsolve(lconv, itercount, ierr)
    !! Solves the linear system using fgmres.

    logical, intent(out) :: lconv
    integer, intent(out) :: itercount
    integer, intent(out) :: ierr
    real(rp) :: nrm2crdk, nrm2v, epsn, eps, reps
    integer :: ibeg, iend
    integer :: n, rci_request
    integer, dimension(128):: ipar
    real(rp), dimension(128):: dpar
    character(len=256) :: msg

    lconv = .false.; ierr = 0
    nrm2crdk = norm2(crdk)
    epsn = sqrt(epsilon(0.0_rp)*(1.0_rp + nrm2crdk))
    n = 3*num_atoms
    sol = 1.0_rp; rhs = -fvalk; fgmres_tmp = 0.0_rp

    !Initialize fgmres solver
    call dfgmres_init(n, sol, rhs, rci_request, ipar, dpar, fgmres_tmp)
    if (rci_request < 0) then
        call logger%log_msg('<dfgmres_init> failed.')
    end if
    ipar(2) = logger%fu
    ipar(9) = 1 !Perform residual stopping test
    ipar(10) = 0 !No user-defined stopping tests
    ipar(12) = 1 !Test for zero norm of the next generated vector
    ipar(15) = lmxitr

    dpar(1) = rtol !Relative tolerance
    dpar(2) = atol !Absolute tolerance

    !Preliminary check for fgmres solver
    call dfgmres_check(n, sol, rhs, rci_request, ipar, dpar, fgmres_tmp)
    if (rci_request < 0) then
        call logger%log_msg('<dfgmres_check> failed.')
    end if

    do
        call dfgmres(n, sol, rhs, rci_request, ipar, dpar, fgmres_tmp)
        if (rci_request == 0) then
            !solution found
            call dfgmres_get(n, sol, rhs, rci_request, ipar, dpar, fgmres_tmp,&
                itercount)
            lconv = .true.
            exit
        else if (rci_request == 1) then
            !Perform matrix-vector multiplication
            ibeg = ipar(22); iend = ibeg + n - 1
            nrm2v = norm2(fgmres_tmp(ibeg:iend))
            eps = epsn/nrm2v; reps = 1.0_rp/eps
            !print*, eps
            pvcrd = crdk + eps*fgmres_tmp(ibeg:iend)
            call se_calc_fval(ierr)
            ibeg = ipar(23); iend = ibeg + n - 1
            fgmres_tmp(ibeg:iend) = (fval-fvalk)*reps
    !print*, fgmres_tmp(ibeg:iend)
    !stop
        else
            !fgmres failed
            write(msg,'(a,i0)') '<dfgmres> failed, err code = ', rci_request
            call logger%log_msg(msg)
            ierr = 1
            exit
        end if
    end do

    end subroutine

!******************************************************************************

subroutine integrate_rk(ierr)
    !! Performs one step of BD integration using explicit RK scheme.

    integer, intent(out) :: ierr
    real(rp), parameter :: sixth = 1.0_rp/6.0_rp
    real(rp), parameter :: two_third = 2.0_rp/3.0_rp
    integer :: j

    ierr = 0; diffusion = 0.0_rp

    !Early return for no HI
    if (.not. lhdia) then
        !Calculate standard Gaussian vectors
        call brn_calc_dw(diffusion)
        diffusion = sqrt(tim_stp)*diffusion
        diffusion(:,nts_mobsam+1:) = diffusion(:,nts_mobsam+1:)/math_sqrt3
        diffusion(:,nts_mobsam+1:) = ( diffusion(:,nts_mobsam+1:) &
                + diffusion(:,1:nts_mobsam) )*0.5_rp*tim_stp
        return
    end if

    !Diffusion term with HI
    if (cntr_mobsam == 0) then
        !Calculate the mobility matrix
        call calc_rpy_tensor(mob)
        !Calculate standard Gaussian vectors
        call brn_calc_dw(diffusion)
        diffusion = sqrt(tim_stp)*diffusion
        diffusion(:,nts_mobsam+1:) = diffusion(:,nts_mobsam+1:)/math_sqrt3
        diffusion(:,nts_mobsam+1:) = ( diffusion(:,nts_mobsam+1:) &
                + diffusion(:,1:nts_mobsam) )*0.5_rp*tim_stp
        !Calculate the Brownian vectors
        call brn_calc_bdw(mob, diffusion, ierr)
        if (ierr /= 0) return
        !Cholesky decomposition overwrites the upper triangular part of mob
        !with the factorization result U. Putting the diagonal elements
        !back to unity (for RPY) as U is no longer needed. The lower
        !triangular part of mob will be used in calculating drift.
        do j = 1, 3*num_atoms
            mob(j,j) = 1.0_rp
        end do
    end if
    cntr_mobsam = cntr_mobsam + 1

    h2 = pvcrd; h3 = pvcrd; crdn = pvcrd
    call calc_drift(ierr)
    if (ierr /= 0) return
    h2 = h2 + drift; h3 = h3 + 0.25_rp*drift; crdn = crdn + sixth*drift

    pvcrd = h2
    call calc_drift(ierr)
    if (ierr /= 0) return
    h3 = h3 + 0.25_rp*drift; crdn = crdn + sixth*drift

    h3 = h3 + diffusion(:,nts_mobsam+cntr_mobsam)*sqrt_two_dt*1.5_rp
    pvcrd = h3
    call calc_drift(ierr)
    if (ierr /= 0) return
    crdn = crdn + two_third*drift + sqrt_two_dt*diffusion(:,cntr_mobsam)
    pvcrd = crdn

    if (lhdia) then
        if (cntr_mobsam == nts_mobsam) cntr_mobsam = 0
    end if

    end subroutine

!******************************************************************************

subroutine calc_diffusion(ierr)
    !! Calculates the diffusion term of the SDE. Updates module variables
    !!`diffusion` and `cntr_mobsam`.

    integer, intent(out) :: ierr
    integer :: j, f
    logical :: lconv
    character(len=256) :: msg

    diffusion = 0.0_rp; ierr = 0

    !Early return for no HI
    if (.not. lhdia) then
        !Calculate standard Gaussian vectors
        call brn_calc_dw(diffusion)
        diffusion = sqrt_two_dt*diffusion
        return
    end if

    !Diffusion term with HI
    if (cntr_mobsam == 0) then
        !Calculate the mobility matrix
        call calc_rpy_tensor(mob)
        !Calculate standard Gaussian vectors
        call brn_calc_dw(diffusion)
        !Calculate the Brownian vectors
        if (mob_fctr == 'CHOL') then
            call brn_calc_bdw(mob, diffusion, ierr)
            if (ierr /= 0) return
            !Cholesky decomposition overwrites the upper triangular part of mob
            !with the factorization result U. Putting the diagonal elements
            !back to unity (for RPY) as U is no longer needed. The lower
            !triangular part of mob will be used in calculating drift.
            do j = 1, 3*num_atoms
                mob(j,j) = 1.0_rp
            end do
        else if (mob_fctr == 'LANC') then
            call brn_calc_bdw(mob, diffusion, ierr, lconv, f)
            if (ierr /= 0) return
            if (.not. lconv) then
                write(msg,'(a,i0)') '<brn_calc_bdw> not converged. nts = ', nts
                call logger%log_msg(msg)
            end if
        end if
        diffusion = sqrt_two_dt*diffusion
    end if
    cntr_mobsam = cntr_mobsam + 1

    end subroutine

!******************************************************************************

subroutine calc_drift(ierr)
    !! Calculates the drift term of the SDE.

    integer, intent(out) :: ierr
    integer :: i

    drift = 0.0_rp; ierr = 0
    
    !Update forces
    call ia_calc_forces(ierr)
    if (ierr /= 0) return

    !Calculate ambient flow velocity at the atom locations.
    select case(flow_style)
    case(1)
        !Steady simple shear: Flow along x, gradient along y
        do i = 1, num_atoms
            drift(3*i-2) = flow_params(1)*coordinates(2,i)
        end do
    case(2)
        !Steady planar extension in x-y plane
        do i = 1, num_atoms
            drift(3*i-2) =  flow_params(1)*coordinates(1,i)
            drift(3*i-1) = -flow_params(1)*coordinates(2,i)
        end do
    case(3)
        !Steady uniaxial extension along x
        do i = 1, num_atoms
            drift(3*i-2) =  flow_params(1)*coordinates(1,i)
            drift(3*i-1) = -0.5_rp*flow_params(1)*coordinates(2,i)
            drift(3*i)   = -0.5_rp*flow_params(1)*coordinates(3,i)
        end do
    case default
        continue
    end select

    !drift <- mob * forces + drift. Only the lower triangular part of
    !mob is accessed.
    if (lhdia) then
        call dsymv('L', 3*num_atoms, 1.0_rp, mob, 3*num_atoms, &
                    pvfrc, 1, 1.0_rp, drift, 1)
    else
        drift = drift + pvfrc
    end if

    drift = tim_stp * drift

    end subroutine

!******************************************************************************

subroutine calc_rpy_tensor(mob)
    !! Calculates the RPY approximation to the mobility tensor.

    real(rp), dimension(:,:), intent(out) :: mob
        !! (3*num_atoms*, 3*num_atoms*) matrix; stores the mobility tensor.
    real(rp), dimension(3)   :: ri, rj, rij
    real(rp), dimension(3,3) :: matrpy
    real(rp) :: rijm, irijm, irijm2
    real(rp) :: c1, c2, consij
    integer :: i, j

    mob = 0.0_rp

    !Calculate the RPY tensor (in strictly upper triangular form)
    do j = 2, num_atoms
        rj = coordinates(:,j)
        do i = 1, (j-1)
            ri = coordinates(:,i)
            rij = rj - ri
            rijm = norm2(rij)
            irijm = 1.0_rp/rijm
            irijm2 = irijm*irijm

            if (rijm >= 2.0_rp) then
              C1 =  1.0_rp + (2.0_rp/3.0_rp)*irijm2
              C2 =  1.0_rp - 2.0_rp*irijm2
              consij = 0.75_rp*irijm
            else
              C1 = 1.0_rp - 9.0_rp*rijm/(32.0_rp)
              C2 = 3.0_rp*rijm/(32.0_rp)
              consij = 1.0_rp
            end if

            matrpy(1,1) = consij*( C1 + C2*rij(1)*rij(1)*irijm2 )
            matrpy(2,1) = consij*(      C2*rij(2)*rij(1)*irijm2 )
            matrpy(3,1) = consij*(      C2*rij(3)*rij(1)*irijm2 )

            matrpy(1,2) = consij*(      C2*rij(1)*rij(2)*irijm2 )
            matrpy(2,2) = consij*( C1 + C2*rij(2)*rij(2)*irijm2 )
            matrpy(3,2) = consij*(      C2*rij(3)*rij(2)*irijm2 )

            matrpy(1,3) = consij*(      C2*rij(1)*rij(3)*irijm2 )
            matrpy(2,3) = consij*(      C2*rij(2)*rij(3)*irijm2 )
            matrpy(3,3) = consij*( C1 + C2*rij(3)*rij(3)*irijm2 )

            mob(3*i-2:3*i,3*j-2:3*j) = matrpy
        end do
    end do
   
    !Copy strictly upper triangular part to strictly lower triangular part.
    do j = 2, 3*num_atoms
        do i = 1, j-1
            mob(j,i) = mob(i,j)
        end do
    end do

    !Put one on the diagonal.
    do j = 1, 3*num_atoms
       mob(j,j) = 1.0_rp
    end do

    end subroutine

!******************************************************************************

end module m_bd_solver
