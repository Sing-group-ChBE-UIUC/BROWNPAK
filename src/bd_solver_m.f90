module bd_solver_m
    !! Routines implementing Brownian Dynamics (BD) solver.

use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
use mkl_blas
!use mkl_rci

use constants_m
use logger_m
use control_m
use atmcfg_m
use interaction_m
use brown_m
use stats_m
use config_io_m
use trajectory_m

implicit none

private
public :: bds_init, bds_run, bds_finish, se_fval, se_jacv

logical  :: lhdia, write_traj
integer  :: lanc_mxitr, nlmxitr, kdmax, nts_mobsam
real(rp) :: tim_stp, lanc_tol, ftol, stptol
integer(ip_long) :: nts, nts_sim, nts_dump, nts_samp, nts_log
character(len=4) :: mob_fctr, bdintg
character(len=:), allocatable :: fn_revive

type(smbx_t), pointer :: psmbx => null()
type(atmcfg_t), pointer :: patc => null()

integer  :: cntr_mobsam
real(rp) :: sqrt_two_dt

real(rp), dimension(:  ), allocatable :: drift
real(rp), dimension(:,:), allocatable :: diffusion
real(rp), dimension(:,:), allocatable :: mob

real(rp), dimension(:), allocatable :: crd0
real(rp), dimension(:), allocatable :: sol
real(rp), dimension(:), allocatable :: nitsol_rwork
real(rp), dimension(:), pointer :: pvcrd => null()
real(rp), dimension(:), pointer :: pvfrc => null()

contains
 
!*******************************************************************************
 
subroutine bds_init(cpar, num_atoms, job_tag, ierr)
    !! Initializes the BD solver.

    type(ctrlpar_t),  intent(in) :: cpar
    integer, intent(in) :: num_atoms
    character(len=*), intent(in) :: job_tag
    integer,         intent(out) :: ierr
    integer :: lrwrk

    ierr = 0

    lhdia = cpar%lhdia;       write_traj = cpar%write_traj
    nts_sim = cpar%nts_sim;   nts_dump = cpar%nts_dump
    nts_samp = cpar%nts_samp; nts_log = cpar%nts_log
    mob_fctr = cpar%mob_fctr; bdintg = cpar%bdintg
    nts_mobsam = cpar%nts_mobsam
    tim_stp = cpar%tim_stp

    if (mob_fctr == 'LANC') then
        lanc_mxitr = cpar%lanc_mxitr
        lanc_tol = cpar%lanc_tol
    end if

    if (bdintg == 'SE') then
        nlmxitr = cpar%se_nlmxitr; kdmax = cpar%se_kdmax
        ftol = cpar%se_tol(1); stptol = cpar%se_tol(2)
    end if

    fn_revive = cpar%fn_revive//job_tag
    sqrt_two_dt = sqrt(2.0_rp*tim_stp)

    !Allocate memory for BD drift-diffusion SDE
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

    if (bdintg == 'SE') then
        allocate(crd0(3*num_atoms))
        allocate(sol(3*num_atoms))
        lrwrk = 3*num_atoms*(kdmax+5)+kdmax*(kdmax+3)
        allocate(nitsol_rwork(lrwrk))
    end if

    end subroutine

!******************************************************************************

subroutine bds_finish()
    !! Clears up memory allocated in [[bds_init]].

    call brn_finish()

    if (allocated(mob)) deallocate(mob)
    if (allocated(drift)) deallocate(drift)
    if (allocated(diffusion)) deallocate(diffusion)

    if (allocated(crd0)) deallocate(crd0)
    if (allocated(sol)) deallocate(sol)
    if (allocated(nitsol_rwork)) deallocate(nitsol_rwork)

    pvcrd => null(); pvfrc => null()
    psmbx => null(); patc => null()

    end subroutine

!******************************************************************************

subroutine bds_run(nts_beg, simbox, atc, ierr)
    !! Driver for BD integrator.
    !!
    !! Repeatedly calls [[integrate_em]] or [[integrate_se]] to update atom
    !! positions.

    integer(ip_long), intent(in) :: nts_beg
    type(smbx_t), intent(in), target :: simbox
    type(atmcfg_t), intent(in out), target :: atc
    integer, intent(out) :: ierr
        !! Error flag
    integer(ip_long) :: nts
    real(rp), dimension(3) :: com
    type(c_ptr) :: loc
    logical :: limud

    nts = nts_beg; ierr = 0

    !Assign pointers
    psmbx => simbox; patc => atc
    loc = c_loc(patc%coordinates)
    call c_f_pointer(loc, pvcrd, [3*patc%num_atoms])
    loc = c_loc(patc%forces)
    call c_f_pointer(loc, pvfrc, [3*patc%num_atoms])

    !Is this an isolated untethered molecule in an unbounded domain?
    if ( (patc%num_tethers == 0) .and. (psmbx%imcon == 0) ) then
        limud = .true.
    else
        limud = .false.
    end if

    !For isolated untethered molecule, ensure c.o.m. is at the origin.
    if (limud) call psmbx%to_origin(atc%coordinates)

    !Log & dump starting configuration
    call logger%log_msg('nts = '//str_from_num(nts))
    call write_dump(nts, psmbx, patc, fn_revive)

    do 
        if (nts > nts_sim) exit
        if (bdintg=='EM') then
            call integrate_em(ierr)
        else if (bdintg=='SE') then
            call integrate_se(ierr)
        end if
        if (ierr /= 0) return

        nts = nts + 1

        !Apply boundary conditions: For isolated untethered molecule,
        !update molecule_com & bring c.o.m. to the center of the box.
        if ( limud ) then
            call psmbx%to_origin(patc%coordinates, com)
            patc%molecule_com = patc%molecule_com + com
        end if

        !For PBC, wrap atom positions
        if (psmbx%imcon /= 0) call psmbx%wrap_all(patc%coordinates)

        !Accumulate statistics
        call stats_accumulate(psmbx, patc)

        !Logging
        if (mod(nts,nts_log) == 0) then
            call logger%log_msg('nts = '//str_from_num(nts))
        end if

        !Dump revive file
        if (mod(nts,nts_dump)==0) then
            call write_dump(nts, psmbx, patc, fn_revive)
        end if

        if (mod(nts,nts_samp)==0) then
            !Write stats
            call stats_write(nts, psmbx, atc)
            !Write traj
            if (write_traj) call traj%append_frame(nts, patc%coordinates)
        end if
    end do

    !Dump the final configuration
    call write_dump(nts, psmbx, patc, fn_revive)

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
    real(rp), dimension(0) :: rpar !Dummy arg, 0-dim array
    integer, dimension(0)  :: ipar !Dummy arg, 0-dim array
    integer, dimension(10) :: input
    integer, dimension(6) :: info
    integer :: n
    character(len=256) :: msg
    !integer :: iplvl, ipunit
    !common /nitprint/ iplvl, ipunit
    external nitsol

    ierr = 0; n = 3*patc%num_atoms
    rpar = 0.0_rp; ipar = 0
    input = 0
    input(1) = nlmxitr; input(4) = kdmax; input(5) = 0
    input(7) = 1
    !iplvl = 0; ipunit = logger%fu

    !Calculate the diffusion term
    call calc_diffusion(ierr)
    if (ierr /= 0) return

    !Evaluate the nonlinear function arising out of discretizing the SDE in
    !semi-implicit form
    crd0 = pvcrd !Coordinates at the end of last time step
    sol = pvcrd  !Initial guess

    call nitsol(n, sol, se_fval, se_jacv, ftol, stptol, input, &
        info, nitsol_rwork, rpar, ipar, ierr, ddot, dnrm2)

    msg = ''
    select case (ierr)
    case(:-1)
        write(msg,'(a,i0,a)') 'illegal value in input(', abs(ierr), ')'
    case(1)
        write(msg,'(a)') 'max nonlinear iterations reached'
    case(2)
        write(msg,'(a)') 'failed to evaluate nonlinear function'
    case(3)
        write(msg,'(a)') 'failed to evaluate J*v'
    case(4)
        write(msg,'(a)') 'failed to evaluate P^{-1}*v'
    case(5)
        write(msg,'(a)') 'insufficient initial norm reduction'
    case(6)
        write(msg,'(a)') 'failed to reach an acceptable step by backtracking'
    end select
    if (len_trim(msg) > 0) then
        call logger%log_msg(msg)
        return
    end if
    !Comment out for verbose output from nitsol.
    !write(msg, '(4(a,i0))') 'nfe = ', info(1), ', nli = ', info(4), &
    !    ', nni = ', info(5), ', nbktrk = ', info(6)
    !call logger%log_msg(msg)

    pvcrd = sol

    if (lhdia) then
        if (cntr_mobsam == nts_mobsam) cntr_mobsam = 0
    end if

    end subroutine

!******************************************************************************

subroutine se_fval(n, xcur, fcur, rpar, ipar, ierr)
    !! Calculates the nonlinear function.

    integer, intent(in) :: n
    real(rp), dimension(*), intent(in) :: xcur
    real(rp), dimension(*), intent(out) :: fcur
    real(rp), dimension(*), intent(in) :: rpar
    integer, dimension(*), intent(in) :: ipar
    integer, intent(out) :: ierr

    ierr = 0; fcur(1:n) = 0.0_rp; pvcrd = xcur(1:n)

    call calc_drift(ierr)
    if (ierr /= 0) return

    if (lhdia) then
        fcur(1:n) = pvcrd - drift - crd0 - diffusion(:,cntr_mobsam)
    else
        fcur(1:n) = pvcrd - drift - crd0 - diffusion(:,1)
    end if

    end subroutine

!******************************************************************************

subroutine se_jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, ierr)
    !! Calculates jacobian times vector product. This is a dummy subroutine.

    integer, intent(in) :: n
    real(rp), dimension(:), intent(in) :: xcur
    real(rp), dimension(:), intent(in) :: fcur
    integer, intent(in) :: ijob
    real(rp), dimension(:), intent(in) :: v
    real(rp), dimension(:), intent(out) :: z
    real(rp), dimension(:), intent(in) :: rpar
    integer, dimension(:), intent(in) :: ipar
    integer, intent(out) :: ierr

    ierr = 0

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
        call calc_rpy_tensor(patc%coordinates, mob)
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
            do j = 1, size(mob,2)
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
    call ia_calc_forces(psmbx, patc, ierr)
    if (ierr /= 0) return

    !Calculate ambient flow velocity at the atom locations.
    select case(patc%flow_style)
    case(1)
        !Steady simple shear: Flow along x, gradient along y
        do i = 1, patc%num_atoms
            drift(3*i-2) = patc%flow_params(1)*patc%coordinates(2,i)
        end do
    case(2)
        !Steady planar extension in x-y plane
        do i = 1, patc%num_atoms
            drift(3*i-2) =  patc%flow_params(1)*patc%coordinates(1,i)
            drift(3*i-1) = -patc%flow_params(1)*patc%coordinates(2,i)
        end do
    case(3)
        !Steady uniaxial extension along x
        do i = 1, patc%num_atoms
            drift(3*i-2) =  patc%flow_params(1)*patc%coordinates(1,i)
            drift(3*i-1) = -0.5_rp*patc%flow_params(1)*patc%coordinates(2,i)
            drift(3*i)   = -0.5_rp*patc%flow_params(1)*patc%coordinates(3,i)
        end do
    case default
        continue
    end select

    !drift <- mob * forces + drift. Only the lower triangular part of
    !mob is accessed.
    if (lhdia) then
        call dsymv('L', 3*patc%num_atoms, 1.0_rp, mob, 3*patc%num_atoms, &
                    pvfrc, 1, 1.0_rp, drift, 1)
    else
        drift = drift + pvfrc
    end if

    drift = tim_stp * drift

    end subroutine

!******************************************************************************

subroutine calc_rpy_tensor(coordinates, mob)
    !! Calculates the RPY approximation to the mobility tensor.

    real(rp), dimension(:,:), intent(in) :: coordinates
        !! (3, num_atoms*) matrix; stores the atom positions.
    real(rp), dimension(:,:), intent(out) :: mob
        !! (3*num_atoms*, 3*num_atoms*) matrix; stores the mobility tensor.
    real(rp), dimension(3)   :: ri, rj, rij
    real(rp), dimension(3,3) :: matrpy
    real(rp) :: rijm, irijm, irijm2
    real(rp) :: c1, c2, consij
    integer :: i, j, num_atoms

    num_atoms = size(coordinates,2)
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

end module bd_solver_m
