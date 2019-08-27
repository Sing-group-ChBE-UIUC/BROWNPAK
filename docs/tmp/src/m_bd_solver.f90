module m_bd_solver
    !! Routines implementing Brownian Dynamics (BD) solver.

use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

use blas95
use lapack95
use f95_precision

use m_precision
use m_constants_math
use m_strings
use m_ran_num
use m_globals
use m_interaction, only: ia_calc_forces
use m_stats_io
use m_config_io
use m_logger, only: logger => master_logger

implicit none

private

public :: bds_init, bds_run, bds_finish

real(rp), dimension(:),   allocatable :: drift
real(rp), dimension(:),   allocatable :: diffusion
real(rp), dimension(:,:), allocatable :: mobility

contains
 
!*******************************************************************************
 
subroutine bds_init(ierr)
    !! Initializes the BD solver.

    integer, intent(out) :: ierr
    logical :: lexists

    ierr = 0

    !Allocate memory for BD integration
    allocate(drift (3*num_atoms))
    drift = 0.0_rp

    allocate(diffusion(3*num_atoms))
    diffusion = 0.0_rp

    if (lhdia) then
        allocate(mobility (3*num_atoms,3*num_atoms))
        mobility = 0.0_rp
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
                call traj%create(fn_traj//trim(adjustl(job_tag)), num_atoms, 0, traj_frmcmp)
            end if
        else
            !Create new trajectory file
            call traj%create(fn_traj//trim(adjustl(job_tag)), num_atoms, 0, traj_frmcmp)
        end if
    end if

    end subroutine

!******************************************************************************

subroutine bds_finish()
    !! Clears up memory allocated in [[bds_init]].

    call traj%close()
    if (allocated(drift)) deallocate(drift)
    if (allocated(diffusion)) deallocate(diffusion)
    if (allocated(mobility)) deallocate(mobility)

    end subroutine

!******************************************************************************

subroutine bds_run()
    !! Driver for BD integrator.
    !!
    !! Repeatedly calls [[bds_integrate_fd]] or [[bds_integrate_hi]] to update
    !! atom positions.

    real(rp), dimension(3) :: com
    integer :: flow_style_
    integer :: ierr

    ierr = 0

    !For isolated untethered molecule, ensure c.o.m. is at the center of the
    !box.
    if ( (num_tethers == 0) .and. (imcon == 0) ) then
        call simbox%to_center(coordinates)
    end if

    !Log & dump starting configuration
    call logger%info('bds_run', ' nts: '//str_from_num(nts))
    call write_dump(fn_revive//trim(adjustl(job_tag)))

    !Equilibration run. This is done under free-draining condition.
    if (leql) then
        !Turn off flow during equilibration
        flow_style_ = flow_style
        flow_style = 0

        do while (nts < nts_eql)
            call bds_integrate_fd(ierr)
            if (ierr /= 0) return
            nts = nts + 1

            !Apply boundary conditions: For isolated untethered molecule,
            !bring c.o.m. to the origin. No need to update molecule_com during
            !equilibration.
            if ( (num_tethers == 0) .and. (imcon == 0) ) then
                call simbox%to_center(coordinates)
            end if

            !For PBC, wrap atom positions
            if (imcon /= 0) call simbox%wrap_all(coordinates)

            !Logging
            if (mod(nts,nts_log) == 0) then
                call logger%info('bds_run', ' nts: '//str_from_num(nts))
            end if

            !Dump revive file
            if (mod(nts,nts_dump)==0) then
                call write_dump(fn_revive//trim(adjustl(job_tag)))
            end if

            !Equilibration stats
            if (write_eql_stats) then
                if (mod(nts,nts_eql_samp)==0) call stats_write()
            end if
        end do
        call logger%info('finish_equil_run', 'Equilibration completed')
        call stats_finish()
        call write_dump(fn_revive//trim(adjustl(job_tag)))
        leql = .false.
        nts = 0
        !Turn flow back on
        flow_style = flow_style_
        !If continue on to production run, create new stats file.
        if (nts_sim > 0) call stats_init()
    end if

    !Production run.
    do while (nts < nts_sim)
        if (lhdia) then
            call bds_integrate_hi(ierr)
        else
            call bds_integrate_fd(ierr)
        end if
        if (ierr /= 0) return
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
            call logger%info('bds_run', ' nts: '//str_from_num(nts))
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
        end if
    end do
    !Dump the final configuration
    call write_dump(fn_revive//trim(adjustl(job_tag)))

    end subroutine

!******************************************************************************

subroutine bds_integrate_fd(ierr)
    !! Performs one step of free-draining BD integeration.

    integer, intent(out) :: ierr
    real(rp), dimension(:), pointer :: fptr_forces
    type(c_ptr) :: cptr_forces
    real(rp) :: sqrt_two_dt
    integer :: nb
    integer :: i

    ierr = 0
    drift = 0.0_rp; diffusion = 0.0_rp
    sqrt_two_dt = sqrt(2*tim_stp)
    nb = size(coordinates,2)

    !Update forces
    call ia_calc_forces(ierr)
    if (ierr /= 0) return

    !Calculate ambient flow velocity at the atom locations and store it in
    !module variable `drift`.
    call calc_ambient_velocity()

    !Calculate mobility matrix, take its dot product with force, & add to drift.
    !Using MKL routine 'symv' here.

    cptr_forces = c_loc(forces)
    call c_f_pointer(cptr_forces, fptr_forces, [3*nb])

    !No hydrodynamic interaction 
    drift = drift + fptr_forces

    !Calculate random force & store it in module variable diffusion.
    call get_rv_gaussian(0.0_rp, 1.0_rp, diffusion, 1000000) 

    !Update positions
    do i = 1, nb
        coordinates(:,i) = coordinates(:,i) + tim_stp*drift(3*i-2:3*i) &
           + sqrt_two_dt*diffusion(3*i-2:3*i)
    end do

    end subroutine

!******************************************************************************

subroutine bds_integrate_hi(ierr)
    !! Performs one step of BD integeration including HI.

    integer, intent(out) :: ierr
    real(rp), dimension(:), pointer :: fptr_forces
    type(c_ptr) :: cptr_forces
    real(rp) :: sqrt_two_dt
    integer :: nb, info
    integer :: i

    ierr = 0
    drift = 0.0_rp; diffusion = 0.0_rp
    sqrt_two_dt = sqrt(2*tim_stp)
    nb = size(coordinates,2)

    !Update forces
    call ia_calc_forces(ierr)
    if (ierr /= 0) return

    !Calculate ambient flow velocity at the atom locations and store it in
    !module variable `drift`.
    call calc_ambient_velocity()

    !Calculate mobility matrix, take its dot product with force, & add to drift.
    !Using MKL routine 'symv' here.

    cptr_forces = c_loc(forces)
    call c_f_pointer(cptr_forces, fptr_forces, [3*nb])

    !With hydrodynamic interaction 
    call calc_rpy_tensor()

    if (flow_style == 0) then
        call symv(mobility, fptr_forces, drift, 'U')
    else
        call symv(mobility, fptr_forces, drift, 'U', 1.0_rp, 1.0_rp)
    end if

    !Calculate random force & store it in module variable diffusion.
    call get_rv_gaussian(0.0_rp, 1.0_rp, diffusion, 1000000) 

    !Performs the operation: x = A.x, where A is the square root of mobility.
    !Only the upper triangular part of the square root matrix is
    !considered. For hi_method='cholesky', the square root matrix is calculated
    !before trmv.
    !Trmv operation on the sqrt of mobility matrix.

    if (mob_fctr == 'CHOL') then
        call potrf(mobility, 'U', info)
        if (info /= 0) then
            ierr = 1
            call logger%fatal( 'potrf:', 'err code '//str_from_num(info) )
            return
        end if
        ! Default for trmv: access the upper triangular part only
        call trmv(mobility, diffusion, 'U', 'N', 'U')
    else if (mob_fctr == 'KRYL') then
        call calc_Bdw_kryl (mobility, diffusion, ierr)
        if (ierr /= 0) return
    end if

    !Update positions
    do i = 1, nb
        coordinates(:,i) = coordinates(:,i) + tim_stp*drift(3*i-2:3*i) &
            + sqrt_two_dt*diffusion(3*i-2:3*i)
    end do

    end subroutine

!******************************************************************************

subroutine calc_ambient_velocity()
    !! Calculates ambient velocity at the atom positions and store it in module
    !! variable `drift`.

    integer :: i, nb

    nb = size(coordinates,2)

    select case(flow_style)
    case(1)
        !Steady simple shear: Flow along x, gradient along y
        do i = 1, nb
            drift(3*i-2) = flow_params(1)*coordinates(2,i)
        end do
    case(2)
        !Steady planar extension in x-y plane
        do i = 1, nb
            drift(3*i-2) =  flow_params(1)*coordinates(1,i)
            drift(3*i-1) = -flow_params(1)*coordinates(2,i)
        end do
    case(3)
        !Steady uniaxial extension along x
        do i = 1, nb
            drift(3*i-2) =  flow_params(1)*coordinates(1,i)
            drift(3*i-1) = -0.5_rp*flow_params(1)*coordinates(2,i)
            drift(3*i)   = -0.5_rp*flow_params(1)*coordinates(3,i)
        end do
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine calc_rpy_tensor()
    !! Calculates the RPY approximation to the mobility tensor. Overwrites
    !! the upper triangular part of `mobility`.

    real(rp), dimension(3)   :: ri, rj, rij
    real(rp), dimension(3,3) :: matrpy
    real(rp) :: rijm, irijm, irijm2
    real(rp) :: c1, c2, consij
    integer :: nb
    integer :: i, j

    nb = size(coordinates,2)

    !Reset mobility matrix to zero.
    mobility = 0.0_rp

    !Calculate the RPY tensor (in strictly upper triangular form)
    do j = 2, nb
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

            mobility(3*i-2:3*i,3*j-2:3*j) = matrpy
        end do
    end do
   
    !Copy upper triangular part to lower triangular part. This is not necessary
    !if the lower triangular part is not accessed in subsequent calculations.

    !do j = 2, (3*nb-1)
    !    do i = 1, j-1
    !        rpy_tensor(j,i) = rpy_tensor(i,j)
    !    end do
    !end do

    !Put one on the diagonal.
    do j = 1,3*nb
       mobility(j,j) = 1.0_rp
    end do

    end subroutine

!******************************************************************************

subroutine calc_Bdw_kryl (D, Bdw, ierr)
    !! Calculates *B.dW* using Krylov subspace method.

    real(rp), dimension(:,:), intent(in)  :: D
        !! (3N,3N) symmetric positive definite matrix.
        !!
        !! Diffusivity matrix. Only the upper triangular part of the matrix is used.
    real(rp), dimension(:),   intent(in out)  :: Bdw
        !! (3N,) vector. On entry contains vector **dW**. On return, contains 
        !! **B.dW**, where **B** is the square root matrix of **D**.
    integer,                  intent(out) :: ierr
        !! Error flag
    integer, parameter     :: dimk = 12
        ! Dimension of Krylov subspace
    real(rp), parameter    :: Ek_thres = 0.01_rp
    real(rp), dimension(size(D,2))       :: Bdwold
    real(rp), dimension(size(D,2))       :: w
    real(rp), dimension(size(D,2), dimk) :: V
    real(rp), dimension(dimk, dimk)      :: P
    real(rp), dimension(dimk)            :: diagH     
    real(rp), dimension(dimk)            :: diagHsav  
    real(rp), dimension(dimk)            :: subH      
    real(rp), dimension(dimk)            :: subHsav   
    real(rp), dimension(dimk)            :: sqrtHcol1 
    real(rp), dimension(dimk)            :: Ptranscol1
    real(rp) :: normz
    real(rp) :: Ek
    integer  :: ncol
    integer  :: j
    integer  :: info
    
    ierr = 0
    
    ncol = size(D,2)
    diagH  = 0.0_rp
    subH   = 0.0_rp
    Bdwold = 0.0_rp
    
    normz = norm2(Bdw)
    
    V(:,1) = Bdw/normz
    
    do j = 1, dimk
       call symv (D, V(:,j), w) !Default using upper triangular part
    
       if ( j > 1 ) then
           w = w - subH(j-1)*V(:,j-1)
       end if
    
       diagH(j) = dot_product(w, V(:,j))
    
       if ( j < dimk) then
           w = w - diagH(j)*V(:,j)
           subH(j) = norm2(w)
           V(:,j+1) = w/subH(j)
       end if
    
       !Update Bdw
    
       if ( j >= 5 ) then
            diagHsav(1:j) = diagH(1:j)
            subHsav(1:j) = subH(1:j)
    
            call stevd (diagHsav(1:j), subHsav(1:j), P(1:j,1:j), info)
      
            if ( info /= 0 ) then
                write(*,*) 'calc_Bdw_kryl: info = ', info
                ierr = 1
                return
            end if
    
            diagHsav(1:j) = sqrt( diagHsav(1:j) )
            
            Ptranscol1(1:j) = P(1,1:j)
      
            Ptranscol1(1:j) = diagHsav(1:j)*Ptranscol1(1:j)
      
            call gemv ( P(1:j,1:j), Ptranscol1(1:j), sqrtHcol1(1:j) )
      
            call gemv ( V(:,1:j), sqrtHcol1(1:j), Bdw, normz, 0.0_rp, 'N')
    
            !Error calculation
    
            if ( j <= 5) then !Check only after several iterations, say 5
                 Ek = 1.0_rp
            else
                 Ek = norm2(Bdw - Bdwold)/norm2(Bdwold)
            end if
    
            if ( Ek <= Ek_thres ) then
                 exit
            end if
    
            Bdwold = Bdw
    
       end if
    
    end do
    
    end subroutine

!******************************************************************************

end module m_bd_solver
