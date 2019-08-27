module m_mpcd
    !! Routines implementing the MPCD solver.

use f95_precision
use lapack95

use m_precision
use m_constants_math
use m_utils_math
use m_strings
use m_ran_num
use m_globals
use m_cell_list
use m_interaction, only: ia_calc_forces
use m_stats_io
use m_config_io
use m_logger, only: logger => master_logger

implicit none

private

public :: mpcd_init, mpcd_run, mpcd_finish

real(rp), dimension(:), allocatable :: mass
    !! *(num_atoms_tot,)* array. Stores mass of each atom (including MPCD
    !! atoms). The mass of each MPCD atom is taken as unity.

real(rp), dimension(:), allocatable :: buf_aic
    !! Buffer for atoms in cell storing mass, coordinates, & velocities.
    !! For *nc* atoms, the first *nc* elements stores the mass, the next
    !! *3\*nc* elements stores the coordinates, and the next *3\*nc* elements
    !! store the velocities. The velocities part of this buffer may be overwritten
    !! for relative velocity calculations.

contains
 
!*******************************************************************************
 
subroutine mpcd_init(ierr)
    !! Initializes the MPCD solver. 

    integer, intent(out) :: ierr
    integer :: i, at_i
    logical :: lexists

    ierr = 0

    !Allocate memory for velocities
    if ( .not. allocated(velocities) ) then
        allocate( velocities(3,num_atoms_tot) )
        velocities = 0.0_rp
    end if

    !Mass of each atom
    if ( .not. allocated(mass) ) then
        allocate( mass(num_atoms_tot) )
        mass = 1.0_rp
        do i = 1, num_atoms
            at_i = atoms(i)
            mass(i) = atom_mass(at_i)
        end do
    end if

    !In case of new simulation, need to assign MPCD atom positions.
    if (.not. lrevive) then
        call simbox%get_rnd_points( coordinates(:,num_atoms+1:num_atoms_tot) )
    end if

    !Allocate buffer for atoms in a cell. Setting it to thrice the average
    !number of MPCD atoms per cell. Each atom has 7 quantities (mass,
    !coordinates, & velocities) associated with it. The buffer will be expanded
    !if necessary.
    allocate( buf_aic(7*(3*mpcd_avnc)) )
    buf_aic = 0.0_rp

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
                if (traj_wmpcd) then
                    !With MPCD data
                    call traj%create(fn_traj//trim(adjustl(job_tag)), &
                        num_atoms, num_mpcd_atoms, traj_frmcmp)
                else
                    !Without MPCD data
                    call traj%create(fn_traj//trim(adjustl(job_tag)), &
                        num_atoms, 0, traj_frmcmp)
                end if
            end if
        else
            !Create new trajectory file
            if (traj_wmpcd) then
                !With MPCD data
                call traj%create(fn_traj//trim(adjustl(job_tag)), &
                    num_atoms, num_mpcd_atoms, traj_frmcmp)
            else
                !Without MPCD data
                call traj%create(fn_traj//trim(adjustl(job_tag)), &
                    num_atoms, 0, traj_frmcmp)
            end if
        end if
    end if

    end subroutine

!******************************************************************************

subroutine mpcd_finish()
    !! Clean up MPCD solver.

    if (allocated(mass)) deallocate(mass)
    if (allocated(buf_aic)) deallocate(buf_aic)
    call traj%close()

    end subroutine

!******************************************************************************

subroutine mpcd_run()
    !! Driver for MPCD integrator.
    !!
    !! Repeatedly calls [[mpcd_stream]] and [[mpcd_collide]] to update atom
    !! positions.
    !!
    !! Handling of flow field is not imlemented yet. If a non-zero `flow_style`
    !! is specified, it will be set to zero.

    integer :: ierr

    ierr = 0
    if (flow_style /= 0) flow_style = 0

    !Log & dump starting configuration
    call logger%info('mpcd_run', ' nts: '//str_from_num(nts))
    call write_dump(fn_revive//trim(adjustl(job_tag)))

    !Equilibration run
    if (leql) then
        do while (nts < nts_eql)
            !Streaming
            call mpcd_stream(ierr)
            if (ierr /= 0) return

            !Collision
            call mpcd_collide()

            nts = nts + 1

            !Logging
            if (mod(nts,nts_log) == 0) then
                call logger%info('mpcd_run', ' nts: '//str_from_num(nts))
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
        !If continue on to production run, create new stats file.
        if (nts_sim > 0) call stats_init()
    end if

    !Production run
    do while (nts < nts_sim)
        call mpcd_stream(ierr)
        if (ierr /= 0) return

        call mpcd_collide()

        nts = nts + 1

        !Logging
        if (mod(nts,nts_log) == 0) then
            call logger%info('mpcd_run', ' nts: '//str_from_num(nts))
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

subroutine mpcd_stream(ierr)
    !! Performs one step of streaming.

    integer, intent(out) :: ierr
    real(rp) :: dt_md
    real(rp) :: half_dt
    integer :: i, istp

    dt_md = tim_stp/nts_md
    half_dt = 0.5_rp*dt_md
    istp = 0

    !Set cell size to collision cells
    if (use_cell_list) call cl_set_cell_size(rcutoff)

    !Setting these to zero for the case num_atoms = 0.
    energy_tot = 0.0_rp
    stress = 0.0_rp

    !Update non-MPCD atoms via MD steps
    if (num_atoms > 0) then
        do while (istp < nts_md)
            !Update velocities over a half-step
            do i = 1, num_atoms
                velocities(:,i) = velocities(:,i) + half_dt*forces(:,i)/mass(i)
            end do

            !Update position over full step
            do i = 1, num_atoms
                coordinates(:,i) = coordinates(:,i) + dt_md*velocities(:,i)
            end do

            !Wrap around periodic boundaries
            call simbox%wrap_all(coordinates(:,1:num_atoms))

            !Update forces over full step
            call ia_calc_forces(ierr)
            if (ierr /= 0) return

            !Update velocities over full step
            do i = 1, num_atoms
                velocities(:,i) = velocities(:,i) + half_dt*forces(:,i)/mass(i)
            end do

            istp = istp + 1
        end do
    end if

    !Update MPCD particles over a streaming step
    do i = num_atoms+1, num_atoms_tot
        coordinates(:,i) = coordinates(:,i) + tim_stp*velocities(:,i)
    end do

    !Wrap MPCD particles around periodic boundaries
    call simbox%wrap_all(coordinates(:,num_atoms+1:))

    end subroutine

!******************************************************************************

subroutine mpcd_collide()
    !! Performs one step of MPCD collision.

    real(rp), dimension(3) :: gsv
    real(rp), dimension(3) :: ri, vi, com, vcom, vrcom
    real(rp), dimension(3) :: angmom, angvel
    real(rp), dimension(3,3) :: mat_mi, inv_mat_mi
    real(rp), dimension(3,3) :: stress_kin
    integer, dimension(:), pointer :: aic => null()
    real(rp) :: mass_tot, massi, var
    integer :: num_cells, na_cell
    integer :: ofstc, ofstv
    integer :: i, iatm, icell, ibeg, iend

    !Grid shift vector
    call get_rv_uniform(-0.5_rp, 0.5_rp, gsv)

    !Forward shift
    do i = 1, num_atoms_tot
        coordinates(:,i) = coordinates(:,i) + gsv
    end do

    !Wrap around periodic boundaries
    call simbox%wrap_all(coordinates)

    !Set cell size to collision cells
    call cl_set_cell_size(1.0_rp)
    call cl_build(coordinates)

    num_cells = cl_get_num_cells()

    !Momentum exchange via collision, velocity update (MPC-AT+a version)

    !Loop over all collision cells
    do icell = 0, (num_cells-1)
        !Get cell contents
        call cl_get_contents(icell, aic)

        !Number of atoms in cell
        na_cell = size(aic)

        !For fewer than 2 atoms, no need for relative velocities.
        if (na_cell < 2) cycle

        !Check if cell velocity buffer needs to be expanded
        if ( size(buf_aic) < 7*na_cell ) then
            deallocate(buf_aic)
            allocate( buf_aic(7*na_cell) )
        end if
        
        !Populate buf_aic
        ofstc = na_cell; ofstv = na_cell + 3*na_cell
        do i = 1, size(aic)
            iatm = aic(i)
            buf_aic(i) = mass(iatm)

            ibeg = ofstc+3*i-2; iend = ofstc+3*i
            buf_aic(ibeg:iend) = coordinates(:,iatm)

            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            buf_aic(ibeg:iend) = velocities(:,iatm)
        end do

        !Total mass, center-of-mass, & c.o.m. velocity
        com = 0.0_rp; vcom = 0.0_rp; mass_tot = 0.0_rp
        do i = 1, na_cell
            mass_tot = mass_tot + buf_aic(i)

            ibeg = ofstc+3*i-2; iend = ofstc+3*i
            com = com + buf_aic(i)*buf_aic(ibeg:iend)

            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            vcom = vcom + buf_aic(i)*buf_aic(ibeg:iend)
        end do
        com = com/mass_tot
        vcom = vcom/mass_tot

        !Draw relative velocities & store in the buf_aic. This overwrites the
        !existing velocities.
        ibeg = ofstv+1; iend = ofstv+3*na_cell
        call get_rv_gaussian( 0.0_rp, 1.0_rp, buf_aic(ibeg:iend) )
        !Multiply by the variance, i.e. kT/m
        do i = 1, na_cell
            var = 1.0_rp/buf_aic(i)
            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            buf_aic(ibeg:iend) = var*buf_aic(ibeg:iend)
        end do

        !Calcuate the mean of the relative velocities.
        vrcom = 0.0_rp
        do i = 1, na_cell
            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            vrcom = vrcom + buf_aic(i)*buf_aic(ibeg:iend)
        end do
        vrcom = vrcom/mass_tot

        !Subtract off the c.o.m. position & vrcom
        do i = 1, na_cell
            ibeg = ofstc+3*i-2; iend = ofstc+3*i
            buf_aic(ibeg:iend) = buf_aic(ibeg:iend) - com

            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            buf_aic(ibeg:iend) = buf_aic(ibeg:iend) - vrcom
        end do

        !Add the relative velocities to the c.o.m. velocity. This produces new
        !velocities that do not conserve angular momentum, i.e. this is
        !MPC-AT-a.
        do i = 1, na_cell
            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            buf_aic(ibeg:iend) =  buf_aic(ibeg:iend) + vcom
        end do

        !For only two atoms in the cell, the above is equivalent to MPC-AT+a.
        !Just copy the velocities back and move on to the next cell.
        if (na_cell == 2) then
            !Explicitly unrolling the loop over i as i = 1 & i = 2.
            i = 1; iatm = aic(i); ibeg = ofstv+3*i-2; iend = ofstv+3*i
            velocities(:,iatm) = buf_aic(ibeg:iend)
            i = 2; iatm = aic(i); ibeg = ofstv+3*i-2; iend = ofstv+3*i
            velocities(:,iatm) = buf_aic(ibeg:iend)
            cycle
        end if

        !Calculate moment of intertia matrix (symmetric) about the c.o.m. and angular
        !momentum.
        mat_mi = 0.0_rp; angmom = 0.0_rp
        do i = 1, na_cell
            massi = buf_aic(i)
            ri = buf_aic(ofstc+3*i-2:ofstc+3*i)
            mat_mi(1,1) = mat_mi(1,1) + massi*( ri(2)**2 + ri(3)**2 )
            mat_mi(2,1) = mat_mi(2,1) - massi*( ri(2)*ri(1) )
            mat_mi(3,1) = mat_mi(3,1) - massi*( ri(3)*ri(1) )
            mat_mi(2,2) = mat_mi(2,2) + massi*( ri(1)**2 + ri(3)**2 )
            mat_mi(3,2) = mat_mi(3,2) - massi*( ri(3)*ri(2) )
            mat_mi(3,3) = mat_mi(3,3) + massi*( ri(1)**2 + ri(2)**2 )

            vi = buf_aic(ofstv+3*i-2:ofstv+3*i)
            call cross(ri, vi, vrcom) !Using vrcom as a buffer here
            angmom = angmom + massi*vrcom
        end do
        !Update the symmetric elements
        mat_mi(1,2) = mat_mi(2,1)
        mat_mi(1,3) = mat_mi(3,1)
        mat_mi(2,3) = mat_mi(3,2)

        !Solve for angular velocity. If the mat_mi is singular, do a least
        !squares solutions. Else, get the explicit inverse and multiply.
        !Checking determinant against a small value 1E-8.
        if (det(mat_mi) > 1.0E-10_rp) then
            !Invert the moment of inertia matrix
            call invert_mat33(mat_mi, inv_mat_mi)
            !Solve for angular velocity
            angvel = matmul(inv_mat_mi, -angmom)
        else
            !TODO: The following must be replaced with a more efficient SVD
            !routine for 3 x 3 matrices. As of now, calling LAPACK for
            !least-squares solution.
            angmom = -angmom
            call gelsy(mat_mi, angmom)
            angvel = angmom
        end if

        !Correct velocities to enforce angular mometum conservation
        do i = 1, na_cell
            iatm = aic(i)
            massi = buf_aic(i)
            ri = buf_aic(ofstc+3*i-2:ofstc+3*i)
            call cross(angvel, ri, vi) !Using vi as a buffer here
             
            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            buf_aic(ibeg:iend) = buf_aic(ibeg:iend) + vi
        end do

        !Copy back
        do i = 1, na_cell
            iatm = aic(i)
            ibeg = ofstv+3*i-2; iend = ofstv+3*i
            velocities(:,iatm) = buf_aic(ibeg:iend)
        end do
    end do

    !Reverse shift
    do i = 1, num_atoms_tot
        coordinates(:,i) = coordinates(:,i) - gsv
    end do

    !Wrap around periodic boundaries
    call simbox%wrap_all(coordinates)

    !Update energy and stress
    energy_kin = 0.0_rp
    stress_kin = 0.0_rp
    stress_slvnt = 0.0_rp

    mass_tot = 0.0_rp; com = 0.0_rp; vcom = 0.0_rp
    do i = 1, num_atoms_tot
        massi = mass(i)
        mass_tot = mass_tot + massi
        com = com + massi*coordinates(:,i)
        vcom = vcom + massi*velocities(:,i)
    end do
    com = com/mass_tot
    vcom = vcom/mass_tot

    do i = 1, num_atoms
        massi = mass(i)
        energy_kin = energy_kin + 0.5_rp*massi*sum(velocities(:,i)**2)

        vi = velocities(:,i) - vcom
        stress_kin(1,1) = stress_kin(1,1) + massi*vi(1)*vi(1)
        stress_kin(2,1) = stress_kin(2,1) + massi*vi(2)*vi(1)
        stress_kin(3,1) = stress_kin(3,1) + massi*vi(3)*vi(1)

        stress_kin(1,2) = stress_kin(1,2) + massi*vi(1)*vi(2)
        stress_kin(2,2) = stress_kin(2,2) + massi*vi(2)*vi(2)
        stress_kin(3,2) = stress_kin(3,2) + massi*vi(3)*vi(2)

        stress_kin(1,3) = stress_kin(1,3) + massi*vi(1)*vi(3)
        stress_kin(2,3) = stress_kin(2,3) + massi*vi(2)*vi(3)
        stress_kin(3,3) = stress_kin(3,3) + massi*vi(3)*vi(3)
    end do

    do i = num_atoms+1, num_atoms_tot
        !Mass of MPCD particles is 1 (non-dimensionalized form).
        energy_kin = energy_kin + 0.5_rp*sum(velocities(:,i)**2)

        vi = velocities(:,i) - vcom
        stress_slvnt(1,1) = stress_slvnt(1,1) + vi(1)*vi(1)
        stress_slvnt(2,1) = stress_slvnt(2,1) + vi(2)*vi(1)
        stress_slvnt(3,1) = stress_slvnt(3,1) + vi(3)*vi(1)

        stress_slvnt(1,2) = stress_slvnt(1,2) + vi(1)*vi(2)
        stress_slvnt(2,2) = stress_slvnt(2,2) + vi(2)*vi(2)
        stress_slvnt(3,2) = stress_slvnt(3,2) + vi(3)*vi(2)

        stress_slvnt(1,3) = stress_slvnt(1,3) + vi(1)*vi(3)
        stress_slvnt(2,3) = stress_slvnt(2,3) + vi(2)*vi(3)
        stress_slvnt(3,3) = stress_slvnt(3,3) + vi(3)*vi(3)
    end do
    energy_tot = energy_tot + energy_kin

    stress_kin = stress_kin/simbox%volume
    stress = stress + stress_kin
    stress_slvnt = stress_slvnt/simbox%volume

    end subroutine

!******************************************************************************

end module m_mpcd
