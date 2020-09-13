module m_stats_io
    !! Computes and writes properties calculated during simulation.

use m_precision
use m_constants_math
use m_strings
use m_globals

implicit none

private

public :: stats_init, stats_finish, stats_write

integer :: fu_stats
    !! Unit number of file `fn_stats`.
real(rp), dimension(3) :: span
    !! Span of a chain molecule. For multiple molecules this is averaged.
real(rp) :: reedsq
    !! Mean squared end-to-end distance of a linear chain molecule.
    !! For multiple molecules this is averaged. Not defined for rings.
real(rp) :: rgsq
    !! Mean squared gyration radius of a chain/ring molecule.
    !! For multiple molecules this is averaged.
real(rp), dimension(3) :: reev
    !! End-to-end vector. Not defined for rings.
real(rp) :: asph
    !! Asphericity of the entire molecule
real(rp) :: prol
    !! Prolateness of the entire molecule
real(rp) :: rgsq_bbone
    !! Mean squared gyration radius of the backbone of a branched
    !! chain/ring molecule.
real(rp) :: rgsq_sc
    !! Mean squared gyration radius of the side chains of a branched
    !! chain/ring molecule.
real(rp) :: reedsq_sc
    !! Mean squared end-to-end distance of the side chains of a branched
    !! chain/ring molecule.
real(rp) :: asph_sc
    !! Mean asphericity of the side chains of a branched chain/ring molecule.
real(rp) :: prol_sc
    !! Mean prolateness of the side chains of a branched chain/ring molecule.
real(rp), dimension(:,:), allocatable :: molbuf
    !! Buffer for the largest molecule in the system

contains

!******************************************************************************

subroutine stats_init()
    !! Set up for stats collection

    logical :: lexists
    integer :: n

    if (lrevive) then
        !Restarting simulation
        !Check if file exists
        inquire(file=fn_stats//trim(adjustl(job_tag)), exist=lexists)
        if (lexists) then
            !Open existing file for appending statistics
            open(newunit=fu_stats, file=fn_stats//trim(adjustl(job_tag)), &
                action='write', position='append', status='old')
        else
            !Open new file for writing statistics
            open(newunit=fu_stats, file=fn_stats//trim(adjustl(job_tag)), &
                action='write', status='new')
            !Write header
            call stats_write_hdr()
        end if
    else
        !New simulation
        !Open file for writing production (or relaxation) statistics
        open(newunit=fu_stats, file=fn_stats//trim(adjustl(job_tag)), &
            action='write', status='replace')
        !Write header
        call stats_write_hdr()
    end if

    !Allocate space for unwrapping molecule under PBC
    if ( (imcon /= 0) .and. (num_molecules > 0) ) then
        n = maxval(molecules(2,:))
        if (n > 1) allocate( molbuf(3,n) )
    end if

    end subroutine

!******************************************************************************

subroutine stats_finish()
    !! Closes any files opened in `stats_init`.

    integer :: fu
    logical :: lopened

    inquire(file=fn_stats//trim(adjustl(job_tag)), number=fu, &
        opened=lopened)
    if (lopened) close(fu)

    if (allocated(molbuf)) deallocate(molbuf)

    end subroutine

!******************************************************************************

subroutine stats_write_hdr()
    !! Driver for writing header of file `fn_stats`. This is called only after
    !! appropriate files have been opened.

    write(fu_stats,'(a12)', advance='no') 'nts'

    if (num_bonds > 0) then
        write(fu_stats,'(*(2x,a16))', advance='no') &
            'bndlen_avg', 'bndlen_min', 'bndlen_max', 'energy_bond'
    end if

    if (num_angles > 0) then
        write(fu_stats,'(2x,a16)', advance='no') 'energy_angle'
    end if

    if (num_dihedrals > 0) then
        write(fu_stats,'(2x,a16)', advance='no') 'energy_dihedral'
    end if

    if (lvdw) then
        write(fu_stats,'(2x,a16)', advance='no') 'energy_vdw'
    end if

    if (num_tethers > 0) then
        write(fu_stats,'(2x,a16)', advance='no') 'energy_tether'
    end if

    if (num_externals > 0) then
        write(fu_stats,'(2x,a16)', advance='no') 'energy_external'
    end if

    if (sim_style == 2) then
        write(fu_stats,'(2x,a16)', advance='no') 'energy_kin'
    end if

    write(fu_stats,'(2x,a16)', advance='no') 'energy_tot'

    !No further properties will be calculated for structure relaxation
    if (sim_style == 0) then
        write(fu_stats,*)
        return
    end if

    !Stress (written only during production run)
    if (num_atoms > 0) then
        write(fu_stats,'(*(2x,a16))', advance='no') 'sxx', 'syx', 'szx', &
            'sxy', 'syy', 'szy', 'sxz', 'syz', 'szz'
    end if
    if (sim_style == 2) then
        write(fu_stats,'(*(2x,a16))', advance='no') 'slv_sxx', 'slv_syx', 'slv_szx', &
            'slv_sxy', 'slv_syy', 'slv_szy', 'slv_sxz', 'slv_syz', 'slv_szz'
    end if

    !Structural properties: Size (for multiple molecules these will be
    !averaged quantities)
    if ( num_bonds > 0 ) then
        write(fu_stats,'(*(2x,a16))', advance='no') 'reedsq', 'rgsq',  &
            'spanx', 'spany', 'spanz'
    end if

    !Further structural properties (for single molecule)
    if ( (imcon == 0) .and. (num_bonds > 0) ) then
        write(fu_stats,'(*(2x,a16))', advance='no') 'reevx', 'reevy', &
            'reevz', 'asph', 'prol'
        if (num_branches > 0) then
            write(fu_stats,'(*(2x,a16))', advance='no') 'rgsq_bbone',  &
                'rgsq_sc', 'reedsq_sc', 'asph_sc', 'prol_sc'
        end if
        if (flow_style == 0) then
            write(fu_stats,'(*(2x,a16))', advance='no') 'comx', 'comy', 'comz'
        end if
    end if

    write(fu_stats,*)

    end subroutine

!******************************************************************************

subroutine stats_write()
    !! Writing statistics

    write(fu_stats,'(i12)', advance='no') nts

    if (num_bonds > 0) then
        write(fu_stats,'(*(2x,es16.7))', advance='no') &
            bndlen, bndlen_min, bndlen_max, energy_bond
    end if

    if (num_angles > 0) then
        write(fu_stats,'(2x,es16.7)', advance='no') energy_angle
    end if

    if (num_dihedrals > 0) then
        write(fu_stats,'(2x,es16.7)', advance='no') energy_dihedral
    end if

    if (lvdw) then
        write(fu_stats,'(2x,es16.7)', advance='no') energy_vdw
    end if

    if (num_tethers > 0) then
        write(fu_stats,'(2x,es16.7)', advance='no') energy_tether
    end if

    if (num_externals > 0) then
        write(fu_stats,'(2x,es16.7)', advance='no') energy_external
    end if

    if (sim_style == 2) then
        write(fu_stats,'(2x,es16.7)', advance='no') energy_kin
    end if

    write(fu_stats,'(2x,es16.7)', advance='no') energy_tot

    !No further properties will be calculated for structure relaxation
    if (sim_style == 0) then
        write(fu_stats,*)
        return
    end if

    !Stress (calculated only during production run)
    if (num_atoms > 0) then
        write(fu_stats,'(*(2x,es16.7))', advance='no') stress_accu
    end if
    if (sim_style == 2) then
        write(fu_stats,'(*(2x,es16.7))', advance='no') stress_slvnt
    end if

    !Structural properties: Size (for multiple molecules these will be
    !averaged quantities)
    if (imcon == 0) then
        call stats_compute_ic0()
    else
        call stats_compute_ic1()
    end if

    if ( num_bonds > 0 ) then
        write(fu_stats,'(*(2x,es16.7))', advance='no') reedsq, rgsq, span
    end if

    !Further structural properties (for single molecule in unbounded domain)
    if ( (imcon == 0) .and. (num_bonds > 0) ) then
        write(fu_stats,'(*(2x,es16.7))', advance='no') reev, asph, prol
        if (num_branches > 0) then
            write(fu_stats,'(*(2x,es16.7))', advance='no') rgsq_bbone,  &
                rgsq_sc, reedsq_sc, asph_sc, prol_sc
        end if
        if (flow_style == 0) then
            write(fu_stats,'(*(2x,es16.7))', advance='no') molecule_com
        end if
    end if

    write(fu_stats,*)

    end subroutine

!******************************************************************************

subroutine stats_compute_ic1()
    !! Computes statistics for a possibly multiple chains in a periodic domain.

    real(rp), dimension(3) :: com
    real(rp) :: rgsq_im
    integer :: natm, nmol
    integer :: i, imol, ibeg

    rgsq = 0.0_rp; reedsq = 0.0_rp; span = 0.0_rp

    do imol = 1, num_molecules
        natm = molecules(2,imol)
        if (natm > 1) then
            molbuf = 0.0_rp
            ibeg = molecules(3,imol)
            do i = 1, natm
                molbuf(:,i) = coordinates(:,ibeg+i-1)
            end do
            do i = 1, natm
                molbuf(:,i) = molbuf(:,i) - molbuf(:,1)
                call simbox%get_image(molbuf(:,i))
            end do
            span = span + maxval(molbuf(:,1:natm),2) - minval(molbuf(:,1:natm),2)
            reev = molbuf(:,natm)
            reedsq = reedsq + sum(reev**2)
            rgsq_im = 0.0_rp
            com = sum(molbuf(:,1:natm), 2)/natm
            do i = 1, natm
                rgsq_im = rgsq_im + sum((molbuf(:,i) - com)**2)
            end do
            rgsq = rgsq + rgsq_im/natm
        end if
    end do
    nmol = count(molecules(2,:) > 1)
    span = span/nmol; reedsq = reedsq/nmol; rgsq = rgsq/nmol

    end subroutine

!******************************************************************************

subroutine stats_compute_ic0()
    !! Computes statistics for a single chain in unbounded domain.

    real(rp), dimension(3) :: gt_ev, gt_ev_sc, gt_ev_isc
    real(rp), dimension(3,3) :: S
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: com, com_bbone, com_sc
    real(rp), dimension(3) :: r_tp
    real(rp) :: rgsq_isc, asph_isc, prol_isc
    integer :: n_sc
    integer :: na_bbone 
    integer :: na_sc
    integer :: ibr, ia_beg, ia_end
    integer :: i

    !Number of backbone atoms
    if (num_branches == 0) then
        na_bbone = num_atoms
        n_sc = 0 !Number of side chains
    else
        na_bbone = branches(2,1)
        n_sc = num_branches - 1  !Number of side chains
    end if

    !Zero out
    span = 0.0_rp; rgsq = 0.0_rp; reedsq = 0.0_rp
    asph = 0.0_rp; prol = 0.0_rp
    rgsq_bbone = 0.0_rp 
    rgsq_sc = 0.0_rp; reedsq_sc = 0.0_rp; asph_sc = 0.0_rp; prol_sc = 0.0_rp
    gt_ev = 0.0_rp;
    gt_ev_sc = 0.0_rp;

    S = 0.0_rp; reev = 0.0_rp

    span = maxval(coordinates,2) - minval(coordinates,2)
    reev = coordinates(:,na_bbone) - coordinates(:,1)
    reedsq = reev(1)*reev(1) + reev(2)*reev(2) + reev(3)*reev(3)

    com = sum(coordinates,2)/num_atoms
    do i = 1, num_atoms
        ri = coordinates(:,i) - com
        S(1,1) = S(1,1) + ri(1)*ri(1)
        S(1,2) = S(1,2) + ri(1)*ri(2)
        S(1,3) = S(1,3) + ri(1)*ri(3)
        S(2,2) = S(2,2) + ri(2)*ri(2)
        S(2,3) = S(2,3) + ri(2)*ri(3)
        S(3,3) = S(3,3) + ri(3)*ri(3)
    end do
    S = S/num_atoms
    rgsq = S(1,1) + S(2,2) + S(3,3)

    call dsyevc3(S, gt_ev)
    call calc_shape(gt_ev(1), gt_ev(3), gt_ev(2), asph, prol) 

    com_bbone = sum(coordinates(:,1:na_bbone),2)/na_bbone
    do i = 1, na_bbone
        rgsq_bbone = rgsq_bbone + sum((coordinates(:,i)-com_bbone)**2)
    end do
    rgsq_bbone = rgsq_bbone/na_bbone

    do ibr = 2, num_branches
        ia_beg = branches(3,ibr)
        ia_end = ia_beg + branches(2,ibr) - 1
        r_tp = coordinates(:,branches(1,ibr))

        !Add end-to-end distance squared of this side chain
        reedsq_sc = reedsq_sc + sum((coordinates(:,ia_end) - r_tp)**2)

        !Number of side chain atoms. Adding 1 for the tether point
        na_sc = branches(2,ibr) + 1

        !Get side chain c.o.m.
        com_sc = r_tp + sum(coordinates(:,ia_beg:ia_end),2)
        com_sc = com_sc/na_sc

        !Get radius of gyration squared of this side chain
        S = 0.0_rp
        ri = r_tp - com_sc
        S(1,1) = S(1,1) + ri(1)*ri(1)
        S(1,2) = S(1,2) + ri(1)*ri(2)
        S(1,3) = S(1,3) + ri(1)*ri(3)
        S(2,2) = S(2,2) + ri(2)*ri(2)
        S(2,3) = S(2,3) + ri(2)*ri(3)
        S(3,3) = S(3,3) + ri(3)*ri(3)

        do i = ia_beg, ia_end
            ri = coordinates(:,i) - com_sc
            S(1,1) = S(1,1) + ri(1)*ri(1)
            S(1,2) = S(1,2) + ri(1)*ri(2)
            S(1,3) = S(1,3) + ri(1)*ri(3)
            S(2,2) = S(2,2) + ri(2)*ri(2)
            S(2,3) = S(2,3) + ri(2)*ri(3)
            S(3,3) = S(3,3) + ri(3)*ri(3)
        end do
        S = S/na_sc
        rgsq_isc = S(1,1) + S(2,2) + S(3,3)

        call dsyevc3(S, gt_ev_isc)
        !Note: gt_ev(2) and gt_ev(3) are flipped so that the three eigen values are
        !in descending order
        call calc_shape(gt_ev_isc(1), gt_ev_isc(3), gt_ev_isc(2), &
            asph_isc, prol_isc) 
        !Add radius of gyration, etc. of this side chain
        rgsq_sc = rgsq_sc + rgsq_isc
        gt_ev_sc = gt_ev_sc + gt_ev_isc
        asph_sc = asph_sc + asph_isc
        prol_sc = prol_sc + prol_isc
    end do

    if (n_sc > 0) then
        rgsq_sc = rgsq_sc/n_sc
        reedsq_sc = reedsq_sc/n_sc
        gt_ev_sc = gt_ev_sc/n_sc
        asph_sc = asph_sc/n_sc
        prol_sc = prol_sc/n_sc
    end if

    end subroutine

!******************************************************************************

subroutine calc_shape(ev1, ev2, ev3, asph, prol)
    !! Given three eigen values of the gyration tensor, calculates asphericity
    !! and prolateness. Note that ev1 >= ev2 >= ev3.

    real(rp), intent(in)  :: ev1
    real(rp), intent(in)  :: ev2
    real(rp), intent(in)  :: ev3
    real(rp), intent(out) :: asph
    real(rp), intent(out) :: prol
    real(rp) :: rgsq
    real(rp) :: ev_av, ev1mav, ev2mav, ev3mav, evmav_sumsq

    rgsq = ev1 + ev2 + ev3
    ev_av = rgsq*math_third
    ev1mav = ev1 - ev_av
    ev2mav = ev2 - ev_av
    ev3mav = ev3 - ev_av
    evmav_sumsq = ev1mav**2 + ev2mav**2 + ev3mav**2

    asph = 1.5_rp*evmav_sumsq/(rgsq*rgsq)
    prol = 4*(ev1mav*ev2mav*ev3mav)/((2*math_third)*evmav_sumsq)**1.5

    end subroutine

!******************************************************************************

subroutine dsyevc3(a, w)

    !! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
    !! analytical algorithm.
    !! Only the diagonal and upper triangular parts of A are accessed. The access
    !! is read-only.
    !!
    !! Copyright (C) 2006  Joachim Kopp
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    ! .. Arguments ..
      REAL(RP), DIMENSION(3,3), INTENT(IN) :: A
      REAL(RP), DIMENSION(3), INTENT(OUT) ::  W(3)

     !.. Local Variables ..
      REAL(RP) ::  M, C1, C0
      REAL(RP) ::  DE, DD, EE, FF
      REAL(RP) ::  P, SQRTP, Q, C, S, PHI
  
     !Determine coefficients of characteristic poynomial. We write
     !      | A   D   F  |
     ! A =  | D*  B   E  |
     !      | F*  E*  C  |

      DE    = A(1,2) * A(2,3)
      DD    = A(1,2)**2
      EE    = A(2,3)**2
      FF    = A(1,3)**2
      M     = A(1,1) + A(2,2) + A(3,3)
      C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) &
               - (DD + EE + FF)
      C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) &
               - 2.0_RP * A(1,3)*DE

      P     = M**2 - 3.0_RP * C1
      Q     = M*(P - (3.0_RP/2.0_RP)*C1) - (27.0_RP/2.0_RP)*C0
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0_RP * ( 0.25_RP * C1**2 * (P - C1) &
                + C0 * (Q + (27.0_RP/4.0_RP)*C0) )
      PHI   = (1.0_RP/3.0_RP) * ATAN2(SQRT(ABS(PHI)), Q)

      C     = SQRTP * COS(PHI)
      S     = (1.0_RP/MATH_SQRT3) * SQRTP * SIN(PHI)

      W(2) = (1.0_RP/3.0_RP) * (M - C)
      W(3) = W(2) + S
      W(1) = W(2) + C
      W(2) = W(2) - S

      END SUBROUTINE

!******************************************************************************

end module m_stats_io
