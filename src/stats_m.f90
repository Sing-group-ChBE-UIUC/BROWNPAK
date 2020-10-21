module stats_m
    !! Computes and writes properties calculated during simulation.

use constants_m
use strings_m
use control_m
use simbox_m
use atmcfg_m

implicit none

private

public :: stats_init, stats_finish, stats_write, stats_accumulate, &
            bndlen, bndlen_min, bndlen_max, energy_bond, energy_angle, &
            energy_dihedral, energy_vdw, energy_tether, energy_external, &
            energy_tot, stress

character(len=:), allocatable :: fn_stats
integer  :: fu_stats
integer  :: bnsiz, bnsiz_stress

integer  :: bnpop, bnpop_stress

real(rp) :: bndlen
real(rp) :: bndlen_min
real(rp) :: bndlen_max

real(rp) :: energy_bond,     energy_bond_accu
real(rp) :: energy_angle,    energy_angle_accu
real(rp) :: energy_dihedral, energy_dihedral_accu
real(rp) :: energy_vdw,      energy_vdw_accu
real(rp) :: energy_tether,   energy_tether_accu
real(rp) :: energy_external, energy_external_accu
real(rp) :: energy_tot,      energy_tot_accu

real(rp), dimension(3,3) :: stress, stress_accu

real(rp) :: rgsq,   rgsq_accu
real(rp) :: reedsq, reedsq_accu
real(rp), dimension(3) :: span, span_accu

real(rp) :: rgsq_bbone, rgsq_bbone_accu
real(rp) :: rgsq_sc,    rgsq_sc_accu
real(rp) :: reedsq_sc,  reedsq_sc_accu
real(rp), dimension(3) :: span_bbone, span_bbone_accu

real(rp), dimension(3) :: reev

real(rp), dimension(:,:), allocatable :: molbuf
    ! Buffer for the largest molecule in the system

logical :: lvdw

contains

!******************************************************************************

subroutine stats_init(cpar, simbox, atc, job_tag)
    !! Set up for stats collection

    type(ctrlpar_t), intent(in) :: cpar
    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in) :: atc
    character(len=*), intent(in) :: job_tag
    logical :: lexists
    integer :: n

    bnsiz = cpar%stats_binsize(1); bnsiz_stress = cpar%stats_binsize(2)
    lvdw = .true.
    if ( (.not. cpar%lvdw) .or. (atc%num_vdw_types==0) ) lvdw = .false.

    !Open/create file for writing statistics
    fn_stats = cpar%fn_stats//job_tag
    if (cpar%lrevive) then
        !Restarting simulation
        !Check if file exists
        inquire(file=fn_stats, exist=lexists)
        if (lexists) then
            !Open existing file for appending statistics
            open(newunit=fu_stats, file=fn_stats, action='write', &
                position='append', status='old')
        else
            !Open new file for writing statistics
            open(newunit=fu_stats, file=fn_stats, action='write', status='new')
            !Write header
            call write_hdr(simbox, atc)
        end if
    else
        !New simulation
        !Open file for writing production (or relaxation) statistics
        open(newunit=fu_stats, file=fn_stats, action='write', status='replace')
        !Write header
        call write_hdr(simbox, atc)
    end if

    !Allocate space for unwrapping molecule under PBC
    if ( (simbox%imcon /= 0) .and. (atc%num_molecules > 0) ) then
        n = maxval(atc%molecules(2,:))
        if (n > 1) allocate( molbuf(3,n) )
    end if

    call zero_out()

    end subroutine

!******************************************************************************

subroutine stats_finish()
    !! Clean up for stats collection.

    logical :: lopened

    if (allocated(fn_stats)) then
        inquire(file=fn_stats, number=fu_stats, opened=lopened)
        if (lopened) close(fu_stats)
        deallocate(fn_stats)
    end if

    if (allocated(molbuf)) deallocate(molbuf)

    call zero_out()

    end subroutine

!******************************************************************************

subroutine zero_out()

    bnpop = 0; bnpop_stress = 0
    bndlen = 0.0_rp
    bndlen_min = 0.0_rp
    bndlen_max = 0.0_rp
    
    energy_bond = 0.0_rp;     energy_bond_accu = 0.0_rp
    energy_angle = 0.0_rp;    energy_angle_accu = 0.0_rp
    energy_dihedral = 0.0_rp; energy_dihedral_accu = 0.0_rp
    energy_vdw = 0.0_rp;      energy_vdw_accu = 0.0_rp
    energy_tether = 0.0_rp;   energy_tether_accu = 0.0_rp
    energy_external = 0.0_rp; energy_external_accu = 0.0_rp
    energy_tot = 0.0_rp;      energy_tot_accu = 0.0_rp
    
    stress = 0.0_rp; stress_accu = 0.0_rp
    
    rgsq = 0.0_rp;   rgsq_accu = 0.0_rp
    reedsq = 0.0_rp; reedsq_accu = 0.0_rp
    span = 0.0_rp; span_accu = 0.0_rp
    
    rgsq_bbone = 0.0_rp; rgsq_bbone_accu = 0.0_rp
    rgsq_sc = 0.0_rp;    rgsq_sc_accu = 0.0_rp
    reedsq_sc = 0.0_rp;  reedsq_sc_accu = 0.0_rp
    span_bbone = 0.0_rp; span_bbone_accu = 0.0_rp
    
    reev = 0.0_rp

    end subroutine

!******************************************************************************

subroutine write_hdr(simbox, atc)
    !! Writes header of the file `fn_stats`.

    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in) :: atc

    write(fu_stats,'(a12)', advance='no') 'nts'

    if (atc%num_bonds > 0) then
        write(fu_stats, 200, advance='no') &
            'bndlen_avg', 'bndlen_min', 'bndlen_max', 'enrg_bnd'
    end if

    if (atc%num_angles > 0) then
        write(fu_stats, 100, advance='no') 'enrg_ang'
    end if

    if (atc%num_dihedrals > 0) then
        write(fu_stats, 100, advance='no') 'enrg_dhd'
    end if

    if (lvdw) then
        write(fu_stats, 100, advance='no') 'enrg_vdw'
    end if

    if (atc%num_tethers > 0) then
        write(fu_stats, 100, advance='no') 'enrg_teth'
    end if

    if (atc%num_externals > 0) then
        write(fu_stats, 100, advance='no') 'enrg_ext'
    end if

    write(fu_stats, 100, advance='no') 'enrg_tot'

    !Stress
    if (atc%num_atoms > 0) then
        write(fu_stats, 200, advance='no') 'sxx', 'syx', 'szx', &
            'syy', 'szy', 'szz'
    end if

    !Structural properties: Size 
    if ( atc%num_bonds > 0 ) then
        write(fu_stats, 200, advance='no') 'rgsq', 'reedsq',  &
            'spanx', 'spany', 'spanz'
    end if

    !Further structural properties (for single molecule in unbounded domain)
    if ( (simbox%imcon == 0) .and. (atc%num_bonds > 0) ) then
        write(fu_stats, 200, advance='no') 'reevx', 'reevy', 'reevz'
        if (atc%num_branches > 0) then
            write(fu_stats, 200, advance='no') 'rgsq_bbone',  &
                'rgsq_sc', 'reedsq_sc', 'span_bbone_x', 'span_bbone_y',&
                'span_bbone_z'
        end if
        if (atc%flow_style == 0) then
            write(fu_stats, 200, advance='no') 'comx', 'comy', 'comz'
        end if
    end if

    write(fu_stats,*)

100 format(1x,a14)
200 format(*(1x,a14))

    end subroutine

!******************************************************************************

subroutine stats_write(nts, simbox, atc)
    !! Writes statistics to `fn_stats`.

    integer(ip_long), intent(in) :: nts
    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in) :: atc

    write(fu_stats,'(i12)', advance='no') nts

    if (atc%num_bonds > 0) then
        write(fu_stats, 200, advance='no') &
            bndlen, bndlen_min, bndlen_max, energy_bond_accu
    end if

    if (atc%num_angles > 0) then
        write(fu_stats, 100, advance='no') energy_angle_accu
    end if

    if (atc%num_dihedrals > 0) then
        write(fu_stats, 100, advance='no') energy_dihedral_accu
    end if

    if (lvdw) then
        write(fu_stats, 100, advance='no') energy_vdw_accu
    end if

    if (atc%num_tethers > 0) then
        write(fu_stats, 100, advance='no') energy_tether_accu
    end if

    if (atc%num_externals > 0) then
        write(fu_stats, 100, advance='no') energy_external_accu
    end if

    write(fu_stats, 100, advance='no') energy_tot_accu

    !Stress (calculated only during production run)
    if (atc%num_atoms > 0) then
        write(fu_stats, 200, advance='no') stress_accu(1,1), &
            stress_accu(2,1), stress_accu(3,1), stress_accu(2,2), &
            stress_accu(3,2), stress_accu(3,3)
    end if

    !Structural properties: Size
    if ( atc%num_bonds > 0 ) then
        write(fu_stats, 200, advance='no') rgsq_accu, &
            reedsq_accu, span_accu
    end if

    !Further structural properties (for single molecule in unbounded domain)
    if ( (simbox%imcon == 0) .and. (atc%num_bonds > 0) ) then
        write(fu_stats, 200, advance='no') reev
        if (atc%num_branches > 0) then
            write(fu_stats, 200, advance='no') rgsq_bbone_accu, &
                rgsq_sc_accu, reedsq_sc_accu, span_bbone_accu
        end if
        if (atc%flow_style == 0) then
            write(fu_stats, 200, advance='no') atc%molecule_com
        end if
    end if

    write(fu_stats,*)

100 format(1x,es14.6)
200 format(*(1x,es14.6))

    end subroutine

!******************************************************************************

subroutine stats_accumulate(simbox, atc)
    !! Accumulates statistics.

    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in) :: atc

    if (simbox%imcon == 0) then
        call compute_ic0(atc)
    else
        call compute_ic1(simbox, atc)
    end if

    if (bnpop == bnsiz) then
        bnpop = 0
        energy_bond_accu = 0.0_rp
        energy_angle_accu = 0.0_rp
        energy_dihedral_accu = 0.0_rp
        energy_vdw_accu = 0.0_rp
        energy_tether_accu = 0.0_rp
        energy_external_accu = 0.0_rp
        energy_tot_accu = 0.0_rp

        rgsq_accu = 0.0_rp
        reedsq_accu = 0.0_rp
        span_accu = 0.0_rp

        rgsq_bbone_accu = 0.0_rp
        rgsq_sc_accu = 0.0_rp
        reedsq_sc_accu = 0.0_rp
        span_bbone_accu = 0.0_rp
    end if

    bnpop = bnpop + 1
    energy_bond_accu = energy_bond_accu &
        + ((energy_bond-energy_bond_accu)/bnpop)
    energy_angle_accu = energy_angle_accu &
        + ((energy_angle-energy_angle_accu)/bnpop)
    energy_dihedral_accu = energy_dihedral_accu &
        + ((energy_dihedral-energy_dihedral_accu)/bnpop)
    energy_vdw_accu = energy_vdw_accu &
        + ((energy_vdw-energy_vdw_accu)/bnpop)
    energy_tether_accu = energy_tether_accu &
        + ((energy_tether-energy_tether_accu)/bnpop)
    energy_external_accu = energy_external_accu &
        + ((energy_external-energy_external_accu)/bnpop)
    energy_tot_accu = energy_tot_accu &
        + ((energy_tot-energy_tot_accu)/bnpop)

    rgsq_accu = rgsq_accu + ((rgsq-rgsq_accu)/bnpop)
    reedsq_accu = reedsq_accu + ((reedsq-reedsq_accu)/bnpop)
    span_accu = span_accu + ((span-span_accu)/bnpop)

    rgsq_bbone_accu = rgsq_bbone_accu + ((rgsq_bbone-rgsq_bbone_accu)/bnpop)
    rgsq_sc_accu = rgsq_sc_accu + ((rgsq_sc-rgsq_sc_accu)/bnpop)
    reedsq_sc_accu = reedsq_sc_accu + ((reedsq_sc-reedsq_sc_accu)/bnpop)
    span_bbone_accu = span_bbone_accu + ((span_bbone-span_bbone_accu)/bnpop)

    if (bnpop_stress == bnsiz_stress) then
        bnpop_stress = 0
        stress_accu = 0.0_rp
    end if
    bnpop_stress = bnpop_stress + 1
    stress_accu = stress_accu + ( (stress-stress_accu)/bnpop_stress )

    end subroutine

!******************************************************************************

subroutine compute_ic1(simbox, atc)
    !! Computes statistics for possibly multiple chains in a periodic domain.

    type(smbx_t), intent(in) :: simbox
    type(atmcfg_t), intent(in) :: atc
    real(rp), dimension(3) :: com
    real(rp) :: rgsq_im
    integer :: natm, nmol
    integer :: i, imol, ibeg

    rgsq = 0.0_rp; reedsq = 0.0_rp; span = 0.0_rp

    do imol = 1, atc%num_molecules
        natm = atc%molecules(2,imol)
        if (natm > 1) then
            molbuf = 0.0_rp
            ibeg = atc%molecules(3,imol)
            do i = 1, natm
                molbuf(:,i) = atc%coordinates(:,ibeg+i-1)
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
    nmol = count(atc%molecules(2,:) > 1)
    span = span/nmol; reedsq = reedsq/nmol; rgsq = rgsq/nmol

    end subroutine

!******************************************************************************

subroutine compute_ic0(atc)
    !! Computes statistics for a single chain in unbounded domain.

    type(atmcfg_t), intent(in) :: atc
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: com, com_bbone, com_sc
    real(rp), dimension(3) :: r_tp
    real(rp) :: rgsq_isc
    integer :: n_sc, na_bbone, na_sc
    integer :: i, ibr, ia_beg, ia_end

    !Number of backbone atoms
    if (atc%num_branches == 0) then
        na_bbone = atc%num_atoms
        n_sc = 0 !Number of side chains
    else
        na_bbone = atc%branches(2,1)
        n_sc = atc%num_branches - 1  !Number of side chains
    end if

    !Zero out
    rgsq = 0.0_rp; reedsq = 0.0_rp; span = 0.0_rp 
    rgsq_bbone = 0.0_rp; rgsq_sc = 0.0_rp; reedsq_sc = 0.0_rp
    span_bbone = 0.0_rp
    reev = 0.0_rp

    span = maxval(atc%coordinates,2) - minval(atc%coordinates,2)
    span_bbone = maxval(atc%coordinates(:,1:na_bbone),2) &
        - minval(atc%coordinates(:,1:na_bbone),2)
    reev = atc%coordinates(:,na_bbone) - atc%coordinates(:,1)
    reedsq = sum(reev**2)

    com = sum(atc%coordinates,2)/atc%num_atoms
    do i = 1, atc%num_atoms
        ri = atc%coordinates(:,i) - com
        rgsq = rgsq + sum(ri**2)
    end do
    rgsq = rgsq/atc%num_atoms

    com_bbone = sum(atc%coordinates(:,1:na_bbone),2)/na_bbone
    do i = 1, na_bbone
        rgsq_bbone = rgsq_bbone + sum((atc%coordinates(:,i)-com_bbone)**2)
    end do
    rgsq_bbone = rgsq_bbone/na_bbone

    do ibr = 2, atc%num_branches
        ia_beg = atc%branches(3,ibr)
        ia_end = ia_beg + atc%branches(2,ibr) - 1
        r_tp = atc%coordinates(:,atc%branches(1,ibr))

        !Add end-to-end distance squared of this side chain
        reedsq_sc = reedsq_sc + sum((atc%coordinates(:,ia_end) - r_tp)**2)

        !Number of side chain atoms. Adding 1 for the tether point
        na_sc = atc%branches(2,ibr) + 1

        !Get side chain c.o.m.
        com_sc = r_tp + sum(atc%coordinates(:,ia_beg:ia_end),2)
        com_sc = com_sc/na_sc

        !Get radius of gyration squared of this side chain
        rgsq_isc = 0.0_rp
        ri = r_tp - com_sc
        rgsq_isc = rgsq_isc + sum(ri**2)

        do i = ia_beg, ia_end
            ri = atc%coordinates(:,i) - com_sc
            rgsq_isc = rgsq_isc + sum(ri**2)
        end do
        rgsq_isc = rgsq_isc/na_sc
        !Add radius of gyration of this side chain
        rgsq_sc = rgsq_sc + rgsq_isc
    end do

    if (n_sc > 0) then
        rgsq_sc = rgsq_sc/n_sc
        reedsq_sc = reedsq_sc/n_sc
    end if

    end subroutine

!******************************************************************************

end module stats_m
