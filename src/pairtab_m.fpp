#: include 'asserts.fypp'

module pairtab_m

use constants_m
use vector_m
use table_m
use simbox_m
use aabbtree_m
use cell_list_m
use atmcfg_m
use connectivity_m

implicit none

private
public :: pt_init, pt_delete, pt_build

real(rp) :: rcutoff = 0.0_rp
real(rp) :: tskin = 0.0_rp
real(rp) :: rskin_sq = 0.0_rp
real(rp) :: tskin_sq = 0.0_rp
character(len=:), allocatable :: mth_ptgen
    !! Pair table generation method: {'DIR', 'VER', 'AABBT', 'CL'}
type(itable_t) :: exat_tab
type(aabbtree_t) :: tree
real(rp), dimension(:,:), allocatable :: coordinates_save
    !!  (3, *num_atoms*) array
real(rp), dimension(:,:), allocatable :: coordinates_dr
    !!  (3, *num_atoms*) array

contains

!*******************************************************************************
 
subroutine pt_init(mth, num_atoms, excl_atoms, rcut, tskn, bonds, &
        simbox, pair_tab)
    !! Performs initial setup for building a pair list.

    character(len=*), intent(in) :: mth
    integer, intent(in) :: num_atoms
    integer, intent(in) :: excl_atoms
    real(rp), intent(in) :: rcut
    real(rp), intent(in) :: tskn
    integer, dimension(:,:), intent(in) :: bonds
    type(smbx_t), intent(in) :: simbox
    type(itable_t), intent(in out) :: pair_tab
    type(itable_t) :: atbo_tab
    type(itable_t) :: atat_tab

    mth_ptgen = mth; rcutoff = rcut; tskin = tskn

    !Build atom -> bond table
    call atbo_build(num_atoms, size(bonds,2), bonds, atbo_tab)

    !Build atm -> bonded atoms table
    call atat_build(num_atoms, bonds, atbo_tab, atat_tab)

    !Build excluded atoms table
    call exat_build(num_atoms, excl_atoms, atat_tab, exat_tab)

    !Delete tables no longer needed.
    call atbo_tab%delete()
    call atat_tab%delete()

    !Initialize the pair table
    call itbl_init(pair_tab, num_atoms)
    !The last row of pair_tab will remain empty due to Newton's third law;
    !this will be used in an assert statement.

    select case(mth_ptgen)
    case('VER')
        !Setup the verlet scheme
        rskin_sq = (rcutoff+tskin)**2; tskin_sq = tskin**2
        allocate( coordinates_save(3,num_atoms) )
        allocate( coordinates_dr(3,num_atoms)   )
        !Initializing to zero
        coordinates_save = 0.0_rp
        coordinates_dr = 0.0_rp

    case('AABBT')
        !Setup an AABB tree and insert the atoms
        call tree%init(num_atoms, tskin/rcutoff)

    case('CL')
        !Setup a cell list
        ! For relaxation/BD simulation the cell list is used for short range
        ! force calculation.
        call cl_init(num_atoms, rcutoff, simbox)
        call cl_set_cell_size(rcutoff, simbox)
        call cl_build_cell_nbrs()

    case default
        !Perform a direct N^2 calculation
        rskin_sq = rcutoff**2
    end select

    end subroutine
 
!*******************************************************************************
 
subroutine pt_delete(pair_tab) 
    !! Cleanup for pair list calculation.

    type(itable_t), intent(in out) :: pair_tab

    select case(mth_ptgen)
    case('VER')
        if (allocated(coordinates_save)) deallocate(coordinates_save)
        if (allocated(coordinates_dr)  ) deallocate(coordinates_dr)
        rskin_sq = 0.0_rp; tskin_sq = 0.0_rp
    case('AABBT')
        call tree%delete()
    case('CL')
        call cl_delete()
    case default
        rskin_sq = 0.0_rp
    end select

    rcutoff = 0.0_rp; tskin = 0.0_rp

    call pair_tab%delete()
    call exat_tab%delete()

    end subroutine
 
!*******************************************************************************
 
subroutine pt_build(simbox, coordinates, pair_tab)
    !! Builds a pair table. The resulting table is stored in the module variable
    !! `pair_tab`.

    type(smbx_t), intent(in) :: simbox
    real(rp), dimension(:,:), intent(in) :: coordinates
    type(itable_t), intent(in out) :: pair_tab

    select case(mth_ptgen)
    case('VER')
        call build_pt_verlet(simbox, coordinates, pair_tab)
    case('AABBT')
        call build_pt_aabbtree(simbox, coordinates, pair_tab)
    case('CL')
        call build_pt_cell_list(simbox, coordinates, pair_tab)
    case default
        call build_pt_n2(simbox, coordinates, pair_tab)
    end select

    end subroutine
 
!******************************************************************************

subroutine build_pt_n2 (simbox, coordinates, pair_tab)
    !! Builds a pair table using direct N^2 looping over all pairs.

    type(smbx_t), intent(in) :: simbox
    real(rp), dimension(:,:), intent(in) :: coordinates
    type(itable_t), intent(in out) :: pair_tab
    integer, dimension(:), pointer :: pnbrs => null()
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rj
    real(rp), dimension(3) :: rij
    real(rp) :: rij_sq
    integer :: i, j, num_atoms

    !Clear out the pair table
    call pair_tab%clear()

    num_atoms = size(coordinates,2)

    !Loop over all pairs to build list
    do i = 1, num_atoms-1
        ri = coordinates(:,i)
        do j = i+1, num_atoms
            !Check if atom j is an excluded atom for atom i
            if ( exat_tab%is_in(i,j) ) cycle
            rj = coordinates(:,j)
            rij = rj - ri
            if (simbox%imcon /= 0) call simbox%get_image(rij)
            rij_sq = sum(rij**2)
            if (rij_sq < rskin_sq) call pair_tab%append(i,j)
        end do
    end do

    call pair_tab%get_row(num_atoms,pnbrs)
    @:ASSERT(size(pnbrs) == 0)

    !Release additional memory
    call pair_tab%shrink_to_fit()

    end subroutine

!******************************************************************************

subroutine build_pt_verlet(simbox, coordinates, pair_tab)
    !! Builds a pair table using the Verlet scheme.

    type(smbx_t), intent(in) :: simbox
    real(rp), dimension(:,:), intent(in) :: coordinates
    type(itable_t), intent(in out) :: pair_tab
    integer, dimension(:), pointer :: pnbrs => null()
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rj
    real(rp), dimension(3) :: rij
    real(rp) :: dr_sq_max
    real(rp) :: rij_sq
    logical, save :: first_call = .true.
    integer :: i, j, num_atoms

    num_atoms = size(coordinates,2)

    !On first call no check for rebuilding
    if (.not. first_call) then
        !Check whether rebuilding the list is necessary
        coordinates_dr = coordinates(:,1:num_atoms) - coordinates_save
        if (simbox%imcon /= 0) then
            do i = 1, num_atoms
                call simbox%get_image(coordinates_dr(:,i))
            end do
        end if
        dr_sq_max = maxval(sum(coordinates_dr**2, dim=1))
        if ( 4*dr_sq_max < tskin_sq ) return
    end if

    first_call = .false.

    !Clear out the pair table
    call pair_tab%clear()

    !Loop over all pairs to build list
    do i = 1, num_atoms-1
        ri = coordinates(:,i)
        do j = i+1, num_atoms
            !Check if atom j is an excluded atom for atom i
            if ( exat_tab%is_in(i,j) ) cycle
            rj = coordinates(:,j)
            rij = rj - ri
            if (simbox%imcon /= 0) call simbox%get_image(rij)
            rij_sq = sum(rij**2)
            if (rij_sq < rskin_sq) call pair_tab%append(i,j)
        end do
    end do

    call pair_tab%get_row(num_atoms,pnbrs)
    @:ASSERT(size(pnbrs) == 0)

    !Release additional memory
    call pair_tab%shrink_to_fit()

    !Back up positions
    coordinates_save = coordinates

    end subroutine

!*******************************************************************************

subroutine build_pt_aabbtree(simbox, coordinates, pair_tab)

    type(smbx_t), intent(in) :: simbox
    real(rp), dimension(:,:), intent(in) :: coordinates
    type(itable_t), intent(in out) :: pair_tab
    type(ivector_t) :: nbrs
    integer, dimension(:), pointer :: pnbrs => null()
    integer :: ia, inbr, ia_nbr, num_atoms

    num_atoms = size(coordinates,2)

    !Clear out the pair table
    call pair_tab%clear()

    !Clear out the AABB tree
    call tree%clear()
    !Insert atoms into the tree
    do ia = 1, num_atoms
        call tree%insert(ia, coordinates(:,ia), rcutoff)
    end do

    !Loop over all atoms to find neighbors and insert into pair_tab.
    call ivector_init(nbrs)
    do ia = 1, num_atoms
        call tree%query(ia, nbrs)
        call nbrs%sort()
        call nbrs%get_data(pnbrs)
        do inbr = 1, size(pnbrs)
            ia_nbr = pnbrs(inbr)
            !Check if atom ia_nbr is an excluded atom for atom ia
            if ( exat_tab%is_in(ia,ia_nbr) ) cycle
            if (.not. pair_tab%is_in(ia_nbr,ia)) then
                call pair_tab%append(ia, ia_nbr)
            end if
        end do
    end do

    call pair_tab%get_row(num_atoms,pnbrs)
    @:ASSERT(size(pnbrs) == 0)

    call nbrs%delete()

    !Release additional memory
    call pair_tab%shrink_to_fit()

    end subroutine

!******************************************************************************

subroutine build_pt_cell_list (simbox, coordinates, pair_tab)
    !! Builds a pair table using cell list.

    type(smbx_t), intent(in) :: simbox
    real(rp), dimension(:,:), intent(in) :: coordinates
    type(itable_t), intent(in out) :: pair_tab
    type(ivector_t), dimension(:), allocatable :: ptaov
        !Pair table as an array of vectors.
    integer, dimension(:), pointer :: aic => null()
    integer, dimension(:), pointer :: ainc => null()
    integer, dimension(:), pointer :: nbr_cells => null()
    integer, dimension(:), pointer :: pnbrs => null()
    integer :: num_cells, num_atoms
    integer :: icell, jcell, iatm, jatm
    integer :: i, j, k

    num_atoms = size(coordinates, 2)

    !Clear out the pair table
    call pair_tab%clear()

    !Initialize ptaov
    allocate(ptaov(num_atoms))
    do i = 1, size(ptaov)
        call ivector_init(ptaov(i))
    end do

    call cl_build(coordinates)
    num_cells = cl_get_num_cells()

    do icell = 0, (num_cells-1)
        call cl_get_contents(icell, aic)
        call cl_get_nbr_cells(icell, nbr_cells)

        !Interaction with particles within the cell
        do i = 1, size(aic) - 1
            iatm = aic(i)
            do j = i+1, size(aic)
                jatm = aic(j)
                !Check if atom jatm is an excluded atom for atom iatm
                if ( exat_tab%is_in(iatm,jatm) ) cycle
                call ptaov(iatm)%append(jatm)
            end do
        end do

        !Interaction with particles belonging to neighboring cells
        do i = 1, size(aic)
            iatm = aic(i)
            do j = 1, size(nbr_cells)
                jcell = nbr_cells(j)
                call cl_get_contents(jcell, ainc)
                do k = 1, size(ainc)
                    jatm = ainc(k)
                    !Check if atom jatm is an excluded atom for atom iatm
                    if ( exat_tab%is_in(iatm,jatm) ) cycle
                    call ptaov(iatm)%append(jatm)
                end do
            end do
        end do
    end do

    !Copy the entries in ptaov to pair_tab
    do iatm = 1, num_atoms
        call ptaov(iatm)%get_data(pnbrs)
        do j = 1, size(pnbrs)
            call pair_tab%append(iatm,pnbrs(j))
        end do
    end do

    call pair_tab%get_row(num_atoms,pnbrs)
    @:ASSERT(size(pnbrs) == 0)

    !Release additional memory
    call pair_tab%shrink_to_fit()

    do i = 1, size(ptaov)
        call ptaov(i)%delete()
    end do
    deallocate(ptaov)

    end subroutine

!*******************************************************************************

end module pairtab_m
