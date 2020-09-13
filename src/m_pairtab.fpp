#: include 'asserts.fypp'

module m_pairtab

use m_precision
use m_vector
use m_table
use m_aabbtree
use m_globals
use m_cell_list
use m_connectivity

implicit none

private
public :: pair_tab
public :: pt_init, pt_finish, pt_build

type(itable_t)   :: pair_tab
type(aabbtree_t) :: tree
real(rp) :: rskin_sq = 0.0_rp
real(rp) :: tskin_sq = 0.0_rp
real(rp), dimension(:,:), allocatable :: coordinates_save
    !!  (3, *num_atoms*) array
real(rp), dimension(:,:), allocatable :: coordinates_dr
    !!  (3, *num_atoms*) array

contains

!*******************************************************************************
 
subroutine pt_init() 
    !! Performs initial setup for building a pair list.

    integer :: i

    !Initialize and build atom -> bond table, stored in the module variable
    !atbo_tab in module m_connectivity. The atom -> bond table is required for
    !building the table for excluded atoms.
    call atbo_build()

    !Build excluded atoms table, stored in the module variable exat_tab in
    !module m_connectivity.
    call exat_build()

    !Delete atom -> bond table, it is no longer needed.
    call atbo_tab%delete()

    !Initialize the pair table
    call itbl_init(pair_tab, num_atoms)
    !The last row of pair_tab will remain empty due to Newton's third law;
    !this will be used in an assert statement.

    if (use_verlet_tab) then
        !Setup the verlet scheme
        rskin_sq = (rcutoff+tskin)**2; tskin_sq = tskin**2
        allocate( coordinates_save(3,num_atoms) )
        allocate( coordinates_dr(3,num_atoms)   )
        !Initializing to zero
        coordinates_save = 0.0_rp
        coordinates_dr = 0.0_rp

    else if (use_aabbtree) then
        !Setup an AABB tree and insert the atoms
        call tree%init(num_atoms, tskin/rcutoff)
        !Insert atoms
        do i = 1, num_atoms
            call tree%insert(i, coordinates(:,i), rcutoff)
        end do
        !Build an optimized tree
        call tree%rebuild()

    else if (use_cell_list) then
        !Setup a cell list
        if ( (sim_style == 0) .or. (sim_style == 1) ) then
            ! For relaxation/BD simulation the cell list is used for short range
            ! force calculation.
            call cl_init(num_atoms, rcutoff)
            call cl_set_cell_size(rcutoff)
            call cl_build_cell_nbrs()
        else if (sim_style == 2) then
            ! For MPCD simulation the cell list is always used for sorting. If
            ! use_cell_list == .true., the cell list will also be used for force
            ! calculation.
            call cl_init(num_atoms_tot, 1.0_rp) !Collision cell size is 1.0.
            call cl_set_cell_size(rcutoff) !Set cell size to global cutoff 
            call cl_build_cell_nbrs()
        end if
    else
        !Perform a direct N^2 calculation
        rskin_sq = rcutoff**2
    end if

    end subroutine
 
!*******************************************************************************
 
subroutine pt_finish() 
    !! Cleanup for pair list calculation.

    if (use_verlet_tab) then
        if (allocated(coordinates_save)) deallocate(coordinates_save)
        if (allocated(coordinates_dr)  ) deallocate(coordinates_dr)
        rskin_sq = 0.0_rp; tskin_sq = 0.0_rp
    else if (use_aabbtree ) then
        call tree%delete()
    else if (use_cell_list) then
        call cl_delete()
    end if

    call pair_tab%delete()
    call exat_tab%delete()

    end subroutine
 
!*******************************************************************************
 
subroutine pt_build() 
    !! Builds a pair table. The resulting table is stored in the module variable
    !! `pair_tab`.

    if (use_verlet_tab) then
        call build_pt_verlet()
    else if (use_aabbtree) then
        call build_pt_aabbtree()
    else if (use_cell_list) then
        call build_pt_cell_list()
    else
        call build_pt_n2()
    end if

    end subroutine
 
!******************************************************************************

subroutine build_pt_n2()
    !! Builds a pair table using direct N^2 looping over all pairs.

    integer, dimension(:), pointer :: pnbrs => null()
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rj
    real(rp), dimension(3) :: rij
    real(rp) :: rij_sq
    integer :: i, j

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
            if (imcon /= 0) call simbox%get_image(rij)
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

subroutine build_pt_verlet()
    !! Builds a pair table using the Verlet scheme.

    integer, dimension(:), pointer :: pnbrs => null()
    real(rp), dimension(3) :: ri
    real(rp), dimension(3) :: rj
    real(rp), dimension(3) :: rij
    real(rp) :: dr_sq_max
    real(rp) :: rij_sq
    logical, save :: first_call = .true.
    integer :: i, j

    !On first call no check for rebuilding
    if (.not. first_call) then
        !Check whether rebuilding the list is necessary
        coordinates_dr = coordinates(:,1:num_atoms) - coordinates_save
        if (imcon /= 0) then
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
            if (imcon /= 0) call simbox%get_image(rij)
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

subroutine build_pt_aabbtree()

    type(ivector_t) :: nbrs
    integer, dimension(:), pointer :: pnbrs => null()
    integer :: ia, inbr, ia_nbr
    logical :: lstat

    !Clear out the pair table
    call pair_tab%clear()

    !Update the AABB tree with current atom positions
    !do ia = 1, num_atoms
    !    call tree%update(ia, coordinates(:,ia), rcutoff, lstat)
    !end do
    !Insert atoms
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

subroutine build_pt_cell_list()
    !! Builds a pair table using cell list.

    type(ivector_t), dimension(:), allocatable :: ptaov
        !Pair table as an array of vectors.
    integer, dimension(:), pointer :: aic => null()
    integer, dimension(:), pointer :: ainc => null()
    integer, dimension(:), pointer :: nbr_cells => null()
    integer, dimension(:), pointer :: pnbrs => null()
    integer :: num_cells
    integer :: icell, jcell, iatm, jatm
    integer :: i, j, k

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

end module m_pairtab
