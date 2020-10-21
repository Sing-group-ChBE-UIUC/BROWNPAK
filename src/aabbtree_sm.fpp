#: include 'asserts.fypp'

submodule (aabbtree_m) aabbtree_sm

implicit none

contains

!*******************************************************************************

module subroutine node_init(this)

    class(node_t), intent(out) :: this
    real(rp), dimension(3) :: lbnd, ubnd

    this%next = NULL_NODE
    this%parent = NULL_NODE; this%left = NULL_NODE; this%right = NULL_NODE
    this%height = -1; this%atom = 0

    lbnd = 0.0_rp; ubnd = 0.0_rp
    call this%aabb%init(lbnd, ubnd)

    end subroutine

!*******************************************************************************

module function node_isleaf(this) result(res)
    class(node_t), intent(in) :: this
    logical :: res

    res = (this%left == null_node)

    end function

!*******************************************************************************

module function node_asstr(this, frmt) result(buf)
    class(node_t), intent(in) :: this
    character(len=*), intent(in), optional :: frmt
    character(len=:), allocatable :: buf
    character(len=:), allocatable :: aabb_str

    buf = '{'
    buf = buf // 'next: '//str_from_num(this%next)
    buf = buf // ', parent: '//str_from_num(this%parent)
    buf = buf // ', left: '//str_from_num(this%left)
    buf = buf // ', right: '//str_from_num(this%right)
    buf = buf // ', height: '//str_from_num(this%height)
    buf = buf // ', atom: '//str_from_num(this% atom)

    if (present(frmt)) then
        call this%aabb%print(frmt=frmt, str=aabb_str)
    else
        call this%aabb%print(str=aabb_str)
    end if
    buf = buf // ', aabb: {' // aabb_str // '}'
    buf = buf // '}'

    end function

!*******************************************************************************

module subroutine init(this, natoms, tskin)

    class(aabbtree_t), intent(out) :: this
    integer, intent(in) :: natoms
    real(rp), intent(in) :: tskin
    integer :: i
    
    this%capacity = natoms; this%tskin = tskin
    this%size = 0
    allocate(this%nodes(this%capacity))
    !Initially all nodes belong the free store.
    do i = 1, this%capacity
        call this%nodes(i)%init()
        this%nodes(i)%next = i + 1
    end do
    this%nodes(this%capacity)%next = NULL_NODE
    this%freestore = 1 !Head of the free store points to the first node
    this%root = NULL_NODE
    allocate(this%atnd_tab(natoms)); this%atnd_tab = 0

    !Initialize the stack
    call ivector_init(stack)

    end subroutine

!*******************************************************************************

module recursive subroutine print(this, p)

    class(aabbtree_t), intent(in) :: this
    integer, intent(in), optional :: p
    integer :: p_

    p_ = this%root
    if (present(p)) p_ = p

    !Nothing to do for an empty tree
    if (p_ == NULL_NODE) return
    !If the tree is not empty
    call this%print( this%nodes(p_)%left )
    write(*,'(i0,": ",a)') p_, this%nodes(p_)%asstr()
    call this%print( this%nodes(p_)%right )

    end subroutine

!*******************************************************************************

    module subroutine clear(this)

        class(aabbtree_t), intent(in out) :: this
        integer :: i
    
        !Return all nodes in the tree to the free store
        do i = 1, this%capacity
            if (this%nodes(i)%height < 0) cycle
            call this%fs_return(i)
        end do

        this%atnd_tab = 0; this%tskin = 0.0_rp
        this%root = NULL_NODE

        @: ASSERT( this%size == 0 )

        call stack%clear()

        end subroutine

!*******************************************************************************

    module subroutine delete(this)

        class(aabbtree_t), intent(in out) :: this

        if (allocated(this%nodes)) deallocate(this%nodes)
        if (allocated(this%atnd_tab)) deallocate(this%atnd_tab)

        this%capacity = 0; this%size = 0; this%tskin = 0.0_rp
        this%root = NULL_NODE; this%freestore = NULL_NODE
        call stack%delete()

        end subroutine

!*******************************************************************************

module subroutine insert(this, ia, pos, radius)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: ia
    real(rp), dimension(3), intent(in) :: pos
    real(rp), intent(in) :: radius
    real(rp), dimension(3) :: lbnd, ubnd, extent
    integer :: p

    !If atom already exists, nothing to do.
    if (this%atnd_tab(ia) > 0) return

    !Get a new node from the free store
    p = this%fs_acquire()

    !Prepare the new node
    lbnd = pos - radius; ubnd = pos + radius
    extent = ubnd - lbnd
    !Adjust bounds for a fattened AABB.
    lbnd = lbnd - this%tskin*extent; ubnd = ubnd + this%tskin*extent
    !Initialize AABB with bounds
    call this%nodes(p)%aabb%init(lbnd, ubnd)
    !Set atom
    this%nodes(p)%atom = ia
    !Set node height to zero as this will be a leaf node
    this%nodes(p)%height = 0
    !Insert node as a leaf
    call this%insert_leaf(p)

    !Update atnd_tab
    this%atnd_tab(ia) = p 

    end subroutine

!*******************************************************************************

module subroutine remove(this, ia)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: ia
    integer :: p

    !Get the leaf containing atom ia.
    p = this%atnd_tab(ia)
    !Mark the slot in atnd_tab as vacant.
    this%atnd_tab(ia) = 0

    @: ASSERT(p <= this%capacity)
    @: ASSERT( this%nodes(p)%isleaf() )

    !Remove the leaf from the tree.
    call this%remove_leaf(p)
    !Return the removed node to the free store.
    call this%fs_return(p)

    end subroutine

!*******************************************************************************

module subroutine remove_all(this)

    class(aabbtree_t), intent(in out) :: this
    integer :: ia

    do ia = 1, size(this%atnd_tab)
        if (this%atnd_tab(ia) > 0) call this%remove(ia)
    end do

    end subroutine

!*******************************************************************************

module subroutine update_fatm(this, ia, pos, radius, lstat)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: ia
    real(rp), dimension(3), intent(in) :: pos
    real(rp), intent(in) :: radius
    logical, intent(out) :: lstat
    real(rp), dimension(3) :: lbnd, ubnd

    lbnd = pos - radius; ubnd = pos + radius
    call this%update(ia, lbnd, ubnd, lstat)

    end subroutine

!*******************************************************************************

module subroutine update_fatmaabb(this, ia, lbnd, ubnd, lstat)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: ia
    real(rp), dimension(3), intent(in) :: lbnd
    real(rp), dimension(3), intent(in) :: ubnd
    logical, intent(out) :: lstat
    type(aabb_t) :: aabb
    integer :: p

    lstat = .false.
    !Get the leaf containing atom ia.
    p = this%atnd_tab(ia)
    !Initialize a new AABB
    call aabb%init(lbnd, ubnd)
    !No need to update if the atom is still within its fattened AABB
    if ( this%nodes(p)%aabb%includes(aabb) ) then
        return
    else
        !Remove the current leaf
        call this%remove_leaf(p)
        !Fatten the new AABB
        call aabb%fatten(this%tskin)
        !Assign the new AABB
        this%nodes(p)%aabb = aabb
        !Insert a new leaf node
        call this%insert_leaf(p)
        lstat = .true.
    end if

    end subroutine

!*******************************************************************************

module subroutine query_watm(this, ia, nbrs)

    class(aabbtree_t), intent(in) :: this
    integer, intent(in) :: ia
    type(ivector_t), intent(in out) :: nbrs
    integer :: p

    !Test overlap for atom ia against all other atoms
    p = this%atnd_tab(ia)
    call this%query( ia, this%nodes(p)%aabb, nbrs )

    end subroutine

!*******************************************************************************

module subroutine query_watmaabb(this, ia, aabb, nbrs)

    class(aabbtree_t), intent(in) :: this
    integer, intent(in) :: ia
    type(aabb_t), intent(in) :: aabb
    type(ivector_t), intent(in out) :: nbrs
    integer :: p, ja
    
    call stack%clear()
    call nbrs%clear()

    !First push the root node onto a stack
    call stack%append(this%root)  
    do
        if (stack%len == 0) exit
        p = stack%pop()
        if (p == NULL_NODE) cycle
        !Test for overlap between AABBs
        if ( aabb%overlaps(this%nodes(p)%aabb) ) then
            !Check if this is a leaf node
            if ( this%nodes(p)%isleaf() ) then
                !Can't interact with itself
                ja = this%nodes(p)%atom
                if ( ja /= ia ) call nbrs%append(ja)
            else
                call stack%append(this%nodes(p)%left)
                call stack%append(this%nodes(p)%right)
            end if
        end if
    end do

    end subroutine

!*******************************************************************************

module subroutine query_waabb(this, aabb, atms)

    class(aabbtree_t), intent(in) :: this
    type(aabb_t), intent(in) :: aabb
    type(ivector_t), intent(in out) :: atms

    !If tree is empty return an empty vector
    if (this%size == 0) then
        call atms%clear(); return
    else
        !Test overlap of AABB against all leaf AABBs
        call this%query( huge(0), aabb, atms )
    end if

    end subroutine

!*******************************************************************************

module function get_num_atoms(this) result(na)

    class(aabbtree_t), intent(in) :: this
    integer :: na

    na = count(this%atnd_tab > 0)

    end function

!*******************************************************************************

module function get_aabb(this, ia) result(aabb)

    class(aabbtree_t), intent(in) :: this
    integer, intent(in) :: ia
    type(aabb_t) :: aabb
    integer :: p

    p = this%atnd_tab(ia)
    aabb = this%nodes(p)%aabb

    end function

!*******************************************************************************

module function get_height(this) result(height)

    class(aabbtree_t), intent(in) :: this
    integer :: height

    if (this%root == NULL_NODE) then
        height = 0
    else
        height = this%nodes(this%root)%height
    end if

    end function

!*******************************************************************************

module function get_max_balance(this) result(max_balance)

    class(aabbtree_t), intent(in) :: this
    integer :: max_balance
    integer :: i, balance, left, right

    max_balance = 0
    !Loop over all nodes (including those in the free store)
    do i = 1, this%capacity
        if ( this%nodes(i)%height <= 1 ) cycle
        left = this%nodes(i)%left; right = this%nodes(i)%right
        balance = abs( this%nodes(left)%height - this%nodes(right)%height )
        max_balance = max(max_balance, balance)
    end do

    end function

!*******************************************************************************

module function get_srfarea_ratio(this) result(saratio)

    class(aabbtree_t), intent(in) :: this
    real(rp) :: saratio
    real(rp) :: area_root, area_tot
    integer :: i

    if (this%root == NULL_NODE) then
        saratio = 0.0_rp; return
    end if
    area_root = this%nodes(this%root)%aabb%srfarea
    area_tot = 0.0_rp

    !Loop over all nodes
    do i = 1, this%capacity
        !Ignore nodes in the free store
        if ( this%nodes(i)%height < 0 ) cycle
        area_tot = area_tot + this%nodes(i)%aabb%srfarea
    end do
    saratio = area_tot/area_root

    end function

!*******************************************************************************

module subroutine rebuild(this)

    class(aabbtree_t), intent(in out) :: this
    type(ivector_t) :: node_indices
    type(aabb_t) :: aabbi, aabbj, aabb
    real(rp) :: cost, cost_min
    integer :: counter, i, j, ind, jnd, imin, jmin
    integer :: indx_left, indx_right, p, hl, hr

    counter = 0; call ivector_init(node_indices, this%size )

    !Loop over all nodes and store the leaf node indices, return the rest to the
    !free store.
    do i = 1, this%capacity
        !Ignore nodes in the free store
        if ( this%nodes(i)%height < 0 ) cycle
        if ( this%nodes(i)%isleaf() ) then
            this%nodes(i)%parent = NULL_NODE
            call node_indices%append(i)
            counter = counter + 1
        else
            call this%fs_return(i)
        end if
    end do

    !Rebuild tree from bottom up
    do 
        if (counter <= 1) exit
        cost_min = huge(0.0_rp); imin = 0; jmin = 0
        do i = 1, counter
            ind = node_indices%get_val(i)
            aabbi = this%nodes(ind)%aabb
            do j = (i+1), counter
                jnd = node_indices%get_val(j)
                aabbj = this%nodes(jnd)%aabb
                aabb = aabbi + aabbj 
                cost = aabb%srfarea
                if (cost < cost_min) then
                    imin = i; jmin = j; cost_min = cost
                end if
            end do
        end do

        indx_left = node_indices%get_val(imin)
        indx_right = node_indices%get_val(jmin)

        hl = this%nodes(indx_left)%height
        hr = this%nodes(indx_right)%height

        p = this%fs_acquire()
        this%nodes(p)%left = indx_left; this%nodes(p)%right = indx_right
        this%nodes(p)%height = 1 + max(hl, hr)

        this%nodes(p)%aabb = this%nodes(indx_left)%aabb + this%nodes(indx_right)%aabb
        this%nodes(p)%parent = NULL_NODE
        
        this%nodes(indx_left)%parent = p
        this%nodes(indx_right)%parent = p

        call node_indices%set_val(jmin, node_indices%get_val(counter)) 
        call node_indices%set_val(imin, p)
        counter = counter - 1
    end do

    this%root = node_indices%get_val(1)
    call this%validate()

    end subroutine

!*******************************************************************************

module subroutine validate(this)

    class(aabbtree_t), intent(in) :: this
    integer :: num_free_nodes, p

    call this%validate_structure(this%root)
    call this%validate_metrics(this%root)

    num_free_nodes = 0; p = this%freestore
    do
        if (p == NULL_NODE) exit
        p = this%nodes(p)%next
        num_free_nodes = num_free_nodes + 1
    end do

    @: ASSERT( this%get_height() == this%calc_height() )
    @: ASSERT( (this%size+num_free_nodes) == this%capacity )

    end subroutine

!*******************************************************************************

module subroutine insert_leaf(this, leaf)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: leaf
    type(aabb_t) :: aabb, leaf_aabb, combined_aabb
    real(rp) :: sa, combined_sa, old_area, new_area
    real(rp) :: cost, cost_inheritance, cost_left, cost_right
    integer :: p, sibling, old_parent, new_parent
    integer :: left, right

    !If the tree is empty, insert and make it root
    if (this%root == NULL_NODE) then
        this%root = leaf
        this%nodes(this%root)%parent = NULL_NODE
        return
    end if

    !If the tree is not empty, find the best sibling for the node
    leaf_aabb = this%nodes(leaf)%aabb
    p = this%root
    do
        if ( this%nodes(p)%isleaf() ) exit
        left = this%nodes(p)%left; right = this%nodes(p)%right
        sa = this%nodes(p)%aabb%srfarea

        combined_aabb = this%nodes(p)%aabb + leaf_aabb
        combined_sa = combined_aabb%srfarea

        !Cost of creating a new parent for this node & the new leaf
        cost = 2.0_rp*combined_sa

        !Minimum cost of pushing the leaf further down the tree
        cost_inheritance = 2.0_rp*(combined_sa - sa)

        !Cost of descending to the left
        if ( this%nodes(left)%isleaf() ) then
            aabb = this%nodes(left)%aabb + leaf_aabb
            cost_left = aabb%srfarea + cost_inheritance
        else
            aabb = this%nodes(left)%aabb + leaf_aabb
            old_area = this%nodes(left)%aabb%srfarea
            new_area = aabb%srfarea
            cost_left = (new_area - old_area) + cost_inheritance
        end if

        !Cost of descending to the right
        if ( this%nodes(right)%isleaf() ) then
            aabb = this%nodes(right)%aabb + leaf_aabb
            cost_right = aabb%srfarea + cost_inheritance
        else
            aabb = this%nodes(right)%aabb + leaf_aabb
            old_area = this%nodes(right)%aabb%srfarea
            new_area = aabb%srfarea
            cost_right = (new_area - old_area) + cost_inheritance
        end if

        !Descend according to the minimum cost
        if ( (cost < cost_left) .and. (cost < cost_right) ) exit
        !Descend
        if (cost_left < cost_right) then
            p = left
        else
            p = right
        end if
    end do

    sibling = p
    !Create a new parent
    old_parent = this%nodes(sibling)%parent
    new_parent = this%fs_acquire()
    this%nodes(new_parent)%parent = old_parent
    this%nodes(new_parent)%aabb = this%nodes(sibling)%aabb + leaf_aabb
    this%nodes(new_parent)%height = this%nodes(sibling)%height + 1

    !Sibling was not root
    if (old_parent /= NULL_NODE) then
        if ( this%nodes(old_parent)%left == sibling ) then
            this%nodes(old_parent)%left = new_parent
        else
            this%nodes(old_parent)%right = new_parent
        end if
        this%nodes(new_parent)%left = sibling
        this%nodes(new_parent)%right = leaf
        this%nodes(sibling)%parent = new_parent
        this%nodes(leaf)%parent = new_parent
    else
        !Sibling was the root
        this%nodes(new_parent)%left = sibling
        this%nodes(new_parent)%right = leaf
        this%nodes(sibling)%parent = new_parent
        this%nodes(leaf)%parent = new_parent
        this%root = new_parent
    end if

    !Walk back up the tree fixing heights and AABBs.
    p = this%nodes(leaf)%parent
    do
        if (p == NULL_NODE) exit
        p = this%balance(p)
        left = this%nodes(p)%left; right = this%nodes(p)%right
        this%nodes(p)%height = 1 + &
            max(this%nodes(left)%height, this%nodes(right)%height)
        this%nodes(p)%aabb = this%nodes(left)%aabb + this%nodes(right)%aabb
        p = this%nodes(p)%parent
    end do

    end subroutine

!*******************************************************************************

module subroutine remove_leaf(this, leaf)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: leaf
    integer :: parent, grandparent, sibling
    integer :: p, left, right

    if (leaf == this%root) then
        this%root = NULL_NODE
        return
    end if

    parent = this%nodes(leaf)%parent
    grandparent = this%nodes(parent)%parent

    if ( this%nodes(parent)%left == leaf ) then
        sibling = this%nodes(parent)%right
    else
        sibling = this%nodes(parent)%left
    end if

    !Destroy the parent & connect the sibling to the grandparent
    if (grandparent /= NULL_NODE) then
        if ( this%nodes(grandparent)%left == parent ) then
            this%nodes(grandparent)%left = sibling
        else
            this%nodes(grandparent)%right = sibling
        end if

        this%nodes(sibling)%parent = grandparent
        call this%fs_return(parent)

        !Adjust ancestor bounds
        p = grandparent
        do
            if (p == NULL_NODE) exit
            p = this%balance(p)
            left = this%nodes(p)%left; right = this%nodes(p)%right
            this%nodes(p)%aabb = this%nodes(left)%aabb + this%nodes(right)%aabb
            this%nodes(p)%height = 1 + &
                max(this%nodes(left)%height, this%nodes(right)%height)
            p = this%nodes(p)%parent
        end do
    else
        this%root = sibling
        this%nodes(sibling)%parent = NULL_NODE
        call this%fs_return(parent)
    end if

    end subroutine

!*******************************************************************************

module function balance(this, p) result(q)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: p
    integer :: q
    integer :: left, right, right_left, right_right, left_right, left_left
    integer :: current_balance

    associate (nodes => this%nodes)
    @: ASSERT(p /= NULL_NODE)
    if ( nodes(p)%isleaf() .or. (nodes(p)%height < 2) ) then
        q = p
        return
    end if

    left = nodes(p)%left
    right = nodes(p)%right
    @: ASSERT( left <= this%capacity )
    @: ASSERT( right <= this%capacity )

    current_balance = nodes(right)%height - nodes(left)%height

    !Rotate right branch up
    if (current_balance > 1) then
        right_left = nodes(right)%left
        right_right = nodes(right)%right

        @: ASSERT(right_left  <= this%capacity)
        @: ASSERT(right_right <= this%capacity)

        !Swap node and its right-hand child
        nodes(right)%left = p
        nodes(right)%parent = nodes(p)%parent
        nodes(p)%parent = right

        !The node's old parent should now point to its right-hand child
        if (nodes(right)%parent /= NULL_NODE) then
            if (nodes(nodes(right)%parent)%left == p) then
                nodes(nodes(right)%parent)%left = right
            else
                @: ASSERT(nodes(nodes(right)%parent)%right == p)
                nodes(nodes(right)%parent)%right = right
            end if
        else
            this%root = right
        end if

        !Rotate
        if ( nodes(right_left)%height > nodes(right_right)%height ) then
            nodes(right)%right = right_left
            nodes(p)%right = right_right
            nodes(right_right)%parent = p

            nodes(p)%aabb = nodes(left)%aabb + nodes(right_right)%aabb
            nodes(right)%aabb = nodes(p)%aabb + nodes(right_left)%aabb

            nodes(p)%height = 1 + max(nodes(left)%height, nodes(right_right)%height)
            nodes(right)%height = 1 + max(nodes(p)%height, nodes(right_left)%height)
        else
            nodes(right)%right = right_right
            nodes(p)%right = right_left
            nodes(right_left)%parent = p
            nodes(p)%aabb = nodes(left)%aabb + nodes(right_left)%aabb
            nodes(right)%aabb = nodes(p)%aabb + nodes(right_right)%aabb

            nodes(p)%height = 1 + max(nodes(left)%height, nodes(right_left)%height)
            nodes(right)%height = 1 + max(nodes(p)%height, nodes(right_right)%height)
        end if

        q = right
        return
    end if

    !Rotate left branch up
    if (current_balance < -1) then
        left_left = nodes(left)%left
        left_right = nodes(left)%right
        @: ASSERT(left_left  <= this%capacity)
        @: ASSERT(left_right <= this%capacity)

        !Swap node and its left-hand child
        nodes(left)%left = p
        nodes(left)%parent = nodes(p)%parent
        nodes(p)%parent = left

        !The node's old parent should now point to its left-hand child
        if (nodes(left)%parent /= NULL_NODE) then
            if (nodes(nodes(left)%parent)%left == p) then
                nodes(nodes(left)%parent)%left = left
            else
                @: ASSERT(nodes(nodes(left)%parent)%right == p)
                nodes(nodes(left)%parent)%right = left
            end if
        else
            this%root = left
        end if

        !Rotate
        if ( nodes(left_left)%height > nodes(left_right)%height ) then
            nodes(left)%right = left_left
            nodes(p)%left = left_right
            nodes(left_right)%parent = p

            nodes(p)%aabb = nodes(right)%aabb + nodes(left_right)%aabb
            nodes(left)%aabb = nodes(p)%aabb + nodes(left_left)%aabb

            nodes(p)%height = 1 + max(nodes(right)%height, nodes(left_right)%height)
            nodes(left)%height = 1 + max(nodes(p)%height, nodes(left_left)%height)
        else
            nodes(left)%right = left_right
            nodes(p)%left = left_left
            nodes(left_left)%parent = p
            nodes(p)%aabb = nodes(right)%aabb + nodes(left_left)%aabb
            nodes(left)%aabb = nodes(p)%aabb + nodes(left_right)%aabb

            nodes(p)%height = 1 + max(nodes(right)%height, nodes(left_left)%height)
            nodes(left)%height = 1 + max(nodes(p)%height, nodes(left_right)%height)
        end if

        q = left
        return
    end if

    q = p
    end associate

    end function

!*******************************************************************************

module recursive function calc_height(this, p) result(height)

    class(aabbtree_t), intent(in) :: this
    integer, intent(in), optional :: p
    integer :: height, p_

    if (present(p)) then
        p_ = p
    else
        p_ = this%root
    end if

    if (this%nodes(p_)%isleaf()) then
        height = 0
    else
        height = 1 + max(this%calc_height(this%nodes(p_)%left), &
            this%calc_height(this%nodes(p_)%right))
    end if

    end function

!*******************************************************************************

module function fs_acquire(this) result(p)

    class(aabbtree_t), intent(in out) :: this
    integer :: p
    type(node_t), dimension(:), allocatable :: tmp
    integer :: n, i

    if (this%freestore == NULL_NODE) then
        @: ASSERT(this%size == this%capacity)
        !Double the capacity
        n = this%capacity
        allocate( tmp(2*n) )
        tmp(1:n) = this%nodes
        do i = (n+1), 2*n
            call tmp(i)%init()
            tmp(i)%next = i + 1
        end do
        tmp(n)%next = n + 1
        tmp(2*n)%next = NULL_NODE
        call move_alloc(tmp, this%nodes)
        this%capacity = 2*n
        this%freestore = n + 1
    end if

    !Acquire a new node
    p = this%freestore
    this%freestore = this%nodes(p)%next
    !Increment the size of the tree
    this%size = this%size + 1

    end function

!*******************************************************************************

module subroutine fs_return(this, p)

    class(aabbtree_t), intent(in out) :: this
    integer, intent(in) :: p

    @: ASSERT(p <= this%capacity)
    @: ASSERT(0 < this%size)

    call this%nodes(p)%init()
    this%nodes(p)%next = this%freestore
    this%freestore = p
    !Decrement the size of the tree
    this%size = this%size - 1

    end subroutine

!*******************************************************************************

module recursive subroutine validate_structure(this, p)

    class(aabbtree_t), intent(in) :: this
    integer, intent(in) :: p
    integer :: left, right

    if (p == NULL_NODE) return
    if (p == this%root) then
        @: ASSERT(this%nodes(p)%parent == NULL_NODE)
    end if

    left = this%nodes(p)%left
    right = this%nodes(p)%right

    if ( this%nodes(p)%isleaf() ) then
        @: ASSERT(left == NULL_NODE)
        @: ASSERT(right == NULL_NODE)
        @: ASSERT(this%nodes(p)%height == 0)
        return
    end if

    @: ASSERT( left <= this%capacity )
    @: ASSERT( right <= this%capacity )

    @: ASSERT( this%nodes(left)%parent == p )
    @: ASSERT( this%nodes(right)%parent == p )

    call this%validate_structure(left)
    call this%validate_structure(right) 

    end subroutine

!*******************************************************************************

module recursive subroutine validate_metrics(this, p)

    class(aabbtree_t), intent(in) :: this
    integer, intent(in) :: p
    type(aabb_t) :: aabb
    integer :: left, right, height

    if (p == NULL_NODE) return

    left = this%nodes(p)%left
    right = this%nodes(p)%right

    if ( this%nodes(p)%isleaf() ) then
        @: ASSERT(left == NULL_NODE)
        @: ASSERT(right == NULL_NODE)
        @: ASSERT(this%nodes(p)%height == 0)
        return
    end if

    @: ASSERT( left <= this%capacity )
    @: ASSERT( right <= this%capacity )

    height = 1 + max(this%nodes(left)%height, this%nodes(right)%height)
    @: ASSERT(this%nodes(p)%height == height)

    aabb = this%nodes(left)%aabb + this%nodes(right)%aabb
    @: ASSERT( all(aabb%lbnd == this%nodes(p)%aabb%lbnd) )
    @: ASSERT( all(aabb%ubnd == this%nodes(p)%aabb%ubnd) )

    call this%validate_metrics(left)
    call this%validate_metrics(right) 

    end subroutine

!*******************************************************************************

end submodule aabbtree_sm
