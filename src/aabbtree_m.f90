! This code is adapted from the C++ implementation by 
! [Lester Hedges](http://lesterhedges.net). The original
! C++ implementation can be found [here](https://github.com/lohedges/aabbcc#readme).
!
 !****************************************************************************!
 !                                                                            !
 ! The zlib/libpng License (Zlib)                                             !
 !                                                                            !
 ! Copyright (c) 2009 Erin Catto http://www.box2d.org                         !
 ! Copyright (c) 2016 Lester Hedges <lester.hedges+aabbcc@gmail.com>          !
 !                                                                            !
 ! This software is provided 'as-is', without any express or implied          !
 ! warranty. In no event will the authors be held liable for any damages      !
 ! arising from the use of this software.                                     !
 !                                                                            !
 ! Permission is granted to anyone to use this software for any purpose,      !
 ! including commercial applications, and to alter it and redistribute it     !
 ! freely, subject to the following restrictions:                             !
 !                                                                            !
 ! 1. The origin of this software must not be misrepresented; you must not    !
 !    claim that you wrote the original software. If you use this software    !
 !    in a product, an acknowledgment in the product documentation would be   !
 !    appreciated but is not required.                                        !
 !                                                                            !
 ! 2. Altered source versions must be plainly marked as such, and must not be !
 !    misrepresented as being the original software.                          !
 !                                                                            !
 ! 3. This notice may not be removed or altered from any source distribution. !
 !                                                                            !
 ! This code was adapted from parts of the Box2D Physics Engine,              !
 ! http://www.box2d.org                                                       !
 !                                                                            !
 !****************************************************************************!

module aabbtree_m
!! Implements an axis-aligned bounding box (AABB) tree.
!! This code is adapted from the C++ implementation by 
!! [Lester Hedges](http://lesterhedges.net). The original
!! C++ implementation can be found [here](https://github.com/lohedges/aabbcc#readme).

use constants_m
use strings_m
use vector_m
use aabb_m

implicit none

private

public :: aabbtree_t

integer, parameter :: NULL_NODE = 0

type node_t
    integer :: next = NULL_NODE
    integer :: parent = NULL_NODE
    integer :: left = NULL_NODE
    integer :: right = NULL_NODE
    integer :: height = -1
    integer :: atom = 0
    type(aabb_t) :: aabb
    contains
        procedure :: init => node_init
        procedure :: isleaf => node_isleaf
        procedure :: asstr => node_asstr
end type node_t

type(ivector_t) :: stack

type aabbtree_t
    type(node_t), dimension(:), allocatable :: nodes
    integer, dimension(:), allocatable :: atnd_tab
        !! Atom -> node map.
    integer :: capacity
        !! Maximum number of nodes that the tree can currently handle. This may
        !! increase as more atoms are inserted.
    integer :: freestore
        !! Pointer to head to the free store.
    integer :: size
        !! Number of nodes in the tree.
    integer :: root
        !! Pointer to the tree root.
    real(rp) :: tskin
        !! Thickness of the skin for fattened AABBs, as a fraction of the AABB
        !! base length.
    contains
        procedure :: init
        procedure :: print
        procedure :: clear
        procedure :: delete
        procedure :: insert
        procedure :: remove
        procedure :: remove_all
        procedure :: update_fatm
        procedure :: update_fatmaabb
        generic   :: update => update_fatm, update_fatmaabb
        procedure :: query_watm
        procedure :: query_waabb
        procedure :: query_watmaabb
        generic   :: query => query_watm, query_waabb, query_watmaabb
        procedure :: get_num_atoms
        procedure :: get_aabb
        procedure :: get_height
        procedure :: get_max_balance
        procedure :: get_srfarea_ratio
        procedure :: rebuild
        procedure :: validate
        procedure, private :: insert_leaf
        procedure, private :: remove_leaf
        procedure, private :: balance
        procedure, private :: calc_height
        procedure, private :: fs_acquire
        procedure, private :: fs_return
        procedure, private :: validate_structure
        procedure, private :: validate_metrics
end type aabbtree_t

interface
    module subroutine node_init(this)
        !! Initializes a node.
        class(node_t), intent(out) :: this
            !! A *node_t* instance.
        end subroutine

    module function node_isleaf(this) result(res)
        !! Is this a leaf node? 
        class(node_t), intent(in) :: this
            !! A *node_t* instance.
        logical :: res
            !! *true* if this is a leaf node, *false* otherwise.
        end function

    module function node_asstr(this, frmt) result(buf)
        !! Returns a string representation of a node.
        class(node_t), intent(in) :: this
            !! A *node_t* instance.
        character(len=*), intent(in), optional :: frmt
            !! Format string for real numbers, e.g. '(f15.6)'. Default: *(g0.6)*.
        character(len=:), allocatable :: buf
            !! Return value
        end function

    module subroutine init(this, natoms, tskin)
        !! Initializes an AABB tree.
        class(aabbtree_t), intent(out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: natoms
            !! Estimated number of atoms to be handled by this tree.
        real(rp), intent(in) :: tskin
            !! Thickness of the skin for fattened AABBs, as a fraction of the
            !! AABB base length.
        end subroutine

    module recursive subroutine print(this, p)
        !! Prints a subtree of an AABB tree rooted at `p` in order.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in), optional :: p
            !! Pointer to the root of the subtree. Default is the root of the
            !! whole tree.
        end subroutine

    module subroutine clear(this)
        !! Clears an AABB tree. Associated memory is not deallocated.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        end subroutine

    module subroutine delete(this)
        !! Deletes an AABB tree. All associated memory is deallocated.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        end subroutine

    module subroutine insert(this, ia, pos, radius)
        !! Inserts an atom into an AABB tree.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: ia
            !! Atom index
        real(rp), dimension(3), intent(in) :: pos
            !! Atom position
        real(rp), intent(in) :: radius
            !! Atom radius (or cutoff distance for point particles)
        end subroutine

    module subroutine remove(this, ia)
        !! Removes an atom from an AABB tree.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: ia
            !! Atom index
        end subroutine

    module subroutine remove_all(this)
        !! Removes all atoms from an AABB tree.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        end subroutine

    module subroutine update_fatm(this, ia, pos, radius, lstat)
        !! Updates an AABB tree for the case when an atom leaves its fattened
        !! AABB.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: ia
            !! Atom index
        real(rp), dimension(3), intent(in) :: pos
            !! Atom position
        real(rp), intent(in) :: radius
            !! Atom radius (or cutoff distance for point particles)
        logical, intent(out) :: lstat
            !! *true* if the atom was reinserted, *false* otherwise.
        end subroutine

    module subroutine update_fatmaabb(this, ia, lbnd, ubnd, lstat)
        !! Updates an AABB tree for the case when an atom leaves its fattened
        !! AABB, with the bounds of the new AABB for the atom as input.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: ia
            !! Atom index
        real(rp), dimension(3), intent(in) :: lbnd
            !! Lower bound of atom AABB
        real(rp), dimension(3), intent(in) :: ubnd
            !! Upper bound of atom AABB
        logical, intent(out) :: lstat
            !! *true* if the atom was reinserted, *false* otherwise.
        end subroutine

    module subroutine query_watm(this, ia, nbrs)
        !! Query an AABB tree for a set of potential neighbors of an atom.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: ia
            !! Atom index
        type(ivector_t), intent(in out) :: nbrs
            !! List of potential neighbors.
        end subroutine

    module subroutine query_watmaabb(this, ia, aabb, nbrs)
        !! Query an AABB tree for a set of potential neighbors of an atom & its
        !! bounding AABB.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: ia
            !! Atom index.
        type(aabb_t), intent(in) :: aabb
            !! Bounding AABB for atom with index `ia`.
        type(ivector_t), intent(in out) :: nbrs
            !! List of potential neighbors.
        end subroutine

    module subroutine query_waabb(this, aabb, atms)
        !! Query an AABB tree for the set of atoms whose AABBS overlap with
        !! `aabb`. 
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        type(aabb_t), intent(in) :: aabb
            !! An *aabb_t* instance.
        type(ivector_t), intent(in out) :: atms
            !! List of potential neighbors.
        end subroutine

    module function get_num_atoms(this) result(na)
        !! Returns the number of atoms in an AABB tree.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer :: na
            !! Return value
        end function

    module function get_aabb(this, ia) result(aabb)
        !! Returns a copy of the AABB associated with atom with index `ia`.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: ia
            !! Atom index.
        type(aabb_t) :: aabb
            !! AABB of atom `ia`.
        end function

    module function get_height(this) result(height)
        !! Returns the height of an AABB tree.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer :: height
            !! Return value.
        end function

    module function get_max_balance(this) result(max_balance)
        !! Returns the maximum difference between the height of two children
        !! of a node.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer :: max_balance
            !! Return value.
        end function
        
    module function get_srfarea_ratio(this) result(saratio)
        !! Returns the ratio of the sum of the node surface area to the surface
        !! area of the root node.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        real(rp) :: saratio
            !! Return value.
        end function

    module subroutine rebuild(this)
        !! Rebuilds an optimal AABB tree.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        end subroutine

    module subroutine validate(this)
        !! Validates an AABB tree.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        end subroutine

    module subroutine insert_leaf(this, leaf)
        !! Inserts a leaf node into a tree.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: leaf
            !! Pointer to a leaf node
        end subroutine

    module subroutine remove_leaf(this, leaf)
        !! Removes a leaf node from a tree.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: leaf
            !! Pointer to a leaf node
        end subroutine

    module function balance(this, p) result(q)
        !! Balances an AABB tree.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: p
            !! Node index
        integer :: q

        end function

    module recursive function calc_height(this, p) result(height)
        !! Calculates the height of a subtree.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in), optional :: p
            !! Pointer to the root of the subtree. Default: Root of the whole
            !! tree.
        integer :: height
            !! Return value
        end function

    module function fs_acquire(this) result(p)
        !! Acquires a new node from the free store and returns a pointer to it.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer :: p
            !! Return value
        end function

    module subroutine fs_return(this, p)
        !! Returns a node to the free store.
        class(aabbtree_t), intent(in out) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: p
            !! Node pointer.
        end subroutine

    module recursive subroutine validate_structure(this, p)
        !! Asserts that an AABB subtree has a valid structure.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: p
            !! Pointer to the root of the subtree.
        end subroutine

    module recursive subroutine validate_metrics(this, p)
        !! Asserts that an AABB subtree has a valid metric.
        class(aabbtree_t), intent(in) :: this
            !! An *aabbtree_t* instance.
        integer, intent(in) :: p
            !! Pointer to the root of the subtree.
        end subroutine
end interface

!*******************************************************************************

end module aabbtree_m
