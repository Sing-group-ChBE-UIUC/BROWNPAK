!********************************************************************************!
!                                                                                !
! The MIT License (MIT)                                                          !
!                                                                                !
! Copyright (c) 2020 Sarit Dutta <saridut@gmail.com>                             !
!                                                                                !
! Permission is hereby granted, free of charge, to any person obtaining a copy   !
! of this software and associated documentation files (the "Software"), to deal  !
! in the Software without restriction, including without limitation the rights   !
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      !
! copies of the Software, and to permit persons to whom the Software is          !
! furnished to do so, subject to the following conditions:                       !
!                                                                                !
! The above copyright notice and this permission notice shall be included in all !
! copies or substantial portions of the Software.                                !
!                                                                                !
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     !
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       !
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    !
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         !
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  !
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  !
! SOFTWARE.                                                                      !
!                                                                                !
!********************************************************************************!

module trajectory_m
    !! Routines for reading and writing frames from a trajectory file.
 
use constants_m

implicit none

private
public :: traj

type trajectory_t
    integer :: header_size = 0
    integer :: frame_size = 0
    integer :: num_atoms = 0
    integer, dimension(4) :: frmcmp = 0
    character(len=:), allocatable :: fn
    integer :: file_id = 0
    integer :: num_frames = 0
    logical :: isopen = .false.

    contains
        procedure :: create       => traj_create
        procedure :: open         => traj_open
        procedure :: delete       => traj_delete
        procedure :: close        => traj_close
        procedure :: read         => traj_read
        procedure :: append_frame => traj_append_frame
        procedure :: write_frame  => traj_write_frame
        generic   :: init         => create, open
end type trajectory_t
 
type(trajectory_t) :: traj

contains

!******************************************************************************

subroutine traj_create(this, fn, na, frmcmp)
    !!  Creates a `trajectory_t` object with a new underlying file named `fn`.  If
    !!  `fn` already exists, it will be truncated.  The file `fn` is opened for both
    !!  reading and writing.

    class(trajectory_t), intent(out) :: this
        !! *trajectory_t* instance.
    character(len=*), intent(in) :: fn
        !! Name of the underlying trajectory file.
    integer, intent(in) :: na
        !! Number of atoms
    integer, dimension(4), intent(in) :: frmcmp
        !! Binary flags indicating whether that component is present in a frame.
        !! frmcmp(1): coordinates, frmcmp(2): velocities, frmcmp(3): forces,
        !! frmcmp(4): charges
    integer :: file_id

    this%num_atoms = na
    this%frmcmp = frmcmp
    
    !Frame components: nts, coordinates, velocities, forces, charge
    !nts data
    this%frame_size = sizeof_long_int
    !Coordinate data
    if ( frmcmp(1) /= 0 ) this%frame_size = this%frame_size + 3*na*sizeof_real
    !Velocity data
    if ( frmcmp(2) /= 0 ) this%frame_size = this%frame_size + 3*na*sizeof_real
    !Force data
    if ( frmcmp(3) /= 0 ) this%frame_size = this%frame_size + 3*na*sizeof_real
    !Charge data
    if ( frmcmp(4) /= 0 ) this%frame_size = this%frame_size + na*sizeof_real

    ! Representation of header:
    !     * 1 int : `header_size`
    !     * 1 int : `frame_size`
    !     * 1 int : `num_atoms`
    !     * 4 ints: `frmcmp`

    this%header_size = 7*sizeof_int

    !Create trajectory file
    open(newunit=file_id, file=fn, access='stream', form='unformatted', &
            action='readwrite', status='replace')

    this%fn = fn
    this%file_id = file_id
    this%num_frames = 0
    this%isopen = .true.

    write(this%file_id) this%header_size
    write(this%file_id) this%frame_size
    write(this%file_id) this%num_atoms, this%frmcmp

    end subroutine

!******************************************************************************

subroutine traj_open(this, fn, mode, ierr)
    !!  Creates a `trajectory_t` object with a prexisting underlying file named
    !!  `fn`.  If `fn` does not exist, an error will be generated. If
    !!  `mode` == 'rw', the file `fn` is opened for both reading and writing.
    !!  If `mode` == 'r', the file `fn` is opened only for reading.
    !""
    class(trajectory_t), intent(out) :: this
    character(len=*),    intent(in) :: fn
    character(len=*),    intent(in) :: mode
    integer, intent(out) :: ierr
    integer(ip_long) :: file_size
    
    ierr = 0
    this%fn = fn

    !Readwrite mode
    if (mode == 'rw') then
        open(newunit=this%file_id, file=this%fn, access='stream', &
            form='unformatted', action='readwrite', status='old')
    !Readonly mode
    else if (mode == 'r') then
        open(newunit=this%file_id, file=this%fn, access='stream', &
            form='unformatted', action='read', status='old')
    !Unknown mode
    else
        ierr = 1; return
    end if

    this%isopen = .true.

    !Get size of the file
    inquire(unit=this%file_id, size=file_size)

    read(this%file_id) this%header_size
    read(this%file_id) this%frame_size
    read(this%file_id) this%num_atoms, this%frmcmp

    !Integer division to find `num_frames`
    this%num_frames = int((file_size-this%header_size)/this%frame_size, ip) 

    end subroutine

!******************************************************************************

subroutine traj_delete(this)
   !! After a call to this subroutine, all memory within `this` is deallocated,
   !! all components of `this` are reset to zero, and the underlying file is
   !! closed (if open).

    class(trajectory_t), intent(in out) :: this

    call this%close()

    if (allocated(this%fn)) deallocate(this%fn)

    this%header_size = 0
    this%frame_size = 0
    this%num_atoms = 0
    this%frmcmp = 0
    this%file_id = 0
    this%num_frames = 0

    end subroutine

!******************************************************************************

subroutine traj_close(this)
    !! Closes the underlying file of a `trajectory_t`.

    class(trajectory_t), intent(in out) :: this

    if (this%isopen) then
        close(this%file_id)
        this%isopen = .false.
    end if

    end subroutine

!******************************************************************************

subroutine traj_read(this, iframe, nts, ierr, coordinates, velocities, &
        forces, charge)
    !! Read from an open trajectory.

    class(trajectory_t), intent(in) :: this
    integer, intent(in)  :: iframe
        !! Frame number
    integer(ip_long), intent(out) :: nts
        !! Time step counter
    integer, intent(out) :: ierr
        !! Error flag
    real(rp), dimension(:,:), intent(out), optional :: coordinates
    real(rp), dimension(:,:), intent(out), optional :: velocities
    real(rp), dimension(:,:), intent(out), optional :: forces
    real(rp), dimension(:), intent(out), optional :: charge
    integer(ip_long) :: offset_frm
    integer(ip_long) :: offset_if
    integer(ip_long) :: offset_tot
    integer :: na

    ierr = 0
    na = this%num_atoms

    if (iframe > this%num_frames) then
        write(*,*) '`iframe`', iframe, 'out of bounds for `num_frames`', this%num_frames
        ierr = 1
        return
    end if

    if (na == 0) then
        write(*,*) 'No data'; ierr = 1; return
    end if
    offset_frm = this%frame_size
    offset_frm = (iframe-1)*offset_frm + this%header_size + 1

    read(this%file_id, pos=offset_frm) nts
    if (present(coordinates)) then
        if (this%frmcmp(1) /= 1) then
            write(*,*) 'No coordinate data in frame'; ierr = 1; return
        else
            offset_if = sizeof_long_int
            offset_tot = offset_frm + offset_if
            read(this%file_id, pos=offset_tot) coordinates(:,1:na)
        end if
    end if

    if (present(velocities)) then
        if (this%frmcmp(2) /= 1) then
            write(*,*) 'No velocity data in frame'; ierr = 1; return
        else
            offset_if = sizeof_long_int + 3*na*sizeof_real*this%frmcmp(1)
            offset_tot = offset_frm + offset_if
            read(this%file_id, pos=offset_tot) velocities(:,1:na)
        end if
    end if

    if (present(forces)) then
        if (this%frmcmp(3) /= 1) then
            write(*,*) 'No force data in frame'; ierr = 1; return
        else
            offset_if = sizeof_long_int + 3*na*sizeof_real*this%frmcmp(1) &
                        + 3*na*sizeof_real*this%frmcmp(2)
            offset_tot = offset_frm + offset_if
            read(this%file_id, pos=offset_tot) forces(:,1:na)
        end if
    end if

    if (present(charge)) then
        if (this%frmcmp(4) /= 1) then
            write(*,*) 'No charge data in frame'; ierr = 1; return
        else
            offset_if = sizeof_long_int + 3*na*sizeof_real*this%frmcmp(1) &
                        + 3*na*sizeof_real*this%frmcmp(2) &
                        + 3*na*sizeof_real*this%frmcmp(3)
            offset_tot = offset_frm + offset_if
            read(this%file_id, pos=offset_tot) charge(1:na)
        end if
    end if

    end subroutine

!******************************************************************************

subroutine traj_append_frame(this, nts, coordinates)
    !! Write a frame to an open trajectory

    class(trajectory_t),  intent(in out) :: this
    integer(ip_long),         intent(in) :: nts
    real(rp), dimension(:,:), intent(in) :: coordinates
    integer :: iframe

    iframe = this%num_frames + 1
    call this%write_frame(iframe, nts, coordinates)

    end subroutine

!******************************************************************************

subroutine traj_write_frame(this, iframe, nts, coordinates)
    !! Write a frame to an open trajectory

    class(trajectory_t),  intent(in out) :: this
    integer,                  intent(in) :: iframe
    integer(ip_long),         intent(in) :: nts
    real(rp), dimension(:,:), intent(in) :: coordinates
    integer(ip_long) :: offset
    integer :: na

    na = this%num_atoms

    ! Added one to handle appending a frame
    if (iframe > this%num_frames+1) then
        write(*,*) '`iframe`', iframe, 'out of bounds for `num_frames`', &
            (this%num_frames+1)
        return
    end if

    offset = this%frame_size
    offset = (iframe-1)*offset + this%header_size + 1

    write(this%file_id, pos=offset) nts

    if (this%frmcmp(1) == 1) then
        if (this%num_atoms /= 0) then
            write(this%file_id) coordinates(:,1:na)
        end if
    end if

    if (iframe == this%num_frames + 1) this%num_frames = this%num_frames + 1

    end subroutine

!******************************************************************************

end module trajectory_m
