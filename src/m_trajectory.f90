module m_trajectory
    !! Routines for reading and writing frames from a trajectory file.
 
use m_precision

implicit none

type trajectory_t
    integer :: header_size = 0
    integer :: frame_size = 0
    integer :: num_atoms = 0
    integer :: num_mpcd_atoms = 0
    integer :: num_atoms_tot = 0
    integer, dimension(4) :: frmcmp = 0
    character(len=:), allocatable :: fn
    integer :: file_id = 0
    integer :: num_frames = 0
    logical :: isopen = .false.

    contains
        procedure :: create       => traj_create
        procedure :: open         => traj_open
        procedure :: clear        => traj_clear
        procedure :: close        => traj_close
        procedure :: read         => traj_read
        procedure :: append_frame => traj_append_frame
        procedure :: write_frame  => traj_write_frame
        generic   :: init         => create, open
end type trajectory_t
 
contains

!******************************************************************************

subroutine traj_create(this, fn, na, nam, frmcmp)
    !!  Creates a `trajectory_t` object with a new underlying file named `fn`.  If
    !!  `fn` already exists, it will be truncated.  The file `fn` is opened for both
    !!  reading and writing.

    class(trajectory_t), intent(out) :: this
        !! *trajectory_t* instance.
    character(len=*), intent(in) :: fn
        !! Name of the underlying trajectory file.
    integer, intent(in) :: na
        !! Number of non-MPCD atoms
    integer, intent(in) :: nam
        !! Number of MPCD atoms, pass zero to indicate no data from MPCD
        !! particles are present.
    integer, dimension(4), intent(in) :: frmcmp
        !! Binary flags indicating whether that component is present in a frame.
        !! frmcmp(1): coordinates, frmcmp(2): velocities, frmcmp(3): forces,
        !! frmcmp(4): charges
    integer :: nat, file_id

    nat = na + nam
    this%num_atoms = na; this%num_mpcd_atoms = nam; this%num_atoms_tot = nat
    this%frmcmp = frmcmp
    
    !Frame components: nts, coordinates, velocities, forces, charge
    !nts data
    this%frame_size = sizeof_long_int
    !Coordinate data
    if ( frmcmp(1) /= 0 ) this%frame_size = this%frame_size + 3*nat*sizeof_real
    !Velocity data
    if ( frmcmp(2) /= 0 ) this%frame_size = this%frame_size + 3*nat*sizeof_real
    !Force data
    if ( frmcmp(3) /= 0 ) this%frame_size = this%frame_size + 3*na*sizeof_real
    !Charge data
    if ( frmcmp(4) /= 0 ) this%frame_size = this%frame_size + na*sizeof_real

    ! Representation of header:
    !     * 1 int for `header_size`
    !     * 1 int for `frame_size`
    !     * 2 ints:  `num_atoms`, `num_mpcd_atoms`
    !     * 4 ints:  `frmcmp`

    this%header_size = 8*sizeof_int

    !Create trajectory file
    open(newunit=file_id, file=fn, access='stream', form='unformatted', &
            action='readwrite', status='replace')

    this%fn = fn
    this%file_id = file_id
    this%num_frames = 0
    this%isopen = .true.

    write(this%file_id) this%header_size
    write(this%file_id) this%frame_size
    write(this%file_id) this%num_atoms, this%num_mpcd_atoms, this%frmcmp

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
    read(this%file_id) this%num_atoms, this%num_mpcd_atoms, this%frmcmp

    this%num_atoms_tot = this%num_atoms + this%num_mpcd_atoms

    !Integer division to find `num_frames`
    this%num_frames = int((file_size-this%header_size)/this%frame_size, ip) 

    end subroutine

!******************************************************************************

subroutine traj_clear(this)
   !! After a call to this subroutine, all memory within `this` is deallocated,
   !! all components of `this` are reset to zero, and the underlying file is
   !! closed (if open).

    class(trajectory_t), intent(in out) :: this

    call this%close()

    if (allocated(this%fn)) deallocate(this%fn)

    this%header_size = 0
    this%frame_size = 0
    this%num_atoms = 0
    this%num_mpcd_atoms = 0
    this%num_atoms_tot = 0
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

subroutine traj_read(this, iframe, nts, ierr, mflag, &
        coordinates, velocities, forces, charge)
    !! Read from an open trajectory.

    class(trajectory_t), intent(in) :: this
    integer, intent(in)  :: iframe
        !! Frame number
    integer(ip_long), intent(out) :: nts
        !! Time step counter
    integer, intent(out) :: ierr
        !! Error flag
    integer, intent(in), optional  :: mflag
        !! 1: Only non-MPCD data, 2: Only MPCD data, 3: All data
    real(rp), dimension(:,:), intent(out), optional :: coordinates
    real(rp), dimension(:,:), intent(out), optional :: velocities
    real(rp), dimension(:,:), intent(out), optional :: forces
    real(rp), dimension(:), intent(out), optional :: charge
    integer(ip_long) :: offset_frm
    integer(ip_long) :: offset_if
    integer(ip_long) :: offset_tot
    integer :: na, nat, nam

    ierr = 0
    na = this%num_atoms; nam = this%num_mpcd_atoms; nat = this%num_atoms_tot

    if (iframe > this%num_frames) then
        write(*,*) '`iframe`', iframe, 'out of bounds for `num_frames`', this%num_frames
        ierr = 1
        return
    end if

    select case(mflag)
    case(1)
        if (na == 0) then
            write(*,*) 'No data for non-MPCD atoms'; ierr = 1; return
        end if
    case(2)
        if (nam == 0) then
            write(*,*) 'No data for MPCD atoms'; ierr = 1; return
        end if
    case(3)
        if (nat == 0) then
            write(*,*) 'No data for any atoms'; ierr = 1; return
        end if
    case default
        write(*,*) 'Unknown mflag = ', mflag; ierr = 1; return
    end select

    offset_frm = this%frame_size
    offset_frm = (iframe-1)*offset_frm + this%header_size + 1

    read(this%file_id, pos=offset_frm) nts
    if (present(coordinates)) then
        if (this%frmcmp(1) /= 1) then
            write(*,*) 'No coordinate data in frame'; ierr = 1; return
        else
            select case(mflag)
            case(1)
                offset_if = sizeof_long_int
                offset_tot = offset_frm + offset_if
                read(this%file_id, pos=offset_tot) coordinates(:,1:na)
            case(2)
                offset_if = sizeof_long_int + 3*na*sizeof_real
                offset_tot = offset_frm + offset_if
                read(this%file_id, pos=offset_tot) coordinates(:,1:nam)
            case(3)
                offset_if = sizeof_long_int
                offset_tot = offset_frm + offset_if
                read(this%file_id, pos=offset_tot) coordinates(:,1:nat)
            end select
        end if
    end if

    if (present(velocities)) then
        if (this%frmcmp(2) /= 1) then
            write(*,*) 'No velocity data in frame'; ierr = 1; return
        else
            select case(mflag)
            case(1)
                offset_if = sizeof_long_int + 3*nat*sizeof_real*this%frmcmp(1)
                offset_tot = offset_frm + offset_if
                read(this%file_id, pos=offset_tot) velocities(:,1:na)
            case(2)
                offset_if = sizeof_long_int + 3*nat*sizeof_real*this%frmcmp(1) &
                            + 3*na*sizeof_real*this%frmcmp(2)
                offset_tot = offset_frm + offset_if
                read(this%file_id, pos=offset_tot) velocities(:,1:nam)
            case(3)
                offset_if = sizeof_long_int + 3*nat*sizeof_real*this%frmcmp(1)
                offset_tot = offset_frm + offset_if
                read(this%file_id, pos=offset_tot) velocities(:,1:nat)
            end select
        end if
    end if

    if (present(forces)) then
        if (this%frmcmp(3) /= 1) then
            write(*,*) 'No force data in frame'; ierr = 1; return
        else
            offset_if = sizeof_long_int + 3*nat*sizeof_real*this%frmcmp(1) &
                        + 3*nat*sizeof_real*this%frmcmp(2)
            offset_tot = offset_frm + offset_if
            read(this%file_id, pos=offset_tot) forces(:,1:na)
        end if
    end if

    if (present(charge)) then
        if (this%frmcmp(4) /= 1) then
            write(*,*) 'No charge data in frame'; ierr = 1; return
        else
            offset_if = sizeof_long_int + 3*nat*sizeof_real*this%frmcmp(1) &
                        + 3*nat*sizeof_real*this%frmcmp(2) &
                        + 3*na*sizeof_real*this%frmcmp(3)
            offset_tot = offset_frm + offset_if
            read(this%file_id, pos=offset_tot) charge(1:na)
        end if
    end if

    end subroutine

!******************************************************************************

subroutine traj_append_frame(this, nts, coordinates, velocities, forces, charge)
    !! Write a frame to an open trajectory

    class(trajectory_t),  intent(in out) :: this
    integer(ip_long),         intent(in) :: nts
    real(rp), dimension(:,:), intent(in) :: coordinates
    real(rp), dimension(:,:), intent(in) :: velocities
    real(rp), dimension(:,:), intent(in) :: forces
    real(rp), dimension(:), intent(in)   :: charge
    integer :: iframe

    iframe = this%num_frames + 1
    call this%write_frame(iframe, nts, coordinates, velocities, forces, charge)

    end subroutine

!******************************************************************************

subroutine traj_write_frame(this, iframe, nts, coordinates, velocities, forces, charge)
    !! Write a frame to an open trajectory

    class(trajectory_t),  intent(in out) :: this
    integer,                  intent(in) :: iframe
    integer(ip_long),         intent(in) :: nts
    real(rp), dimension(:,:), intent(in) :: coordinates
    real(rp), dimension(:,:), intent(in) :: velocities
    real(rp), dimension(:,:), intent(in) :: forces
    real(rp), dimension(:), intent(in)   :: charge
    integer(ip_long) :: offset
    integer :: na, nat, nam

    na = this%num_atoms; nam = this%num_mpcd_atoms; nat = this%num_atoms_tot

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
        if (this%num_atoms == 0) then
            write(this%file_id) coordinates(:,na+1:)
        else if (this%num_mpcd_atoms == 0) then
            write(this%file_id) coordinates(:,1:na)
        else
            write(this%file_id) coordinates(:,1:nat)
        end if
    end if

    if (this%frmcmp(2) == 1) then
        if (this%num_atoms == 0) then
            write(this%file_id) velocities(:,na+1:)
        else if (this%num_mpcd_atoms == 0) then
            write(this%file_id) velocities(:,1:na)
        else
            write(this%file_id) velocities(:,1:nat)
        end if
    end if

    if (this%frmcmp(3) == 1) write(this%file_id) forces

    if (this%frmcmp(4) == 1) write(this%file_id) charge

    if (iframe == this%num_frames + 1) this%num_frames = this%num_frames + 1

    end subroutine

!******************************************************************************

end module m_trajectory
