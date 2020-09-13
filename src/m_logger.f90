module m_logger
    !! Implements a basic logger.

use iso_fortran_env, only: output_unit

implicit none

private
public :: logger_init, logger 

type logger_t
    character(len=:), allocatable :: fn
        !! Name of the log file
    integer :: fu = huge(0)
        !! Unit number of the log file
    logical :: is_open = .false.
        !! Is the log file open for writing? {T/F}
    contains
        procedure :: finish => logger_finish
        procedure :: log_msg => logger_log_msg
end type logger_t

interface logger_init
    module procedure logger_init
end interface

type(logger_t) :: logger

contains

!********************************************************************************

subroutine logger_init(this, fn, use_stdout)
    !! Initializes a logger.

    type(logger_t), intent(in out) :: this
        !! A `logger_t` instance
    character(len=*), intent(in) :: fn
        !! Name of the log file. If `use_stdout` is true, `fn` is ignored.
    logical, intent(in) :: use_stdout
        !! Write all log messages to stdout rather than a file on disk? {T/F}
    integer :: ierr

    !Create a new log file and open it for writing.
    if (use_stdout) then
        this%fu = output_unit
        this%is_open = .true.
    else
        open(newunit=this%fu, file=fn, action='write', status='replace', &
            iostat=ierr)
        if (ierr /= 0) then
            write(*,'(a,i0)') 'error in opening log file; err code = ', ierr
        else
            this%fn = fn
            this%is_open = .true.
        end if
    end if

    call this%log_msg('start log')

    end subroutine

!********************************************************************************

subroutine logger_finish(this)
    !! Cleanup routine for a `logger_t` instance.

    class(logger_t), intent(in out) :: this
        !! A `logger_t` instance
    
    call this%log_msg('end log')

    if (this%is_open .and. (this%fu /= output_unit)) then
        close(this%fu)
        this%is_open = .false.
        this%fu = huge(0)
    end if
    
    end subroutine

!********************************************************************************

subroutine logger_log_msg(this, msg)
    !! Write a message to the log file.

    class(logger_t), intent(in) :: this
        !! A `logger_t` instance
    character(len=*), intent(in) :: msg
        !! Message to write to the log file
    character(len=40) :: timstmp

    call timestring(timstmp)

    write(this%fu, "('[',a,']',1x,a)") trim(adjustl(timstmp)), &
        trim(adjustl(msg))
    flush(this%fu)

    end subroutine

!********************************************************************************

subroutine timestring ( string )
    !! TIMESTRING writes the current YMDHMS date into a string.
    !
    !  Example:
    !
    !    STRING = '31 May 2001   9:45:54.872 AM'
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, character ( len = * ) STRING, contains the date information.
    !    A character length of 40 should always be sufficient.
    !
    implicit none

    character(len=8) :: ampm
    integer :: d
    integer :: h
    integer :: m
    integer :: mm
    character ( len = 6 ), parameter, dimension(12) :: month = (/ &
      'Jan   ', 'Feb   ', 'March ', 'April ', &
      'May   ', 'June  ', 'July  ', 'Aug   ', &
      'Sept  ', 'Oct   ', 'Nov   ', 'Dec   ' /)
    integer :: n
    integer :: s
    character(len=*) string
    integer :: values(8)
    integer :: y

    call date_and_time ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
      ampm = 'AM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Noon'
      else
        ampm = 'PM'
      end if
    else
      h = h - 12
      if ( h < 12 ) then
        ampm = 'PM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Midnight'
        else
          ampm = 'AM'
        end if
      end if
    end if

    write ( string, '(i2.2,1x,a,1x,i4,1x,i0,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
      d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    return
    end subroutine

!********************************************************************************

end module m_logger
