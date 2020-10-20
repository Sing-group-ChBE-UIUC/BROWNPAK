!********************************************************************************!
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
!********************************************************************************!


module logger_m
    !! Implements a basic logger.

use iso_fortran_env, only: output_unit
use timestamp_m, only: timestring

implicit none

private
public :: logger 

type logger_t
    character(len=:), allocatable :: fn
        !! Name of the log file
    integer :: fu = huge(0)
        !! Unit number of the log file
    logical :: is_open = .false.
        !! Is the log file open for writing? {T/F}
    contains
        procedure :: init => logger_init
        procedure :: finish => logger_finish
        procedure :: log_msg => logger_log_msg
end type logger_t

type(logger_t) :: logger

contains

!********************************************************************************

subroutine logger_init(this, fn, use_stdout)
    !! Initializes a logger.

    class(logger_t), intent(out) :: this
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

end module logger_m
