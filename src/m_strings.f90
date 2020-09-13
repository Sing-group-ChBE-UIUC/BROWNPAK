module m_strings
!! summary: Routines for string handling
!! author:  George Benthien & others
!!
!! Many of these routines were originally written by
!! [George Benthien](https://gbenthien.net/strings/str-index.html), some have
!! been written by Sarit Dutta or gathered/modified from other authors (appropriately credited).
!!
!! These routines were developed primarily to aid in the reading and manipulation
!! of input data from an ASCII text file. Accordingly, it is assumed that all
!! characters to be processed are ASCII characters.  

use m_precision

implicit none

private :: str_from_inum, str_from_ilnum, str_from_dnum

interface str_from_num
!!  Generic  interface for writing a number to a string. The calling syntax is 
!!  `str_from_num(num, frmt)` where `number` is a real number or an integer,
!!  `format` is the format desired, e.g., *e15.6*, *i5*, etc.
    module procedure str_from_inum
    module procedure str_from_ilnum
    module procedure str_from_dnum
end interface str_from_num

contains

!******************************************************************************

pure function str_is_letter(str) result(res)
    !! Returns `.true.` if `str` contains only letters (*a--z* or *A--Z*) and
    !! `.false.` otherwise.

    character(len=*), intent(in) :: str
    logical :: res
    character(len=1) :: ch
    integer :: i
    
    res = .true.

    do i = 1, len(str)
        ch = str(i:i)

        select case(ch)
        case ('A':'Z')
            cycle
        case ('a':'z')
            cycle
        case default
            res = .false.
            return
        end select

    end do
    
    end function

!******************************************************************************

pure function str_is_digit(str) result(res)
    !! Returns `.true.` if `str` contains only digits (0,1,...,9) and 
    !!`.false.` otherwise.

    character(len=*), intent(in) :: str
    logical :: res
    character(len=1) :: ch
    integer :: i

    res = .true.

    do i = 1, len(str)
        ch = str(i:i)

        select case(ch)
        case('0':'9')
            cycle
        case default
            res = .false.
            return
        end select

    end do
    
    end function

!*******************************************************************************

pure function str_is_space(str) result(res)
    !! Returns `.true.` if `str` is non-empty and contains only whitespace
    !! characters (tab or blankspace). Otherwise `.false.` is returned.
    !!
    !! *Note*: This function will return `.false.` for an empty string.

    character(len=*), intent(in) :: str
    logical :: res
    integer :: lenstr
    integer :: ich
    integer :: i

    lenstr = len(str)

    if (lenstr == 0) then
        res = .false.
        return
    end if

    res = .true.

    do i = 1, lenstr
        ich = iachar(str(i:i))
        select case(ich)
        case (9,32)
            res = .true.
        case default
            res = .false.
            return
        end select
    end do

    end function

!*******************************************************************************

pure function str_is_comment(line, comment_str) result(res)
    !!  Returns `.true.` if `line` is a comment, `.false.` other wise.
    !!
    !!  `line` is a comment if `comment_str` is its first non-blank character
    !!  sequence. If `line` is an empty string or contains only blankspaces, the
    !!  return value is `.false.` If `comment_str` is empty, the return value is
    !!  `.true.`.


    character(len=*), intent (in) :: line
        !!  Input string
    character(len=*), intent (in) :: comment_str
        !!  String marking the beginning of a comment.
    logical :: res

    res = .false.

    if ( index(adjustl(line),comment_str) == 1 ) then
         res = .true.
    end if

    end function

!*******************************************************************************

pure function str_compact(str) result(ostr)
    !! Returns a copy of `str` with multiple spaces and tabs converted to
    !! single spaces, control characters deleted, and leading and trailing
    !! spaces removed.
    
    character(len=*), intent(in)  :: str
    character(len=:), allocatable :: ostr
    character(len=len(str)) :: buf
    character(len=1) :: ch
    integer :: lenstr
    integer :: isp, ich
    integer :: i, k
    
    lenstr = len(str)
    isp = 0; k = 0
    buf = ''
    
    do i = 1, lenstr
        ch = str(i:i)
        ich = iachar(ch)
        
        select case(ich)
            case(9,32)     ! space or tab character
                if (isp==0) then
                    k = k+1
                    buf(k:k) = ' '
                end if
                isp = 1
            case(33:)      ! not a space, quote, or control character
                k = k + 1
                buf(k:k) = ch
                isp = 0
        end select
    end do
    
    ostr = trim(adjustl(buf))
    
    end function

!******************************************************************************

pure function str_remove_stcc(str) result(ostr)
    !! Returns a copy of the string `str` with spaces, tabs, and
    !! control characters removed.

    character(len=*), intent(in)  :: str
    character(len=:), allocatable :: ostr
    character(len=len(str)) :: buf
    character(len=1) :: ch
    integer :: lenstr
    integer :: ich
    integer :: k
    integer :: i
    
    lenstr = len(str)
    k = 0
    
    do i = 1, lenstr
        ch = str(i:i)
        ich = iachar(ch)

        select case(ich)    
            case(0:32)  ! space, tab, or control character
                cycle       
            case(33:)  
                k = k+1
                buf(k:k) = ch
        end select
    end do
    
    ostr = trim(adjustl(buf))
    
    end function

!******************************************************************************

pure function str_to_upper(str) result(ucstr)
    !!  This function returns a string that is like the string `str` with all characters
    !!  that are not between a pair of quotes ("..." or '...') converted to uppercase.
    
    character (len=*), intent(in) :: str
    character (len=len(str)) :: ucstr
    integer :: ilen
    integer :: ioffset
    integer :: iquote
    integer :: iav
    integer :: iqc
    integer :: i
    
    ilen = len(str)
    ioffset = iachar('A') - iachar('a')     
    iquote = 0
    ucstr = str

    do i = 1, ilen
        iav = iachar(str(i:i))

        if ( (iquote==0) .and. ((iav==34) .or. (iav==39)) ) then
            iquote = 1
            iqc = iav
            cycle
        end if

        if ( (iquote==1) .and. (iav==iqc) ) then
            iquote = 0
            cycle
        end if

        if (iquote==1) cycle

        if ( (iav >= iachar('a')) .and. (iav <= iachar('z')) ) then
            ucstr(i:i) = achar(iav+ioffset)
        end if
    end do
    
    end function

!**********************************************************************

pure function str_to_lower(str) result(lcstr)
    !!  This function returns a string that is like the string `str` with all characters
    !!  that are not between a pair of quotes ("..." or '...') converted to lowercase.

    character (len=*), intent(in):: str
    character (len=len(str))     :: lcstr
    integer :: ilen
    integer :: ioffset
    integer :: iquote
    integer :: iav
    integer :: iqc
    integer :: i
    
    ilen = len(str)
    ioffset = iachar('A')-iachar('a')
    iquote = 0
    lcstr = str
    
    do i = 1, ilen
        iav = iachar(str(i:i))

        if ( (iquote==0) .and. ((iav==34) .or. (iav==39)) ) then
            iquote = 1
            iqc = iav
            cycle
        end if

        if ( (iquote==1) .and. (iav==iqc) ) then
            iquote=0
            cycle
        end if

        if (iquote==1) cycle

        if ( (iav >= iachar('A')) .and. (iav <= iachar('Z')) ) then
            lcstr(i:i) = achar(iav-ioffset)
        end if
    end do
    
    end function

!******************************************************************************

subroutine str_shift(str, n)
    !! Shifts characters in `str` by `n` positions (positive values
    !! denote a right shift and negative values denote a left shift). Characters
    !! that are shifted off the end are lost. Positions opened up by the shift 
    !! are replaced by spaces.

    character(len=*), intent(in out):: str
    integer, intent(in) :: n
    integer :: lenstr
    integer :: nabs

    lenstr = len(str)
    nabs = iabs(n)

    if (nabs >= lenstr) then
        str = repeat(' ', lenstr)
        return
    end if

    if (n < 0) str = str(nabs+1:)//repeat(' ',nabs)  ! shift left

    if (n > 0) str = repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
    
    end subroutine

!******************************************************************************

subroutine str_insert(str, substr, loc)
    !! Inserts the string `substr` into the string `str` at position `loc`. 
    !! Characters in `str` starting at position `loc` are shifted right to
    !! make room for the inserted string. Trailing spaces of `substr` are 
    !! removed prior to insertion.
    
    character(len=*), intent(in out) :: str
    character(len=*), intent(in) :: substr
    integer, intent(in) :: loc
    character(len=len(str))  ::tempstr
    integer :: len_substr
    
    len_substr = len_trim(substr)
    tempstr = str(loc:)

    call str_shift(tempstr, len_substr)

    tempstr(1:len_substr) = substr(1:len_substr)
    str(loc:) = tempstr
    
    end subroutine

!******************************************************************************

subroutine str_del(str, substr, n)
    !! Deletes first `n` occurrences of substring `substr` from string `str` and
    !! shifts characters left to fill hole. If `n < 0`, all occurances are
    !! deleted.  If `n` is not explicitly provided, it defaults to removing the
    !! first occurrence. Trailing spaces or blanks are not considered part of
    !! `substr`.
    
    character(len=*), intent(in out) :: str
    character(len=*),     intent(in) :: substr
    integer, optional,    intent(in) :: n
    integer :: n_
    integer :: lensubstr
    integer :: ipos
    integer :: cntr
    
    n_ = 1
    if (present(n)) n_ = n

    lensubstr = len_trim(substr)
    cntr = 0
    
    do
        if ((n_ > 0) .and. (cntr > n_)) exit

        ipos = index(str,substr)

        if (ipos == 0) exit

        str = str(:ipos-1)//str(ipos+lensubstr:)
        cntr = cntr + 1
    end do   

    end subroutine

!**********************************************************************

subroutine str_strip_comment(str, comment_str)
    !!  Strips trailing comment from a string.
    !!
    !!  The comment is assumed to begin with the sequence of characters in
    !!  `comment_str`. If the sequence `comment_str` is not found within `str`,
    !!  no changes are made.

    character(len=*), intent (in out) :: str
    !!  Input string
    character(len=*), intent (in) :: comment_str
    !!  String indicating beginning of a comment.
    integer :: ipos

    ipos = index(adjustl(str),comment_str)

    if (ipos /= 0) then
        str = str(1:(ipos-1))
    end if

    end subroutine

!**********************************************************************

subroutine str_get_keyval(str, key, val, delimiter)
    !! Split a string `str` into two strings, `key` and `val` based on space
    !! delimiter.
    !!
    !! A non-empty non-comment string should be passed to this subroutine.
    !! Keys can have corresponding empty values, but keys must always be present

    character (len=*), intent (in) :: str
    character (len=:), allocatable, intent (out) :: key
    character (len=:), allocatable, intent (out) :: val
    character (len=*), intent (in), optional :: delimiter
    character (len=:), allocatable :: delimiter_
    character (len=:), allocatable :: str_just
    integer :: m
    integer :: n
    
    !blankspace is represented as the integer 32 in ascii chart.
    delimiter_ = achar(32)
    if (present(delimiter)) delimiter_ = delimiter

    str_just = trim(adjustl(str))
    n = len(str_just)

    m = index(str_just, delimiter_)

    if (m == 0) then
        key = str_just
        val = ''
    else
        key = trim(str_just(1:m-1))
        val  = str_just(m+len_trim(delimiter_):n)
    end if

    val = trim(adjustl(val))

    end subroutine

!******************************************************************************

subroutine str_match(str, ipos, imatch)
    !! This routine finds the delimiter in string `str` that matches the delimiter
    !! in position `ipos` of `str`. The argument `imatch` contains the position of
    !! the matching delimiter. Allowable delimiters are (), [], {}, <>.
    
    character(len=*), intent(in) :: str
    integer,          intent(in) :: ipos
    integer,         intent(out) :: imatch
    character(len=1) :: delim1
    character(len=1) :: delim2
    character(len=1) :: ch
    integer :: lenstr
    integer :: istart
    integer :: iend
    integer :: inc
    integer :: idelim2
    integer :: isum
    integer :: i
    
    lenstr = len_trim(str)

    delim1 = str(ipos:ipos)

    select case(delim1)
        case('(')
            idelim2=iachar(delim1)+1
            istart=ipos+1
            iend=lenstr
            inc=1
        case(')')
            idelim2=iachar(delim1)-1
            istart=ipos-1
            iend=1
            inc=-1
        case('[','{','<')
            idelim2=iachar(delim1)+2
            istart=ipos+1
            iend=lenstr
            inc=1
        case(']','}','>')
            idelim2=iachar(delim1)-2
            istart=ipos-1
            iend=1
            inc=-1
        case default
            write(*,*) delim1,' is not a valid delimiter'
            return
    end select

    if (istart < 1 .or. istart > lenstr) then
        write(*,*) delim1,' has no matching delimiter'
        return
    end if
    delim2=achar(idelim2) ! matching delimiter
    
    isum = 1
    do i = istart, iend, inc
        ch=str(i:i)
        if (ch /= delim1 .and. ch /= delim2) cycle
        if (ch == delim1) isum = isum+1
        if (ch == delim2) isum = isum-1
        if (isum == 0) exit
    end do

    if(isum /= 0) then
        write(*,*) delim1,' has no matching delimiter'
        return
    end if   

    imatch = i
    
    end subroutine

!**********************************************************************

pure function str_from_inum(num, frmt) result(str)

    integer, intent(in) :: num
    character(len=:), allocatable :: str
    character(len=*), optional, intent(in) :: frmt
    character(len=:), allocatable :: frmt_
    character(len=24) :: buf

    frmt_ = '(i0)'
    if (present(frmt)) frmt_ = frmt

    write(buf, frmt_) num

    str = trim(adjustl(buf))

    end function

!**********************************************************************

pure function str_from_ilnum(num, frmt) result(str)

    integer(ip_long), intent(in) :: num
    character(len=:), allocatable :: str
    character(len=*), optional, intent(in) :: frmt
    character(len=:), allocatable :: frmt_
    character(len=24) :: buf

    frmt_ = '(i0)'
    if (present(frmt)) frmt_ = frmt

    write(buf, frmt_) num

    str = trim(adjustl(buf))

    end function

!**********************************************************************

pure function str_from_dnum(num, frmt) result(str)

    real(rp), intent(in) :: num
    character(len=:), allocatable :: str
    character(len=*), optional, intent(in) :: frmt
    character(len=:), allocatable :: frmt_
    character(len=32) :: buf

    frmt_ = '(g0.15)'
    if (present(frmt)) frmt_ = frmt

    write(buf, frmt_) num

    str = str_trimzero(buf)

    end function

!******************************************************************************

subroutine str_compact_rlstr(str)
    !! author: Izaak Beekman
    !! date: 02/24/2015
    !!
    !! Compact a string representing a real number, so that the same value is
    !! displayed with fewer characters.

    character(len=*),intent(in out) :: str
        !! string representation of a real number.
    character(len=len(str)) :: significand
    character(len=len(str)) :: expnt
    character(len=2) :: separator
    integer :: exp_start
    integer :: decimal_pos
    integer :: sig_trim
    integer :: exp_trim
    integer :: i  !! counter

    str = adjustl(str)
    exp_start = scan(str, 'eEdD')
    if (exp_start == 0) exp_start = scan(str, '-+', back=.true.)
    decimal_pos = scan(str, '.')
    if (exp_start /= 0) separator = str(exp_start:exp_start)

    if ( exp_start < decimal_pos ) then !possibly signed, exponent-less float
        significand = str
        sig_trim = len_trim(significand)
        do i = len_trim(significand),decimal_pos+2,-1 !look from right to left at 0s
                                                       !but save one after the decimal place
            if (significand(i:i) == '0') then
                sig_trim = i-1
            else
                exit
            end if
        end do
        str = trim(significand(1:sig_trim))

    else if (exp_start > decimal_pos) then !float has exponent
        significand = str(1:exp_start-1)
        sig_trim = len_trim(significand)

        do i = len_trim(significand),decimal_pos+2,-1 !look from right to left at 0s
            if (significand(i:i) == '0') then
                sig_trim = i-1
            else
                exit
            end if
        end do

        expnt = adjustl(str(exp_start+1:))
        if (expnt(1:1) == '+' .or. expnt(1:1) == '-') then
            separator = trim(adjustl(separator))//expnt(1:1)
            exp_start = exp_start + 1
            expnt     = adjustl(str(exp_start+1:))
        end if

        exp_trim = 1

        do i = 1,(len_trim(expnt)-1) !look at exponent leading zeros saving last
            if (expnt(i:i) == '0') then
                exp_trim = i+1
            else
                exit
            end if
        end do
        str = trim(adjustl(significand(1:sig_trim)))// &
              trim(adjustl(separator))// &
              trim(adjustl(expnt(exp_trim:)))

    !else
        ! mal-formed real, BUT this code should be unreachable
    end if

    end subroutine

!*******************************************************************************

pure function str_trimzero(str) result(res)
    !! Deletes nonsignificant trailing zeroes from number string str. If number
    !! string ends in a decimal point, one trailing zero is added.

    character(len=*), intent(in) :: str
    character(len=:), allocatable :: res
    character(len=len(str)) :: buf
    character(len=10) :: sexp
    character(len=1) :: ch
    integer :: ipos
    integer :: lbuf
    integer :: i

    buf = str
    ipos = scan(str,'eE')

    if (ipos > 0) then
       sexp = buf(ipos:)
       buf = buf(1:ipos-1)
    endif

    lbuf = len_trim(buf)

    do i = lbuf, 1, -1
        ch = buf(i:i)
        if (ch == '0') cycle          
        if (ch == '.') then
            buf = buf(1:i)//'0'
            if (ipos > 0) buf = trim(buf)//trim(sexp)
            exit
        endif
        buf = buf(1:i)
        exit
    end do

    if(ipos > 0) buf = trim(buf)//trim(sexp)

    res = trim(adjustl(buf))

    end function

!***********************************************************************

pure function str_to_d(str) result(res)

    character(len=*), intent(in) :: str
    real(rp) :: res

    read(str,*) res

    end function

!***********************************************************************

pure function str_to_i(str) result(res)

    character(len=*), intent(in) :: str
    integer :: res

    read(str,*) res

    end function

!***********************************************************************

pure function str_strip(str, chars, ends) result (ostr)
    !! Returns a copy of string `str` with the leading and trailing characters
    !! removed. The `chars` argument is a string specifying the set of characters to
    !! be removed.  The `chars` argument is not a prefix or suffix; rather, all
    !! combinations of its values are stripped. If `ends = 'l'`, only leading
    !! characters are removed, if `ends = 'r'`, only trailing characters are
    !! removed, and if `ends = 'b'` both leading and trailing characters are
    !! removed.

    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: chars
    character(len=1), intent(in) :: ends
        !! {'l', 'r', 'b'} 
    character(len=:), allocatable :: ostr
    integer :: lenstr
    integer :: ibeg
    integer :: iend

    lenstr = len(str)

    select case (ends)
    case ('l')
        ibeg = verify(str, chars)
        iend = lenstr
    case ('r')
        ibeg = 1
        iend = verify(str, chars, .true.)
    case ('b')
        ibeg = verify(str, chars)
        iend = verify(str, chars, .true.)
    case default
        ibeg = 1
        iend = lenstr
    end select

    if ((ibeg==0) .or. (iend==0)) then
        ostr = ''
    else
        ostr = str(ibeg:iend)
    end if

    end function

!**********************************************************************

subroutine str_split(str, delimiter, before)
    !! Routine finds the first instance of a character from `delims` in the the
    !! string `str`. The characters before the found delimiter are output in
    !! `before`. The characters after the found delimiter are output in `str`.
    !! Repeated applications of this routine can be used to parse a string into its
    !! component parts. Multiple whitespaces of `str` are compacted into a single
    !! whitespace before splitting begins. If either `str` or `delimiter` is
    !! empty, an empty string is retured in `before` and `str` remains
    !! unchanged.
    
    character(len=*), intent(in out) :: str
    character(len=*), intent(in) :: delimiter
    character(len=:), allocatable, intent(out) :: before
    integer :: lenstr
    integer :: lendelim
    integer :: ipos
    
    str = str_compact(str)
    lenstr = len(str)
    lendelim = len(delimiter)

    if ( (lenstr == 0) .or. (lendelim == 0) ) then
        ! `str` or `delimiter` is empty
        before = ''
        return        
    end if

    ipos = index(str, delimiter)

    if (ipos == 0) then 
        ! string does not contain any delimiter
        before = ''
        return
    else
        before = str(1:(ipos-1))
        str = str((ipos+lendelim-1):)
    end if
    
    end subroutine

!**********************************************************************

subroutine str_append(dest, source, sep)
    !! Appends a copy of the `source` string to the `dest` string, with 
    !! optional string `sep` in between. It is assumed that `dest` is long
    !! enough to hold the result, otherwise an error will be generated.
    character(len=*),       intent(in out) :: dest
    character(len=*),           intent(in) :: source
    character(len=*), optional, intent(in) :: sep
    character(len=:), allocatable :: sep_
    integer :: len_dest
    integer :: len_source
    integer :: len_sep_
    integer :: ipos

    sep_ = ''
    if (present(sep)) sep_ = sep

    len_dest = len_trim(dest)
    len_source = len_trim(adjustl(source))
    len_sep_ = len(sep_)

    ipos  = len_dest + 1

    if (len_sep_ > 0) then
        dest( ipos:(ipos+len_sep_) ) = sep_
        ipos  = ipos + len_sep_ + 1
    end if

    dest( ipos:(ipos+len_source) ) = trim(adjustl(source))

    end subroutine

!**********************************************************************

pure function str_startswith(str, prefix, start, finish) result(res)
    !! Returns `.true.` if the string `str` starts with `prefix`, otherwise
    !! returns `.false.`. With optional `start`, test beginning at that position.
    !! With optional `finish`, stop comparing beyond that position.

    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: prefix
    integer, optional, intent(in) :: start
    integer, optional, intent(in) :: finish
    logical :: res
    integer :: ibeg
    integer :: iend

    ibeg = 1; iend = len(str)

    if (present(start)) ibeg = start
    if (present(finish)) iend = finish

    !Return .false. if prefix is longer than str(ibeg:iend)
    if (len(prefix) > (iend-ibeg+1)) then
        res = .false.
        return
    end if

    if (index(str(ibeg:iend), prefix) == 1) then
        res = .true.
    else
        res = .false.
    end if
    
    end function

!**********************************************************************

pure function str_endswith(str, suffix, start, finish) result(res)
    !! Returns `.true.` if the string `str` ends with `suffix`, otherwise
    !! return `.false.`. With optional `start`, test beginning at that position.
    !! With optional `finish`, stop comparing beyond that position.

    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: suffix
    integer, optional, intent(in) :: start
    integer, optional, intent(in) :: finish
    logical :: res
    integer :: ibeg
    integer :: iend
    integer :: iloc

    ibeg = 1; iend = len(str)

    if (present(start)) ibeg = start
    if (present(finish)) iend = finish

    !Return .false. if suffix is longer than str(ibeg:iend)
    if (len(suffix) > (iend-ibeg+1)) then
        res = .false.
        return
    end if

    ! Getting last occurrance of suffix
    iloc = index(str(ibeg:iend), suffix, back=.true.) 

    if ( (iloc+len(suffix)-1) == iend) then
        res = .true.
    else
        res = .false.
    end if

    end function

!******************************************************************************

subroutine readline(nunitr, line, comment_str, ios)
    !!  Reads a line from unit=nunitr, ignoring blank lines
    !!  and deleting comments

    integer, intent(in) :: nunitr
    character(len=*), intent(in out):: line
    character(len=*), intent(in) :: comment_str
    integer, intent(out) :: ios
    
    do  
        read(nunitr,'(a)', iostat=ios) line      ! read input line

        if(ios /= 0) return

        if ((len_trim(line) /= 0) .and. (.not. str_is_comment(line, comment_str))) then
            call str_strip_comment(line, comment_str)
            exit
        end if

    end do
    
    end subroutine

!******************************************************************************

end module m_strings
