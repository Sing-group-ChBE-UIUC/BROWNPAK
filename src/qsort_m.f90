!*************************************************************************!
!                                                                         !
! This is free and unencumbered software released into the public domain. !
!                                                                         !
! Anyone is free to copy, modify, publish, use, compile, sell, or         !
! distribute this software, either in source code form or as a compiled   !
! binary, for any purpose, commercial or non-commercial, and by any       !
! means.                                                                  !
!                                                                         !
! In jurisdictions that recognize copyright laws, the author or authors   !
! of this software dedicate any and all copyright interest in the         !
! software to the public domain. We make this dedication for the benefit  !
! of the public at large and to the detriment of our heirs and            !
! successors. We intend this dedication to be an overt act of             !
! relinquishment in perpetuity of all present and future rights to this   !
! software under copyright law.                                           !
!                                                                         !
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         !
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      !
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  !
! IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR       !
! OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   !
! ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR   !
! OTHER DEALINGS IN THE SOFTWARE.                                         !
!                                                                         !
! For more information, please refer to <http://unlicense.org>            !
!                                                                         !
!*************************************************************************!

MODULE qsort_m
    !!  Implements quicksort for a sequence of integers and reals, in combination with
    !!  insertion sort for very short sequences.
    !!
    !! - Quick sort routine from Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990)
    !! "Programmer's Guide to Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.               
    !! - Modified by Alan Miller to include an associated integer array which gives
    !!   the positions of the elements in the original order.  
    !! - Modified for integer array by Sarit Dutta                                  

USE constants_m

IMPLICIT NONE

CONTAINS

!******************************************************************************

RECURSIVE SUBROUTINE iqsort(list, order)
    !!  Sorts a sequence of integers

IMPLICIT NONE

INTEGER, DIMENSION (:), INTENT(IN OUT)  :: list
    !!  Sequence of integers to be sorted
INTEGER, DIMENSION (:), INTENT(OUT), OPTIONAL  :: order
    !!  Indices of the sorted sequence

! Local variable
INTEGER :: i

IF (PRESENT(order)) THEN
    DO i = 1, SIZE(list)
      order(i) = i
    END DO
END IF

CALL quick_sort_1(1, SIZE(list))

CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER             :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 8
    
    IF (right_end < left_end + max_simple_sort_size) THEN
        ! Use interchange sort for small lists
        CALL interchange_sort(left_end, right_end)
    
    ELSE
        ! Use partition ("quick") sort
        reference = list((left_end + right_end)/2)
        i = left_end - 1; j = right_end + 1
    
        DO
            ! Scan list from left end until element >= reference is found
            DO
                i = i + 1
                IF (list(i) >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
                j = j - 1
                IF (list(j) <= reference) EXIT
            END DO
    
            IF (i < j) THEN
                ! Swap two out-of-order elements
                temp = list(i); list(i) = list(j); list(j) = temp
                IF (PRESENT(order)) THEN
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            ELSE IF (i == j) THEN
                i = i + 1
                EXIT
            ELSE
                EXIT
            END IF
        END DO
    
        IF (left_end < j) CALL quick_sort_1(left_end, j)
        IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
    END SUBROUTINE quick_sort_1


    SUBROUTINE interchange_sort(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER                :: temp
    
    DO i = left_end, right_end - 1
        DO j = i+1, right_end
            IF (list(i) > list(j)) THEN
                temp = list(i); list(i) = list(j); list(j) = temp
                IF (PRESENT(order)) THEN
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            END IF
        END DO
    END DO
    
    END SUBROUTINE interchange_sort

END SUBROUTINE iqsort

!******************************************************************************

RECURSIVE SUBROUTINE dqsort(list, order)
    !!  Sorts a sequence of reals
    !""

IMPLICIT NONE

REAL(RP), DIMENSION (:), INTENT(IN OUT)  :: list
    !!  Sequence of reals to be sorted
INTEGER, DIMENSION (:), INTENT(OUT), OPTIONAL  :: order
    !!  Indices of the sorted sequence

! Local variable
INTEGER :: i

IF (PRESENT(order)) THEN
    DO i = 1, SIZE(list)
      order(i) = i
    END DO
END IF

CALL quick_sort_1(1, SIZE(list))

CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(RP)            :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 8
    
    IF (right_end < left_end + max_simple_sort_size) THEN
        ! Use interchange sort for small lists
        CALL interchange_sort(left_end, right_end)
    
    ELSE
        ! Use partition ("quick") sort
        reference = list((left_end + right_end)/2)
        i = left_end - 1; j = right_end + 1
    
        DO
            ! Scan list from left end until element >= reference is found
            DO
                i = i + 1
                IF (list(i) >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
                j = j - 1
                IF (list(j) <= reference) EXIT
            END DO
    
            IF (i < j) THEN
                ! Swap two out-of-order elements
                temp = list(i); list(i) = list(j); list(j) = temp
                IF (PRESENT(order)) THEN
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            ELSE IF (i == j) THEN
                i = i + 1
                EXIT
            ELSE
                EXIT
            END IF
        END DO
    
        IF (left_end < j) CALL quick_sort_1(left_end, j)
        IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
    END SUBROUTINE quick_sort_1


    SUBROUTINE interchange_sort(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(RP)            :: temp
    
    DO i = left_end, right_end - 1
        DO j = i+1, right_end
            IF (list(i) > list(j)) THEN
                temp = list(i); list(i) = list(j); list(j) = temp
                IF (PRESENT(order)) THEN
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            END IF
        END DO
    END DO
    
    END SUBROUTINE interchange_sort

END SUBROUTINE dqsort

!******************************************************************************

END MODULE qsort_m
