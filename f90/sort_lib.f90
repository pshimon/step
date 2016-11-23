MODULE sort_lib
    integer,parameter::wp=kind(1.0d0)

CONTAINS

  SUBROUTINE sort_quick(a, left, right)
  INTEGER, INTENT(IN) :: left, right
  REAL(wp), INTENT(INOUT) :: a(left:right)
  REAL(wp) :: savedVal

  CALL quicksort(left, right)
  CALL insertion

  CONTAINS
    
    SUBROUTINE insertion
    INTEGER :: i, j

    DO i = left+1, right
      CALL compex(left, i)
    END DO

    DO i = left+2, right
      j = i
      CALL saveValue(i)
      DO WHILE(compareValue(j-1))
        CALL moveValue(j, j-1)
        j = j-1
      END DO

      CALL restoreValue(j)
    END DO

    END SUBROUTINE insertion

    RECURSIVE SUBROUTINE quicksort(left, right)
    INTEGER, INTENT(IN) :: left, right

    INTEGER, PARAMETER:: switchsorts=10
    INTEGER :: i

    IF((right-left) > switchsorts) THEN
      CALL exchange((right+left)/2, (right-1))
      CALL compex(left, right-1)
      CALL compex(right, left)
      CALL compex(right-1, right)
      i = partition(left+1, right-1)
      CALL quicksort(left, i-1)
      CALL quicksort(i+1, right)
    END IF

    END SUBROUTINE quicksort

    FUNCTION partition(left, right) RESULT(i)
    INTEGER, INTENT(IN) :: left, right
    INTEGER :: i, j

    i = left - 1
    j = right
    CALL saveValue(right)

    DO
      DO
        i = i+1
        IF(i>right) EXIT
        IF(compareValue(i)) EXIT
      END DO

      DO
        j = j-1
        IF(.NOT.compareValue(j) .OR. j==left) EXIT
      END DO

      IF(i>= j) EXIT
      CALL exchange(i,j)

    END DO

    CALL exchange(i, right)

    END FUNCTION partition

    
    SUBROUTINE exchange(i,j)
! Exchange the contents of the ith and jth elements
    INTEGER, INTENT(IN) :: i,j
    REAL(wp) :: t
    t = a(i)
    a(i) = a(j)
    a(j) = t
    END SUBROUTINE exchange
  
    LOGICAL FUNCTION compare(i,j)
! Compare the contents of the ith and jth elements
! This determines the final sorting order.
! Code for ascending order.
    INTEGER, INTENT(IN) :: i,j
    compare = a(i) < a(j)
    END FUNCTION compare
  
    SUBROUTINE compex(i,j)
    INTEGER, INTENT(IN) :: i,j
    IF(compare(j,i)) CALL exchange(i,j)
    END SUBROUTINE compex
  
    SUBROUTINE moveValue(i,j)
! Overwrite the contents of the jth element
! with the contents of the ith element
    INTEGER, INTENT(IN) :: i,j
    a(i) = a(j)
    END SUBROUTINE moveValue
  
! The next three subprograms are used to store, 
! compare against and restore a particular element.
    LOGICAL FUNCTION compareValue(j)
    INTEGER, INTENT(IN) :: j
    compareValue = savedVal < a(j)
    END FUNCTION compareValue
  
    SUBROUTINE saveValue(i)
    INTEGER, INTENT(IN) :: i
    savedVal = a(i)
    END SUBROUTINE saveValue
  
    SUBROUTINE restoreValue(i)
    INTEGER, INTENT(IN) :: i
    a(i) = savedVal
    END SUBROUTINE restoreValue
  
  END SUBROUTINE sort_quick

END MODULE sort_lib
