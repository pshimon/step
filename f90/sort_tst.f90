PROGRAM sort_tst
use sort_lib
use rng
INTEGER :: left, right,flag,n,n1,i
integer,parameter::m=4
REAL(F32) :: t1, t2,times(2,m)
REAL(F64), ALLOCATABLE :: a(:)
INTEGER ::seed
INTERFACE
    subroutine sortquick(a,l,r)  bind(c,name='sortSelDbl')
    IMPORT::C_DOUBLE,C_INT
    REAL(C_DOUBLE)::a(*)
    integer(c_int),intent(in),value::l,r
    ENDsubroutine sortquick

ENDINTERFACE



    n=100
    n1=10
    do i=1,m
    left = 1
    ALLOCATE(a(n), STAT=flag)
    IF (flag /= 0) THEN
        print *,'allocation problem'
        stop
    END IF

    left = 1
    right = n
    seed=-1
    call INIT_RANDOM_SEED(seed)
    CALL random_number(a)
    CALL cpu_time(t1)
  !  CALL SORT_QUICK_NR_F64(A,LEFT, RIGHT)
  !  CALL SORT_HEAP_F64(A,LEFT, RIGHT)
 !   CALL SORT_SHELL_F64(A,LEFT, RIGHT)
 !   CALL SORT_INS_F64(A,LEFT, RIGHT)    
 !   CALL SORT_SEL_F64(A,LEFT, RIGHT)
 !   CALL SORT_QUICK_F64(A,LEFT, RIGHT)
   CALL sortquick(a,left-1, right-1)
    CALL cpu_time(t2)
    CALL check
    times(1,i) = t2 - t1
 ! now sort the sorted list
    CALL cpu_time(t1)
!    CALL SORT_QUICK_NR_F64(A,LEFT, RIGHT)
 !   CALL SORT_HEAP_F64(A,LEFT, RIGHT)
!    CALL SORT_SHELL_F64(A,LEFT, RIGHT)    
  !  CALL SORT_INS_F64(A,LEFT, RIGHT)
 !   CALL SORT_SEL_F64(A,LEFT, RIGHT)
 !   CALL SORT_QUICK_F64(A,LEFT, RIGHT)
    CALL sortquick(a,left-1, right-1)
    CALL cpu_time(t2)
    CALL check
    times(2,i) = t2 - t1
    print *,n,times(1,i),times(2,i)
    DEALLOCATE(a)
    n=n*n1
    enddo
! Get ready to go around again
 
CONTAINS
  SUBROUTINE check
  INTEGER :: i
  REAL(F64) :: val
  val = a(left)
  DO i = left+1, right
    IF( val <= a(i)) THEN
      val = a(i)
    ELSE
      WRITE(*, '(''list not sorted: '',i8 ,'' > '',e12.4 &
	       &, '' for i = '',i8)')val, a(i), i
      STOP
    END IF
  END DO

  END SUBROUTINE check

 
  

  ENDPROGRAM sort_tst
