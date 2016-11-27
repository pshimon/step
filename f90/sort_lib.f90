!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define BEFORE(A, B) ((A) < (B))
#define SWAP(A, B)   SWAPTMP = A; A = B; B = SWAPTMP 

MODULE SORT_LIB
    USE DATA_BUF
    IMPLICIT NONE
    INTEGER, PARAMETER:: SWITCHSORT=10
    INTEGER, PARAMETER:: NSTACK=32    
CONTAINS
    SUBROUTINE SORT_SEL_F64(A, L, R)
    INTEGER, INTENT(IN) :: L, R
    REAL(F64), INTENT(INOUT) :: A(:)
    INTEGER::I,J,MINI
    REAL(F64)::SWAPTMP
    DO I=L,R-1
        MINI=I
        DO J=I+1,R
            IF(BEFORE(A(J), A(MINI))) MINI = J
        ENDDO
         SWAP(A(I),A(MINI))
    ENDDO
    ENDSUBROUTINE SORT_SEL_F64

    SUBROUTINE SORT_INS_F64(A, L, R) 
    INTEGER, INTENT(IN) :: L, R
    REAL(F64), INTENT(INOUT) :: A(:)
    INTEGER::I,J,MINI
    REAL(F64)::SWAPTMP,V
    MINI=L
    DO I=L,R !FIND LOCATION OF ELEMENT TO BE THE FIRST
        IF(BEFORE(A(I), A(MINI))) MINI = I
    ENDDO
    IF(MINI>L) then 
        SWAP(A(L), A(MINI)) 
    endif
    ! ABOVE STEP IS NOT NECESSARY, BUT USEFUL: ONLY ONE CONDITION TO CHECK IN WHILE LOOP BELOW 
    DO I=L+2,R
        J = I
        V = A(I)
        DO WHILE (BEFORE(V, A(J-1))) 
            A(J) = A(J-1)
            J=J-1
        ENDDO
        A(J) = V; 
    ENDDO
    ENDSUBROUTINE SORT_INS_F64  

    SUBROUTINE SORT_SHELL_F64(A, L, R) 
    INTEGER, INTENT(IN) :: L, R
    REAL(F64), INTENT(INOUT) :: A(:)
    INTEGER::I,J,H
    REAL(F64)::V
    H=1
    DO WHILE(H<=(R-L)/9)
        H=3*H+1
    ENDDO ! finding initial increment
    DO WHILE(H>0) 
        DO I=L+H,R
            J=I
            V=A(I)
            DO WHILE((J>=L+H).AND.BEFORE(V, A(J-H)))
                A(J)=A(J-H)
                J=J-H
            ENDDO
            A(J)=V
        ENDDO
        H=H/3
    ENDDO
    ENDSUBROUTINE SORT_SHELL_F64

    SUBROUTINE SORT_HEAP_F64(A, L, R)
    INTEGER, INTENT(IN) :: L, R
    REAL(F64), INTENT(INOUT) :: A(:)
    INTEGER::I
    REAL(F64)::SWAPTMP
    DO I=(R+1)/2,L,-1
        CALL SIFT_DOWN(I,R)
    END DO
    DO I=R,L+1,-1
        SWAP(A(L),A(I))
        CALL SIFT_DOWN(L,I-1) 
    END DO
contains
    SUBROUTINE SIFT_DOWN(L,R)
    INTEGER,INTENT(IN)::L,R
    INTEGER:: J,JOLD
    REAL(F64) :: A1
    A1=A(L)
    JOLD=L
    J=L+L
    DO 
        IF (J > R) EXIT
        IF (J < R) THEN
            IF (A(J) < A(J+1)) J=J+1 
        END IF
        IF (A1 >= A(J)) EXIT 
        A(JOLD)=A(J) 
        JOLD=J
        J=J+J
    END DO
    A(JOLD)=A1 
    ENDSUBROUTINE SIFT_DOWN
END SUBROUTINE SORT_HEAP_F64
! adopted from hanson-hopkins
SUBROUTINE SORT_QUICK_F64(A, LEFT, RIGHT)
    INTEGER, INTENT(IN) :: LEFT, RIGHT
    REAL(F64), INTENT(INOUT) :: A(LEFT:RIGHT)
    REAL(F64) :: SAVEDVAL
    CALL QUICKSORT(LEFT, RIGHT)
    CALL SORT_INS_F64(A, LEFT, RIGHT)
  CONTAINS
    RECURSIVE SUBROUTINE QUICKSORT(LEFT, RIGHT)
        INTEGER, INTENT(IN) :: LEFT, RIGHT
        INTEGER :: I
        IF((RIGHT-LEFT) > SWITCHSORT) THEN
          CALL EXCHANGE((RIGHT+LEFT)/2, (RIGHT-1))
          CALL COMPEX(LEFT, RIGHT-1)
          CALL COMPEX(RIGHT, LEFT)
          CALL COMPEX(RIGHT-1, RIGHT)
          I = PARTITION(LEFT+1, RIGHT-1)
          CALL QUICKSORT(LEFT, I-1)
          CALL QUICKSORT(I+1, RIGHT)
        END IF
    END SUBROUTINE QUICKSORT

    FUNCTION PARTITION(LEFT, RIGHT) RESULT(I)
        INTEGER, INTENT(IN) :: LEFT, RIGHT
        INTEGER :: I, J
        I = LEFT - 1
        J = RIGHT
        SAVEDVAL = A(RIGHT)
        DO
            DO
                I = I+1
                IF(I>RIGHT) EXIT
                IF(SAVEDVAL < A(I)) EXIT
            ENDDO
            DO
                J = J-1
                IF(SAVEDVAL >= A(J) .OR. J==LEFT) EXIT
            END DO
            IF(I>= J) EXIT
            CALL EXCHANGE(I,J)
        END DO
        CALL EXCHANGE(I, RIGHT)
    END FUNCTION PARTITION

    SUBROUTINE EXCHANGE(I,J)
    INTEGER, INTENT(IN) :: I,J
    REAL(F64) :: T
    T = A(I)
    A(I) = A(J)
    A(J) = T
    END SUBROUTINE EXCHANGE
  
    LOGICAL FUNCTION COMPARE(I,J)
    INTEGER, INTENT(IN) :: I,J
    COMPARE = A(I) < A(J)
    END FUNCTION COMPARE
  
    SUBROUTINE COMPEX(I,J)
    INTEGER, INTENT(IN) :: I,J
    IF(COMPARE(J,I)) CALL EXCHANGE(I,J)
    END SUBROUTINE COMPEX
ENDSUBROUTINE SORT_QUICK_F64

#define STACK_INIT jstack=0
#define STACK_NOT_EMPTY (jstack>0)
#define PUSH(i) jstack=jstack+1;istack(jstack)=(i)
#define PUSH2(A, B)  PUSH(B); PUSH(A)
#define POP(A) A=istack(jstack);jstack=jstack-1
#define POP2(A,B) POP(A);POP(B)
SUBROUTINE SORT_QUICK_NR_F64(A, LEFT, RIGHT)
    INTEGER, INTENT(IN) :: LEFT, RIGHT
    REAL(F64), INTENT(INOUT) :: A(LEFT:RIGHT)
    CALL QUICKSORTNR(LEFT, RIGHT)
    CALL SORT_INS_F64(A, LEFT, RIGHT)
CONTAINS
    SUBROUTINE QUICKSORTNR(LEFT, RIGHT)
    INTEGER, INTENT(IN) :: LEFT, RIGHT
    INTEGER:: I,IR,J,JSTACK,K,L,ISTACK(NSTACK)
    REAL(F64):: V,SWAPTMP
    STACK_INIT
    L=LEFT
    IR=RIGHT
    PUSH2(L, IR)
    DO WHILE STACK_NOT_EMPTY
        POP2(L,IR)
        IF(IR-L<=SWITCHSORT) CYCLE
        K=(L+IR)/2
        SWAP(A(K),A(L+1))
        IF BEFORE(A(IR),A(L)) THEN
            SWAP(A(L), A(IR))
        ENDIF
        IF BEFORE(A(IR),A(L+1)) THEN
            SWAP(A(L+1), A(IR))
        ENDIF
        IF BEFORE(A(L+1),A(L)) THEN
            SWAP(A(L), A(L+1))
        ENDIF
         I=L+2
        J=IR
        V=A(L+1)
        DO
            DO WHILE BEFORE(A(I),V)
                I=I+1
            ENDDO
            DO WHILE BEFORE(V,A(J))
                J=J-1
            ENDDO
            IF(J<I) EXIT
            SWAP(A(I),A(J))
        ENDDO
        SWAP(A(L+1),A(J))
        IF(IR-I+1.GE.J-L)THEN
            PUSH2(IR,I)
            IR=J-1
        ELSE
            PUSH2(J-1,L)
            L=I
        ENDIF
    ENDDO
    ENDSUBROUTINE QUICKSORTNR

ENDSUBROUTINE SORT_QUICK_NR_F64
ENDMODULE SORT_LIB
