!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define ORDER(A, B) IF (BEFORE(B, A)) THEN SWAP(A, B) ENDIF
#define BEFORE(A, B) ((A) < (B))
#define SWAP(A, B)   SWAPTMP = A; A = B; B = SWAPTMP 

MODULE SORT_LIB
    USE DATA_BUF
    IMPLICIT NONE
    INTEGER, PARAMETER:: SWITCHSORT=10
    
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
    IF(MINI>L) SWAP(A(L), A(MINI));
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

    
ENDMODULE SORT_LIB
