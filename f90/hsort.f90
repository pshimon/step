!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! adopted from numerical recipes
! SORTS ARRAY A IN ASSCENDING ORDER
! IND CONTAINS PERMUTATION
SUBROUTINE HSORT_F32(N,A,IND)
    USE FDEFS
    INTEGER,INTENT(IN)::N
    REAL(F32):: A(N)
    INTEGER::IND(N)
    INTEGER ::I
    REAL(F32)::AA
    INTEGER::II
    DO I=N/2,1,-1
        CALL SIFT_DOWN(I,N)
    END DO
    DO I=N,2,-1
        AA=A(I)
        A(I)=A(1)
        A(1)=AA
        II=IND(I)
        IND(I)=IND(1)
        IND(1)=II 
        CALL SIFT_DOWN(1,I-1) 
    END DO
contains
    SUBROUTINE SIFT_DOWN(L,R)
    INTEGER,INTENT(IN)::L,R
    INTEGER:: J,JOLD,I1
    REAL(F32) :: A1
    A1=A(L)
    I1=IND(L)
    JOLD=L
    J=L+L
    DO 
        IF (J > R) EXIT
        IF (J < R) THEN
            IF (A(J) < A(J+1)) J=J+1 
        END IF
        IF (A1 >= A(J)) EXIT 
        A(JOLD)=A(J) 
        IND(JOLD)=IND(J)
        JOLD=J
        J=J+J
    END DO
    A(JOLD)=A1 
    IND(JOLD)=I1
    ENDSUBROUTINE SIFT_DOWN
END SUBROUTINE HSORT_F32



