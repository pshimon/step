!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations        !
! http://industrialphys.com                                !
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADOPTED FROM NUMERICAL RECIPIES

SUBROUTINE GAUSS_JORDAN(A,N,B,M)
    USE FDEFS
    IMPLICIT NONE
    REAL(F64) :: A(N,N),B(N,M)
    INTEGER::N,M
    INTEGER, DIMENSION(N) :: IPIV,INDXR,INDXC
    INTEGER:: I,ICOL,IROW,J,K,L,LL                            
    REAL(F64):: BIG,DUM,PIVINV 
    IPIV=0
    ICOL=0
    IROW=0
    DO  I=1,N                                 
        BIG=0.0                                     
        DO  J=1,N                            
            IF(IPIV(J).NE.1)THEN                 
                DO  K=1,N
                    IF (IPIV(K).EQ.0) THEN
                        IF (ABS(A(J,K)).GE.BIG)THEN
                            BIG=ABS(A(J,K))
                            IROW=J
                            ICOL=K
                        ENDIF
                    ENDIF
                ENDDO 
            ENDIF
        ENDDO 
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
            DO  L=1,N
                DUM=A(IROW,L)
                A(IROW,L)=A(ICOL,L)
                A(ICOL,L)=DUM
            ENDDO 
            DO  L=1,M
               DUM=B(IROW,L)
               B(IROW,L)=B(ICOL,L)
               B(ICOL,L)=DUM
            ENDDO 
        ENDIF
        INDXR(I)=IROW                         
        INDXC(I)=ICOL                              
        PIVINV=1.0D0/A(ICOL,ICOL)
        A(ICOL,ICOL)=1.0D0
        A(ICOL,:)=A(ICOL,:)*PIVINV
        B(ICOL,:)=B(ICOL,:)*PIVINV
        DO  LL=1,N                          
            IF(LL.NE.ICOL)THEN               
                DUM=A(LL,ICOL)
                A(LL,ICOL)=0.0D0
                A(LL,:)=A(LL,:)-A(ICOL,:)*DUM
                B(LL,:)=B(LL,:)-B(ICOL,:)*DUM
            ENDIF
        ENDDO 
    ENDDO                       
    DO  L=N,1,-1                             
        IF(INDXR(L).NE.INDXC(L))THEN               
            DO  K=1,N                           
                DUM=A(K,INDXR(L))                
                A(K,INDXR(L))=A(K,INDXC(L))
                A(K,INDXC(L))=DUM
          ENDDO 
        ENDIF
    ENDDO 
    RETURN                                     
END SUBROUTINE GAUSS_JORDAN

