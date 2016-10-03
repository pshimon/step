!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! adopted from numerical recipies

SUBROUTINE LUDCMP1(A,N,INDX,VV,FLAG)
    USE FDEFS
    INTEGER:: N,INDX(N),FLAG
    REAL(F64):: A(N,N),VV(N)
    REAL(F64),PARAMETER::SMALL=1.0D-20
    INTEGER:: I,IMAX=0,J,K
    REAL(F64):: AAMAX,DUM,SUM1
    FLAG=0
    DO I=1,N
        AAMAX=0.0D0
        DO  J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        ENDDO
        IF (AAMAX.EQ.0.0D0) THEN
            FLAG=1
            RETURN
        ENDIF
        VV(I)=1./AAMAX
    ENDDO
    DO  J=1,N
        DO  I=1,J-1
            SUM1=A(I,J)
            DO  K=1,I-1
                SUM1=SUM1-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=SUM1
        ENDDO
        AAMAX=0.0D0
        DO  I=J,N
            SUM1=A(I,J)
            DO  K=1,J-1
                SUM1=SUM1-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=SUM1
            DUM=VV(I)*ABS(SUM1)
            IF (DUM.GE.AAMAX) THEN
                IMAX=I
                AAMAX=DUM
            ENDIF
        ENDDO
        IF (J.NE.IMAX)THEN
            DO  K=1,N
                DUM=A(IMAX,K)
                A(IMAX,K)=A(J,K)
                A(J,K)=DUM
            ENDDO
            VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.0D0) A(J,J)=SMALL
        IF(J.NE.N)THEN
            DUM=1./A(J,J)
            DO  I=J+1,N
                A(I,J)=A(I,J)*DUM
            ENDDO
        ENDIF
    ENDDO
ENDSUBROUTINE LUDCMP1

SUBROUTINE LUBKSB1(A,N,INDX,B)
    USE FDEFS
    INTEGER:: N,INDX(N)
    REAL(F64):: A(N,N),B(N)
    INTEGER I,II,J,LL
    REAL(F64):: SUM1
    II=0
    DO  I=1,N
        LL=INDX(I)
        SUM1=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
            DO  J=II,I-1
                SUM1=SUM1-A(I,J)*B(J)
            ENDDO
        ELSE IF (SUM1.NE.0.0D0) THEN
            II=I
        ENDIF
        B(I)=SUM1
    ENDDO
    DO  I=N,1,-1
        SUM1=B(I)
        DO  J=I+1,N
          SUM1=SUM1-A(I,J)*B(J)
        ENDDO
        B(I)=SUM1/A(I,I)
    ENDDO
ENDSUBROUTINE LUBKSB1
      

