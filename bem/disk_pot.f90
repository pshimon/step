!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE   DISK_POT(RES,R,Z,A) 
    USE FDEFS
    IMPLICIT NONE    
    REAL(F64)::R,Z,A,RES,RES2
    REAL(F64),PARAMETER::PHI0=ZERO_F64,PHI1=HALF_F64*PI_F64
    INTEGER,PARAMETER:: JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1
    REAL,PARAMETER::EPS=1.E-6
    CALL QROMB(RES2)
    RES=TWO_F64*RES2
    CONTAINS
    FUNCTION DP(PHI) RESULT(RES)
    REAL(F64)::PHI,RES
    REAL(F64)::AR1,AR2,W2,RC,R0,R1,R2
    RC=R*COS(PHI)
    AR1=A-RC
    AR2=A+RC
    W2=R**2+Z**2
    R0=SQRT(W2)
    R1=SQRT(A**2-TWO_F64*A*RC+W2)
    R2=SQRT(A**2+TWO_F64*A*RC+W2)
    RES=R1+R2-TWO_F64*R0+RC*LOG((R1+AR1)/(R0-RC)*(R0+RC)/(R2+AR2))
    ENDFUNCTION DP
    SUBROUTINE TRAPZD(S,N)
      INTEGER:: N
      REAL(F64):: S
      INTEGER IT,J
      REAL(F64):: DEL,ACC,TNM,X
      IF (N.EQ.1) THEN
        S=0.5*(PHI1-PHI0)*(DP(PHI0)+DP(PHI1))
      ELSE
        IT=2**(N-2)
        TNM=IT
        DEL=(PHI1-PHI0)/TNM
        X=PHI0+0.5*DEL
        ACC=0.
        DO  J=1,IT
          ACC=ACC+DP(X)
          X=X+DEL
        ENDDO
        S=0.5*(S+(PHI1-PHI0)*ACC/TNM)
      ENDIF
      RETURN
    ENDSUBROUTINE TRAPZD
    SUBROUTINE QROMB(SS)
    REAL(F64):: SS
    INTEGER:: J
    REAL(F64):: DSS,H(JMAXP),S(JMAXP)
      H(1)=1.
      DO J=1,JMAX
        CALL TRAPZD(S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.,SS,DSS)
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
    ENDDO
    ENDSUBROUTINE QROMB
ENDSUBROUTINE DISK_POT

FUNCTION  DISK_POT0(Z,A) RESULT(RES)
    USE FDEFS
    REAL(F64)::Z,A,RES
    REAL(F64),PARAMETER::FCT=TWO_F64*PI_F64;
    RES=FCT*(SQRT(A**2+Z**2)-ABS(Z))
ENDFUNCTION DISK_POT0

