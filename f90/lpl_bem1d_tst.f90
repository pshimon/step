!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM LPL_BEM1D_TST
    USE TSURFD
    USE LINALG
    USE DATA_BUF
    USE VEC3D
    IMPLICIT NONE
    TYPE(TSURF_TYPE)::S     
    INTEGER::NT ! NUMBER OF TRIANGLES 
    INTEGER::NV ! NUMBER OF VERTICES
    REAL(F64),ALLOCATABLE,DIMENSION(:,:)::LM1
    REAL(F64),ALLOCATABLE,DIMENSION(:,:)::CP
    REAL(F64),ALLOCATABLE,DIMENSION(:)::QT1,E,EE,XP
    INTEGER,ALLOCATABLE::NCVEC(:),CTVEC(:,:),IPIV(:)
    INTEGER::CAC
    CHARACTER(LEN=100) :: ARG    
    INTEGER::FLAG,I
    REAL(F64)::TM0,TM1
    REAL(F64)::DLPL_GF_L1 
    EXTERNAL DLPL_GF_L1
    REAL(F64)::PHI,FCT,DLT,MDLT
    CAC=COMMAND_ARGUMENT_COUNT()
    IF(CAC/=2) THEN
        PRINT *,"USAGE lpl_bem_tst surf lbl"
        STOP
    ENDIF
    CALL GET_COMMAND_ARGUMENT(1, ARG)
!    CALL READ_TSURF_B(S,TRIM(ARG),FLAG)
!    IF(FLAG/=0) THEN
!        PRINT *,'ERROR READING TSURF',FLAG
!        STOP
!    ENDIF
    CALL ICOSAHEDRON(S)
    nt=s%nt
    nv=s%nv
    PRINT*,NT,NV
    ALLOCATE(NCVEC(NV),CTVEC(MAX_CON,NV),LM1(NV,NV),QT1(NV),E(NV),EE(NV),XP(NV),CP(3,NV),IPIV(NV),STAT=FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'memory allocation problem',FLAG
        STOP
    ENDIF
    CALL MAKE_CON_TRG(NCVEC,CTVEC,S,FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'ERROR IN MAKE_CON_TRG',FLAG
        STOP
    ENDIF
    CALL MK_QTOT1(QT1,NV,NCVEC,CTVEC,S,FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'ERROR IN MK_QTOT1',FLAG
        STOP
    ENDIF

    CP=S%VVEC
    CALL CPU_TIME(TM0)
    CALL MK_SAMAT_L(LM1,NV,NV,NCVEC,CTVEC,S,CP,DLPL_GF_L1,FLAG)
    CALL CPU_TIME(TM1)
    PRINT *,'MK_SAMAT_L with LPL_GF_L1 TAKES ',TM1-TM0
    CALL CPU_TIME(TM0)
    CALL GETRF_F64(NV,LM1,IPIV,FLAG)
    CALL CPU_TIME(TM1)
    IF(FLAG/=0) THEN
        PRINT *,'ERROR IN GETRF_F64',FLAG
        STOP
    ENDIF
    PRINT *,'GETRF_F64 TAKES ',TM1-TM0
    E=ONE_F64
    PHI=ONE_F64;
    DO I=1,NV
        XP(I)=XPOT(CP(3,I),PHI)
    ENDDO
    CALL GETRS_F64(NV,LM1,1,XP,IPIV,FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'ERROR IN GETRS_F64',FLAG
        STOP
    ENDIF
    CALL GETRS_F64(NV,LM1,1,E,IPIV,FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'ERROR IN GETRS_F64',FLAG
        STOP
    ENDIF
    PHI=DOT_PRODUCT(QT1,XP)/DOT_PRODUCT(QT1,E)
    E=PHI*E-XP
    PRINT*,'QTOT=',DOT_PRODUCT(QT1,E)
    FCT=3.0_F64*ONEOVER4PI_F64
    EE=E-CP(3,:)*FCT
    MDLT=MAXVAL(ABS(EE))
    DLT=SQRT(DOT_PRODUCT(EE,EE)/NV)
    PRINT*,'MAXERR',MDLT,'STDERR',DLT
    CALL GET_COMMAND_ARGUMENT(2, ARG)
    CALL WRITE_ARR_2_F64(TRIM(ARG)//'-LM1.BIN',FLAG,LM1,NV,NV)
    CALL WRITE_ARR_1_F64(TRIM(ARG)//'-EE1.BIN',FLAG,EE,NV)

    IF(ALLOCATED(NCVEC)) DEALLOCATE(NCVEC)
    IF(ALLOCATED(CTVEC)) DEALLOCATE(CTVEC)
    IF(ALLOCATED(LM1)) DEALLOCATE(LM1)
    IF(ALLOCATED(QT1)) DEALLOCATE(QT1)
    IF(ALLOCATED(E)) DEALLOCATE(E)
    IF(ALLOCATED(EE)) DEALLOCATE(EE)
    IF(ALLOCATED(XP)) DEALLOCATE(XP)
    IF(ALLOCATED(CP)) DEALLOCATE(CP)
    IF(ALLOCATED(IPIV)) DEALLOCATE(IPIV)
    CALL CLEAN_TSURF(S)
    
CONTAINS
SUBROUTINE MK_SAMAT_L(LM,M,N,NCVEC,CTVEC,S,CP,TPOT,FLAG)
    REAL(F64)::LM(M,N),TPOT
    REAL(F64)::CP(3,M)
    INTEGER::M,N,NCVEC(N),CTVEC(MAX_CON,N),FLAG
    TYPE(TSURF_TYPE)::S
    EXTERNAL TPOT
    INTEGER::I,J,K,V(3)
    REAL(F64)::Q(3)
    FLAG=1
    IF(N/=S%NV) RETURN ! WRONG PARAMETERS
    LM=ZERO_F64
    DO I=1,N
    DO K=1,NCVEC(I)
    V=S%TVEC(:,CTVEC(K,I))
    Q=ZERO_F64; WHERE(V==I) Q=ONE_F64
    DO J=1,M
        LM(J,I)=LM(J,I)+TPOT(CP(:,J),S%VVEC(:,V(1)),S%VVEC(:,V(2)),S%VVEC(:,V(3)),Q(1),Q(2),Q(3))
    ENDDO
    ENDDO
    ENDDO
    FLAG=0
ENDSUBROUTINE MK_SAMAT_L

REAL(F64) FUNCTION XPOT(Z,PHI)
    REAL(F64)::Z
    REAL(F64)::PHI
    XPOT=-Z+PHI
ENDFUNCTION XPOT

SUBROUTINE MK_QTOT1(Q,N,NCVEC,CTVEC,S,FLAG)
    REAL(F64)::Q(N)
    INTEGER::N,NCVEC(N),CTVEC(MAX_CON,N),FLAG
    TYPE(TSURF_TYPE)::S
    INTEGER::I,V(3),K
    REAL(F64)::A,B(3),C(3),W(3)
    FLAG=1
    IF(N/=S%NV) RETURN ! WRONG PARAMETERS
    DO I=1,N
    A=ZERO_F64
    DO K=1,NCVEC(I)
    V=S%TVEC(:,CTVEC(K,I))
    B=S%VVEC(:,V(2))-S%VVEC(:,V(1))
    C=S%VVEC(:,V(3))-S%VVEC(:,V(1))
    W=CROSS_F64(B,C)
    A=SQRT(DOT_F64(W,W))
    ENDDO
    Q(I)=A/6.0_F64
    ENDDO
    FLAG=0
ENDSUBROUTINE MK_QTOT1
ENDPROGRAM LPL_BEM1D_TST


