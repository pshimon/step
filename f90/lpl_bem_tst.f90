!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM LPL_BEM_TST
    USE T_SURF
    IMPLICIT NONE
    TYPE(TSURF_TYPE)::S     
    INTEGER::NT ! NUMBER OF TRIANGLES 
    INTEGER::NV ! NUMBER OF VERTICES
    REAL(F64),ALLOCATABLE,DIMENSION(:,:)::LM1,LM0,CP !
    INTEGER,ALLOCATABLE::NCVEC(:),CTVEC(:,:)
    INTEGER::CAC
    CHARACTER(LEN=100) :: ARG    
    INTEGER::FLAG
    REAL(F32)::TM0,TM1
    REAL(F64)::DLT,LPL_GF_L1 !,lmmax,lmmin
    EXTERNAL LPL_GF_L1

    INTERFACE

    FUNCTION LPL_GF_L2(DST,VRT0,VRT1,VRT2,Q0,Q1,Q2) RESULT (RES)
    IMPORT::F32,F64
    REAL(F32)::DST(3),VRT0(3),VRT1(3),VRT2(3)
    REAL(F64)::Q1,Q2,Q0,RES
    ENDFUNCTION LPL_GF_L2
    FUNCTION lplGfL1(DST,VRT0,VRT1,VRT2,Q0,Q1,Q2) RESULT (RES) bind(c,name='lplGfL1')
    IMPORT::C_FLOAT,C_DOUBLE
    REAL(C_FLOAT)::DST(3),VRT0(3),VRT1(3),VRT2(3)
    REAL(C_DOUBLE)::Q1,Q2,Q0,RES
    ENDFUNCTION lplGfL1

    ENDINTERFACE


    CAC=COMMAND_ARGUMENT_COUNT()
    IF(CAC/=1) THEN
        PRINT *,"USAGE lpl_bem_tst surf"
        STOP
    ENDIF
    CALL GET_COMMAND_ARGUMENT(1, ARG)
    CALL READ_TSURF_B(S,TRIM(ARG),FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'ERROR READING TSURF',FLAG
        STOP
    ENDIF
    nt=s%nt
    nv=s%nv
    PRINT*,NT,NV
    ALLOCATE(NCVEC(NV),CTVEC(MAX_CON,NV),LM1(NV,NV),LM0(NV,NV),CP(3,NV),STAT=FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'memory allocation problem',FLAG
        STOP
    ENDIF
    CALL MAKE_CON_TRG(NCVEC,CTVEC,S,FLAG)
    IF(FLAG/=0) THEN
        PRINT *,'ERROR IN MAKE_CON_TRG',FLAG
        STOP
    ENDIF
    CP=S%VVEC
    CALL CPU_TIME(TM0)
    CALL MK_SAMAT_L(LM1,NV,NV,NCVEC,CTVEC,S,CP,LPL_GF_L1,FLAG)
    CALL CPU_TIME(TM1)
    PRINT *,'LPL_GF_L1 TAKES ',TM1-TM0
    !lmmax=maxval(lm1)
    !lmmin=minval(lm1)
    !print *,'lmmin',lmmin,'lmmax',lmmax
    CALL CPU_TIME(TM0)
    CALL MK_SAMAT_L(LM0,NV,NV,NCVEC,CTVEC,S,CP,LPL_GF_L2,FLAG)
!    CALL MK_SAMAT_L(LM0,NV,NV,NCVEC,CTVEC,S,CP,lplGfL1,FLAG)
    CALL CPU_TIME(TM1)
    PRINT *,'LPL_GF_L2 TAKES ',TM1-TM0
!    PRINT *,'lplGfL1 TAKES ',TM1-TM0
    !lmmax=maxval(lm0)
    !lmmin=minval(lm0)
    !print *,'lmmin',lmmin,'lmmax',lmmax

    DLT=maxval(abs(LM1-LM0))
    print *,'maxdif',dlt
    IF(ALLOCATED(NCVEC)) DEALLOCATE(NCVEC)
    IF(ALLOCATED(CTVEC)) DEALLOCATE(CTVEC)
    IF(ALLOCATED(LM1)) DEALLOCATE(LM1)
    IF(ALLOCATED(LM0)) DEALLOCATE(LM0)
    IF(ALLOCATED(CP)) DEALLOCATE(CP)
    CALL CLEAN_TSURF(S)
    
CONTAINS
SUBROUTINE MK_SAMAT_L(LM,M,N,NCVEC,CTVEC,S,CP,TPOT,FLAG)
    REAL(F64)::LM(M,N),TPOT,CP(3,M)
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

ENDPROGRAM LPL_BEM_TST


