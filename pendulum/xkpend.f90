!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM XKPEND
    USE FDEFS
    USE ARRIO
    IMPLICIT NONE
    INTEGER::M=100,N=10000
    REAL(F32)::ALPHA=0.0,BETA=0.0
    REAL(F32),ALLOCATABLE,DIMENSION(:,:)::LP,RK
    INTEGER::FLAG
    REAL(F32)::PHI0=0.0,PSI0=0.001,T0=0.0
    REAL(F32)::DT=1.0E-4
    REAL(F32)::START_TIME,STOP_TIME

    ALLOCATE(LP(3,N),RK(3,N),STAT=FLAG)
    IF(FLAG.NE.0) THEN
        PRINT *,'PROBLEM ALLOCATING MEM'
        STOP
    ENDIF

    CALL CPU_TIME(START_TIME)
    CALL KP_LP2(LP,M,N,DT,PHI0,PSI0,T0,ALPHA,BETA)
    CALL CPU_TIME(STOP_TIME)
    PRINT *,'LEAPFROG TIME:',STOP_TIME-START_TIME
    CALL WRITE_ARR2_F32(LP,"lp.bin",FLAG);
    IF(FLAG.NE.0) THEN
        PRINT *,'PROBLEM writing file'
        STOP
    ENDIF
    

    CALL CPU_TIME(START_TIME)
    CALL KP_RK4(RK,M,N,DT,PHI0,PSI0,T0,ALPHA,BETA)
    CALL CPU_TIME(STOP_TIME)
    PRINT *,'RUNGE-KUTTA TIME:',STOP_TIME-START_TIME
    CALL WRITE_ARR2_F32(RK,"rk.bin",FLAG);
   IF(FLAG.NE.0) THEN
        PRINT *,'PROBLEM writing file'
        STOP
    ENDIF
    
ENDPROGRAM XKPEND 

