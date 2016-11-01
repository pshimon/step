!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM BEM_LAPLACE
    USE LGF
    use tsurf
    IMPLICIT NONE
 
    INTEGER::NS ! NUMBER OF SURFACES
    ! LABELS OF THE SURFACE
    CHARACTER(LEN=16),ALLOCATABLE,DIMENSION(:):: SLABEL ! (NS)
    TYPE(TSURF_TYPE),ALLOCATABLE,DIMENSION(:)::SRF !(NS) SURFACES
    ! (NN OR NT) POTENTIALS AND CHARGES for any surface
    REAL(F32),ALLOCATABLE,DIMENSION(:,:)::FLD     
    INTEGER::NT ! NUMBER OF TRIANGLES 
    INTEGER::NN ! NUMBER OF NODES
    INTEGER::CAC,J,K,I,N1,N2,N3,M
    CHARACTER(LEN=100) :: ARG 
    CHARACTER(LEN=100) :: FNAME   
    INTEGER::FLAG
    REAL(F32)::TM0,TM1,TMA(4) !TIMES 
    INTEGER,PARAMETER::REP_UNIT=6
    INTEGER,PARAMETER::NP=1000
    REAL(F32)::R(3,6*NP),a,b,Q(3)
    REAL(F32)::P(4,6*NP),dlt(3),aa,EPSA,EPSR
    INTEGER::NSUB
    CAC=COMMAND_ARGUMENT_COUNT()
    IF(CAC/=1) THEN
        PRINT *,"USAGE bem-laplace PARFILE"
        STOP
    ENDIF
    CALL GET_COMMAND_ARGUMENT(1, ARG)
    CALL INIT(TRIM(ARG),FLAG)
    PRINT*,NN,NT
    Q=1.0
    R=0.0
    A=1.0
    B=1.0/NP
    DO M=1,NP
        R(1,M)=A+M*B
        R(2,NP+M)=A+M*B
        R(3,2*NP+M)=A+M*B
        R(1,3*NP+M)=-A-M*B
        R(2,4*NP+M)=-A-M*B
        R(3,5*NP+M)=-A-M*B
    ENDDO
    P=0.0
    TMA=0.0
    EPSA=1.0E-6
    EPSR=1.0E-3
    NSUB=10
    DO I=1,4
        CALL CPU_TIME(TM0)
        DO J=1,NS
            DO K=1,SRF(J)%NT
                N1=SRF(J)%TVEC(1,K)
                N2=SRF(J)%TVEC(2,K)
                N3=SRF(J)%TVEC(3,K)
                DO M=1,6*NP
                    P(I,M)=P(I,M)+TST1(R(:,M), &
                                    & SRF(J)%NVEC(:,N1), &
                                    & SRF(J)%NVEC(:,N2), &
                                    & SRF(J)%NVEC(:,N3), &
                                    & Q,I)
                ENDDO
               ENDDO
        ENDDO
        CALL CPU_TIME(TM1)
        TMA(I)=TMA(I)+TM1-TM0
    ENDDO
    PRINT*,TMA
    DO m=1,6*np
        P(2,M)=(P(2,M)-P(1,M))/P(1,M)
        P(3,M)=(P(3,M)-P(1,M))/P(1,M)
        P(4,M)=(P(4,M)-P(1,M))/P(1,M)
    ENDDO
    DLT=0.0
     DO M=1,6*NP
        DO I=2,4
            AA=ABS(P(I,M))
            IF(AA>DLT(I-1)) DLT(I-1)=AA
        ENDDO
    ENDDO
    PRINT*,"MAX DIFF:",DLT
    
CONTAINS
SUBROUTINE READ_PARAMS(FILENAME,FLAG)
    CHARACTER(*) :: FILENAME
    INTEGER,INTENT(OUT)::FLAG
    INTEGER,PARAMETER::IN_UNIT=11
    CHARACTER::B
    CHARACTER(LEN=80) :: LINE
    CHARACTER(LEN=12) :: KEYWORD
    CHARACTER(LEN=60) :: VALUE
    INTEGER::I
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,IOSTAT=FLAG)
    IF(FLAG/=0) RETURN
    DO
        READ(IN_UNIT,'(a80)',IOSTAT=FLAG) LINE
        IF(FLAG/=0)  CYCLE
        IF(LINE(1:1) .NE. '@') CYCLE
        READ(LINE,*) B,KEYWORD, VALUE
        IF ( KEYWORD .EQ. 'begin_bem' ) I=0
        IF ( KEYWORD .EQ. 'nums' ) THEN
            READ( VALUE, * ) NS
            I=I+1
            ALLOCATE(SLABEL(NS),STAT=FLAG)
            IF(FLAG/=0)  RETURN
        ENDIF
        IF ( KEYWORD .EQ. 'surf' ) THEN
            READ( VALUE, * )  SLABEL(I)
            I=I+1
        ENDIF
    
        IF ( KEYWORD .EQ. 'end_bem' ) EXIT
    ENDDO
    FLAG=I-NS-1
    CLOSE(IN_UNIT)
ENDSUBROUTINE READ_PARAMS
SUBROUTINE INIT(FILENAME,FLAG)
    CHARACTER(*) :: FILENAME
    INTEGER,INTENT(OUT)::FLAG
    INTEGER::T,ERR,I
    CHARACTER(LEN=100) :: SURF_FNAME ! FILE DESCRIBING SURFACE (BINARY)

    call READ_PARAMS(FILENAME,ERR)
    IF(ERR/=0) THEN
        FLAG=1
        RETURN
    ENDIF
    ALLOCATE(SRF(NS),STAT=ERR)
    IF(ERR/=0) THEN
        FLAG=2
        RETURN
    ENDIF
    NT=0;NN=0
    DO I=1,NS
        SURF_FNAME=TRIM(SLABEL(I))//'.tsb' 
        CALL READ_TSURF_B(SRF(I),TRIM(SURF_FNAME),ERR)
            IF(ERR/=0) THEN
                FLAG=3
                RETURN
            ENDIF
        NT=NT+SRF(I)%NT
        NN=NT+SRF(I)%NN
    ENDDO
    ALLOCATE(FLD(NT,NS),STAT=ERR)
    IF(ERR/=0) THEN
        FLAG=4
        RETURN
    ENDIF
    FLAG=0
ENDSUBROUTINE INIT
FUNCTION TST1(R,R1,R2,R3,Q,K) RESULT(RES)
    REAL(F32)::R(3),R1(3),R2(3),R3(3),Q(3)
    REAL(F32)::RES
    INTEGER::K
    SELECT CASE(K)
        CASE(1)
            RES=LGFLT(R,R1,R2,R3,Q)
         CASE(2)
            RES=LGFLTN1(R,R1,R2,R3,Q,NSUB)
        CASE(3)
            RES=LGFLTN3(R,R1,R2,R3,Q)
        CASE(4)
            RES=LGFLTN7(R,R1,R2,R3,Q)
        CASE DEFAULT
            PRINT*,'NO SUCH FUNCTION!,EXITING!'
            STOP
    END SELECT
 
ENDFUNCTION TST1

ENDPROGRAM BEM_LAPLACE


