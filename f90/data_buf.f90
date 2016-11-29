!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!! FLOAT ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE DATA_BUF
    USE COMMON_DEFS
    IMPLICIT NONE
    PUBLIC
    INTEGER,PARAMETER::LBL_I08=0,LBL_F32=16,LBL_F64=32,LBL_C32=48,LBL_C64=64,LBL_I32=80
    INTEGER,PARAMETER::HDR_SIZE=8
    TYPE :: CDATA_F32 
        REAL(F32),ALLOCATABLE,DIMENSION(:)::BUF  
        INTEGER::HDR(HDR_SIZE)
    END TYPE  CDATA_F32 

    TYPE :: CDATA_F64 
        REAL(F64),ALLOCATABLE,DIMENSION(:)::BUF  
        INTEGER::HDR(HDR_SIZE)
    END TYPE  CDATA_F64 
  
CONTAINS

    
LOGICAL FUNCTION CMP_HDR(H1,H2)
    INTEGER::H1(HDR_SIZE),H2(HDR_SIZE)
    CMP_HDR=SUM(ABS(H1-H2))==0
ENDFUNCTION CMP_HDR

SUBROUTINE READ_HDR(FNAME,FLAG,H)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG,H(HDR_SIZE)
    INTEGER::RDUNITB
    OPEN(NEWUNIT=RDUNITB,FILE=FNAME,ACCESS="STREAM",&
        STATUS="OLD",IOSTAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(RDUNITB) H
    CLOSE(RDUNITB)
    FLAG=0
ENDSUBROUTINE READ_HDR


! float data    
SUBROUTINE MK_HDR_1_F32(H,N1)
    INTEGER::H(HDR_SIZE),N1
    H=1
    H(1)=LBL_F32+1
    H(2)=N1
ENDSUBROUTINE MK_HDR_1_F32
SUBROUTINE MK_HDR_2_F32(H,N1,N2)
    INTEGER::H(HDR_SIZE),N1,N2
    H=1
    H(1)=LBL_F32+2
    H(2)=N1
    H(3)=N2
ENDSUBROUTINE MK_HDR_2_F32
SUBROUTINE MK_HDR_3_F32(H,N1,N2,N3)
    INTEGER::H(HDR_SIZE),N1,N2,N3
    H=1
    H(1)=LBL_F32+3
    H(2)=N1
    H(3)=N2
    H(4)=N3
ENDSUBROUTINE MK_HDR_3_F32
SUBROUTINE MK_HDR_4_F32(H,N1,N2,N3,N4)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4
    H=1
    H(1)=LBL_F32+4
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
ENDSUBROUTINE MK_HDR_4_F32
SUBROUTINE MK_HDR_5_F32(H,N1,N2,N3,N4,N5)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4,N5
    H=1
    H(1)=LBL_F32+5
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
    H(6)=N5
ENDSUBROUTINE MK_HDR_5_F32
SUBROUTINE MK_HDR_6_F32(H,N1,N2,N3,N4,N5,N6)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4,N5,N6
    H=1
    H(1)=LBL_F32+6
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
    H(6)=N5
    H(7)=N6
ENDSUBROUTINE MK_HDR_6_F32
SUBROUTINE MK_HDR_7_F32(H,N1,N2,N3,N4,N5,N6,N7)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4,N5,N6,N7
    H=1
    H(1)=LBL_F32+7
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
    H(6)=N5
    H(7)=N6
    H(8)=N7
ENDSUBROUTINE MK_HDR_7_F32


SUBROUTINE READ_ARR_7_F32(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6,N7)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4,N5,N6,N7)
    INTEGER::N1,N2,N3,N4,N5,N6,N7
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_7_F32(H1,N1,N2,N3,N4,N5,N6,N7)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_7_F32

SUBROUTINE WRITE_ARR_7_F32(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6,N7)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4,N5,N6,N7)
    INTEGER::N1,N2,N3,N4,N5,N6,N7
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_7_F32(H1,N1,N2,N3,N4,N5,N6,N7)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_7_F32

SUBROUTINE READ_ARR_6_F32(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4,N5,N6)
    INTEGER::N1,N2,N3,N4,N5,N6
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_6_F32(H1,N1,N2,N3,N4,N5,N6)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_6_F32

SUBROUTINE WRITE_ARR_6_F32(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4,N5,N6)
    INTEGER::N1,N2,N3,N4,N5,N6
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_6_F32(H1,N1,N2,N3,N4,N5,N6)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_6_F32

SUBROUTINE READ_ARR_5_F32(FNAME,FLAG,A,N1,N2,N3,N4,N5)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4,N5)
    INTEGER::N1,N2,N3,N4,N5
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_5_F32(H1,N1,N2,N3,N4,N5)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_5_F32

SUBROUTINE WRITE_ARR_5_F32(FNAME,FLAG,A,N1,N2,N3,N4,N5)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4,N5)
    INTEGER::N1,N2,N3,N4,N5
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_5_F32(H1,N1,N2,N3,N4,N5)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_5_F32

SUBROUTINE READ_ARR_4_F32(FNAME,FLAG,A,N1,N2,N3,N4)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4)
    INTEGER::N1,N2,N3,N4
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_4_F32(H1,N1,N2,N3,N4)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_4_F32

SUBROUTINE WRITE_ARR_4_F32(FNAME,FLAG,A,N1,N2,N3,N4)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3,N4)
    INTEGER::N1,N2,N3,N4
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_4_F32(H1,N1,N2,N3,N4)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_4_F32

SUBROUTINE READ_ARR_3_F32(FNAME,FLAG,A,N1,N2,N3)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3)
    INTEGER::N1,N2,N3
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_3_F32(H1,N1,N2,N3)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_3_F32

SUBROUTINE WRITE_ARR_3_F32(FNAME,FLAG,A,N1,N2,N3)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2,N3)
    INTEGER::N1,N2,N3
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_3_F32(H1,N1,N2,N3)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_3_F32

SUBROUTINE READ_ARR_2_F32(FNAME,FLAG,A,N1,N2)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2)
    INTEGER::N1,N2
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_2_F32(H1,N1,N2)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_2_F32

SUBROUTINE WRITE_ARR_2_F32(FNAME,FLAG,A,N1,N2)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1,N2)
    INTEGER::N1,N2
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_2_F32(H1,N1,N2)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_2_F32

SUBROUTINE READ_ARR_1_F32(FNAME,FLAG,A,N1)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1)
    INTEGER::N1
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_1_F32(H1,N1)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_1_F32

SUBROUTINE WRITE_ARR_1_F32(FNAME,FLAG,A,N1)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F32)::A(N1)
    INTEGER::N1
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_1_F32(H1,N1)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_1_F32
SUBROUTINE CLEAN_CDATA_F32(D)
    TYPE(CDATA_F32),INTENT(INOUT) :: D
    IF(ALLOCATED(D%BUF)) DEALLOCATE(D%BUF)
    D%HDR=0
END SUBROUTINE CLEAN_CDATA_F32

SUBROUTINE WRITE_CDATA_F32(D,FNAME,FLAG)
    TYPE(CDATA_F32),INTENT(INOUT):: D
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    INTEGER::WRUNITB
    OPEN(NEWUNIT=WRUNITB,FILE=FNAME,ACCESS="STREAM",&
        STATUS="REPLACE",IOSTAT=FLAG)
    IF(FLAG>0) RETURN
    WRITE(WRUNITB,IOSTAT=FLAG) D%HDR,D%BUF
    CLOSE(WRUNITB)
END SUBROUTINE WRITE_CDATA_F32

SUBROUTINE READ_CDATA_F32(D,FNAME,FLAG)
    TYPE(CDATA_F32),INTENT(INOUT):: D    
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    INTEGER::RDUNITB,N,I
    OPEN(NEWUNIT=RDUNITB,FILE=FNAME,ACCESS="STREAM",&
        STATUS="OLD",IOSTAT=FLAG)
    IF(FLAG/=0) RETURN
    CALL CLEAN_CDATA_F32(D)
    READ(RDUNITB) D%HDR
    N=1 ! COMPUTE TOTAL SIZE
    DO I=2,HDR_SIZE
        N=N*D%HDR(I);
    ENDDO
    ALLOCATE(D%BUF(N),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(RDUNITB,IOSTAT=FLAG) D%BUF
    CLOSE(RDUNITB)    
    
END SUBROUTINE READ_CDATA_F32
        
! double data

SUBROUTINE MK_HDR_1_F64(H,N1)
    INTEGER::H(HDR_SIZE),N1
    H=1
    H(1)=LBL_F64+1
    H(2)=N1
ENDSUBROUTINE MK_HDR_1_F64
SUBROUTINE MK_HDR_2_F64(H,N1,N2)
    INTEGER::H(HDR_SIZE),N1,N2
    H=1
    H(1)=LBL_F64+2
    H(2)=N1
    H(3)=N2
ENDSUBROUTINE MK_HDR_2_F64
SUBROUTINE MK_HDR_3_F64(H,N1,N2,N3)
    INTEGER::H(HDR_SIZE),N1,N2,N3
    H=1
    H(1)=LBL_F64+3
    H(2)=N1
    H(3)=N2
    H(4)=N3
ENDSUBROUTINE MK_HDR_3_F64
SUBROUTINE MK_HDR_4_F64(H,N1,N2,N3,N4)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4
    H=1
    H(1)=LBL_F64+4
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
ENDSUBROUTINE MK_HDR_4_F64
SUBROUTINE MK_HDR_5_F64(H,N1,N2,N3,N4,N5)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4,N5
    H=1
    H(1)=LBL_F64+5
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
    H(6)=N5
ENDSUBROUTINE MK_HDR_5_F64
SUBROUTINE MK_HDR_6_F64(H,N1,N2,N3,N4,N5,N6)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4,N5,N6
    H=1
    H(1)=LBL_F64+6
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
    H(6)=N5
    H(7)=N6
ENDSUBROUTINE MK_HDR_6_F64
SUBROUTINE MK_HDR_7_F64(H,N1,N2,N3,N4,N5,N6,N7)
    INTEGER::H(HDR_SIZE),N1,N2,N3,N4,N5,N6,N7
    H=1
    H(1)=LBL_F64+7
    H(2)=N1
    H(3)=N2
    H(4)=N3
    H(5)=N4
    H(6)=N5
    H(7)=N6
    H(8)=N7
ENDSUBROUTINE MK_HDR_7_F64


SUBROUTINE READ_ARR_7_F64(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6,N7)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4,N5,N6,N7)
    INTEGER::N1,N2,N3,N4,N5,N6,N7
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_7_F64(H1,N1,N2,N3,N4,N5,N6,N7)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_7_F64

SUBROUTINE WRITE_ARR_7_F64(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6,N7)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4,N5,N6,N7)
    INTEGER::N1,N2,N3,N4,N5,N6,N7
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_7_F64(H1,N1,N2,N3,N4,N5,N6,N7)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_7_F64

SUBROUTINE READ_ARR_6_F64(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4,N5,N6)
    INTEGER::N1,N2,N3,N4,N5,N6
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_6_F64(H1,N1,N2,N3,N4,N5,N6)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_6_F64

SUBROUTINE WRITE_ARR_6_F64(FNAME,FLAG,A,N1,N2,N3,N4,N5,N6)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4,N5,N6)
    INTEGER::N1,N2,N3,N4,N5,N6
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_6_F64(H1,N1,N2,N3,N4,N5,N6)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_6_F64
SUBROUTINE READ_ARR_5_F64(FNAME,FLAG,A,N1,N2,N3,N4,N5)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4,N5)
    INTEGER::N1,N2,N3,N4,N5
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_5_F64(H1,N1,N2,N3,N4,N5)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_5_F64

SUBROUTINE WRITE_ARR_5_F64(FNAME,FLAG,A,N1,N2,N3,N4,N5)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4,N5)
    INTEGER::N1,N2,N3,N4,N5
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_5_F64(H1,N1,N2,N3,N4,N5)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_5_F64

SUBROUTINE READ_ARR_4_F64(FNAME,FLAG,A,N1,N2,N3,N4)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4)
    INTEGER::N1,N2,N3,N4
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_4_F64(H1,N1,N2,N3,N4)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_4_F64

SUBROUTINE WRITE_ARR_4_F64(FNAME,FLAG,A,N1,N2,N3,N4)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3,N4)
    INTEGER::N1,N2,N3,N4
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_4_F64(H1,N1,N2,N3,N4)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_4_F64

SUBROUTINE READ_ARR_3_F64(FNAME,FLAG,A,N1,N2,N3)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3)
    INTEGER::N1,N2,N3
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_3_F64(H1,N1,N2,N3)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_3_F64

SUBROUTINE WRITE_ARR_3_F64(FNAME,FLAG,A,N1,N2,N3)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2,N3)
    INTEGER::N1,N2,N3
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_3_F64(H1,N1,N2,N3)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_3_F64

SUBROUTINE READ_ARR_2_F64(FNAME,FLAG,A,N1,N2)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2)
    INTEGER::N1,N2
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_2_F64(H1,N1,N2)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_2_F64

SUBROUTINE WRITE_ARR_2_F64(FNAME,FLAG,A,N1,N2)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1,N2)
    INTEGER::N1,N2
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_2_F64(H1,N1,N2)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_2_F64

SUBROUTINE READ_ARR_1_F64(FNAME,FLAG,A,N1)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1)
    INTEGER::N1
    INTEGER::UNITB,H1(HDR_SIZE),H2(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_1_F64(H1,N1)

    READ(UNITB,IOSTAT=FLAG) H2
    IF(FLAG/=0) THEN
        CLOSE(UNITB)
        RETURN
    ENDIF
    IF(CMP_HDR(H1,H2)) THEN
        READ(UNITB,IOSTAT=FLAG) A
        CLOSE(UNITB)
        RETURN
    ELSE
        FLAG=-1
        CLOSE(UNITB)
        RETURN
    ENDIF
ENDSUBROUTINE READ_ARR_1_F64

SUBROUTINE WRITE_ARR_1_F64(FNAME,FLAG,A,N1)
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    REAL(F64)::A(N1)
    INTEGER::N1
    INTEGER::UNITB,H1(HDR_SIZE)
    OPEN(NEWUNIT=UNITB,FILE=FNAME,ACCESS="STREAM",IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    CALL MK_HDR_1_F64(H1,N1)
    WRITE(UNITB,IOSTAT=FLAG) H1,A
    CLOSE(UNITB)
ENDSUBROUTINE WRITE_ARR_1_F64
SUBROUTINE CLEAN_CDATA_F64(D)
    TYPE(CDATA_F64),INTENT(INOUT) :: D
    IF(ALLOCATED(D%BUF)) DEALLOCATE(D%BUF)
    D%HDR=0
END SUBROUTINE CLEAN_CDATA_F64

SUBROUTINE WRITE_CDATA_F64(D,FNAME,FLAG)
    TYPE(CDATA_F64),INTENT(INOUT):: D
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    INTEGER::WRUNITB
    OPEN(NEWUNIT=WRUNITB,FILE=FNAME,ACCESS="STREAM",&
        STATUS="REPLACE",IOSTAT=FLAG)
    IF(FLAG>0) RETURN
    WRITE(WRUNITB,IOSTAT=FLAG) D%HDR,D%BUF
    CLOSE(WRUNITB)
END SUBROUTINE WRITE_CDATA_F64

SUBROUTINE READ_CDATA_F64(D,FNAME,FLAG)
    TYPE(CDATA_F64),INTENT(INOUT):: D    
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    INTEGER::RDUNITB,N,I
    OPEN(NEWUNIT=RDUNITB,FILE=FNAME,ACCESS="STREAM",&
        STATUS="OLD",IOSTAT=FLAG)
    IF(FLAG/=0) RETURN
    CALL CLEAN_CDATA_F64(D)
    READ(RDUNITB) D%HDR
    N=1 ! COMPUTE TOTAL SIZE
    DO I=2,HDR_SIZE
        N=N*D%HDR(I);
    ENDDO
    ALLOCATE(D%BUF(N),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(RDUNITB,IOSTAT=FLAG) D%BUF
    CLOSE(RDUNITB)    
    
END SUBROUTINE READ_CDATA_F64

ENDMODULE DATA_BUF
