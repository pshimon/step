!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!! FLOAT ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE ARRIO
    USE FDEFS
    IMPLICIT NONE
    PUBLIC
CONTAINS

!   1D ARRAYS
SUBROUTINE  READ_ARR1_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),ALLOCATABLE,DIMENSION(:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1 
    ALLOCATE(A(N1),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR1_F32

SUBROUTINE  WRITE_ARR1_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),DIMENSION(:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR1_F32

!   2D ARRAYS
SUBROUTINE  READ_ARR2_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2 
    ALLOCATE(A(N1,N2),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR2_F32

SUBROUTINE  WRITE_ARR2_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),DIMENSION(:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR2_F32

!   3D ARRAYS
SUBROUTINE  READ_ARR3_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3 
    ALLOCATE(A(N1,N2,N3),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR3_F32

SUBROUTINE  WRITE_ARR3_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),DIMENSION(:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR3_F32

!   4D ARRAYS
SUBROUTINE  READ_ARR4_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3,N4
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3,N4 
    ALLOCATE(A(N1,N2,N3,N4),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR4_F32

SUBROUTINE  WRITE_ARR4_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),DIMENSION(:,:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3,N4
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    N4=SIZE(A,4)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3,N4
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR4_F32

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   DOUBLE ARRAYS!!!!!!!!!!!!!!!!!!!!!!!

!   1D ARRAYS
SUBROUTINE  READ_ARR1_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),ALLOCATABLE,DIMENSION(:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1 
    ALLOCATE(A(N1),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR1_F64

SUBROUTINE  WRITE_ARR1_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),DIMENSION(:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR1_F64

!   2D ARRAYS
SUBROUTINE  READ_ARR2_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2 
    ALLOCATE(A(N1,N2),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR2_F64

SUBROUTINE  WRITE_ARR2_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),DIMENSION(:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR2_F64

!   3D ARRAYS
SUBROUTINE  READ_ARR3_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3 
    ALLOCATE(A(N1,N2,N3),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR3_F64

SUBROUTINE  WRITE_ARR3_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),DIMENSION(:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR3_F64

!   4D ARRAYS
SUBROUTINE  READ_ARR4_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3,N4
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3,N4 
    ALLOCATE(A(N1,N2,N3,N4),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR4_F64

SUBROUTINE  WRITE_ARR4_F64(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F64),DIMENSION(:,:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3,N4
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    N4=SIZE(A,4)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3,N4
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR4_F64

!!!!!!!!!!!!!!!!!!!!!!!!!! COMPLEX FLOAT ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!

!   1D ARRAYS
SUBROUTINE  READ_ARR1_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),ALLOCATABLE,DIMENSION(:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1 
    ALLOCATE(A(N1),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR1_C32

SUBROUTINE  WRITE_ARR1_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),DIMENSION(:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR1_C32

!   2D ARRAYS
SUBROUTINE  READ_ARR2_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2 
    ALLOCATE(A(N1,N2),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR2_C32

SUBROUTINE  WRITE_ARR2_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),DIMENSION(:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR2_C32

!   3D ARRAYS
SUBROUTINE  READ_ARR3_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3 
    ALLOCATE(A(N1,N2,N3),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR3_C32

SUBROUTINE  WRITE_ARR3_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),DIMENSION(:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR3_C32

!   4D ARRAYS
SUBROUTINE  READ_ARR4_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3,N4
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3,N4 
    ALLOCATE(A(N1,N2,N3,N4),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR4_C32

SUBROUTINE  WRITE_ARR4_C32(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F32),DIMENSION(:,:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3,N4
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    N4=SIZE(A,4)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3,N4
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR4_C32

!!!!!!!!!!!!!!!!!COMPLEX   DOUBLE ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   1D ARRAYS
SUBROUTINE  READ_ARR1_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),ALLOCATABLE,DIMENSION(:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1 
    ALLOCATE(A(N1),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR1_C64

SUBROUTINE  WRITE_ARR1_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),DIMENSION(:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR1_C64

!   2D ARRAYS
SUBROUTINE  READ_ARR2_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2 
    ALLOCATE(A(N1,N2),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR2_C64

SUBROUTINE  WRITE_ARR2_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),DIMENSION(:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR2_C64

!   3D ARRAYS
SUBROUTINE  READ_ARR3_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3 
    ALLOCATE(A(N1,N2,N3),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR3_C64

SUBROUTINE  WRITE_ARR3_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),DIMENSION(:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR3_C64

!   4D ARRAYS
SUBROUTINE  READ_ARR4_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3,N4
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3,N4 
    ALLOCATE(A(N1,N2,N3,N4),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR4_C64

SUBROUTINE  WRITE_ARR4_C64(A,FILENAME,FLAG)
    USE FDEFS
    COMPLEX(F64),DIMENSION(:,:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3,N4
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    N4=SIZE(A,4)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3,N4
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR4_C64


!!!!!!!!!!!!!!!!!!!!!!!!!! INTEGER ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   1D ARRAYS
SUBROUTINE  READ_ARR1_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1 
    ALLOCATE(A(N1),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR1_INT

SUBROUTINE  WRITE_ARR1_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,DIMENSION(:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR1_INT

!   2D ARRAYS
SUBROUTINE  READ_ARR2_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2 
    ALLOCATE(A(N1,N2),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR2_INT

SUBROUTINE  WRITE_ARR2_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,DIMENSION(:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR2_INT

!   3D ARRAYS
SUBROUTINE  READ_ARR3_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3 
    ALLOCATE(A(N1,N2,N3),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR3_INT

SUBROUTINE  WRITE_ARR3_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,DIMENSION(:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR3_INT

!   4D ARRAYS
SUBROUTINE  READ_ARR4_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,N1,N2,N3,N4
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) N1,N2,N3,N4 
    ALLOCATE(A(N1,N2,N3,N4),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR4_INT

SUBROUTINE  WRITE_ARR4_INT(A,FILENAME,FLAG)
    USE FDEFS
    INTEGER,DIMENSION(:,:,:,:),INTENT(IN) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,N1,N2,N3,N4
    PARAMETER(OUT_UNIT=12)
    N1=SIZE(A,1)
    N2=SIZE(A,2)
    N3=SIZE(A,3)
    N4=SIZE(A,4)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) N1,N2,N3,N4
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR4_INT
!   TEXT 2D ARRAYS
SUBROUTINE  READ_TXT_ARR2_F32(A,FILENAME,FLAG)
    USE FDEFS
    REAL(F32),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,n1,n2,i,j
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,*) N1,N2 
    ALLOCATE(A(N1,N2),STAT=FLAG)
    IF(FLAG/=0) RETURN
    do i=1,n1
        READ(IN_UNIT,*) (A(i,j),j=1,n2)
    enddo
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_TXT_ARR2_F32

end MODULE ARRIO


