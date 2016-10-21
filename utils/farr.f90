!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!! FLOAT ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE FARR
    USE FDEFS
    IMPLICIT NONE
    PUBLIC
    CONTAINS

SUBROUTINE  READ_ARR7_F32(A,FILENAME,FLAG)
    REAL(F32),ALLOCATABLE,DIMENSION(:,:,:,:,:,:,:) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER ::IN_UNIT,R,SHP(7)
    PARAMETER(IN_UNIT=11)
    OPEN(UNIT=IN_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="OLD")
    IF(FLAG/=0) RETURN
    READ(IN_UNIT) R,SHP
    ALLOCATE(A(SHP(1),SHP(2),SHP(3),SHP(4),SHP(5),SHP(6),SHP(7)),STAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(IN_UNIT,IOSTAT=FLAG) A
    IF(FLAG/=0) THEN
            DEALLOCATE(A)
            RETURN
    END IF
    CLOSE(IN_UNIT)
END SUBROUTINE  READ_ARR7_F32

SUBROUTINE  WRITE_ARR7_F32(A,FILENAME,FLAG)
    REAL(F32),DIMENSION(:,:,:,:,:,:,:) ::A
    INTEGER,INTENT(OUT) ::FLAG
    CHARACTER(*) :: FILENAME
    INTEGER::OUT_UNIT,R
    PARAMETER(OUT_UNIT=12,R=7)
    OPEN(UNIT=OUT_UNIT,FILE=FILENAME,ACCESS="STREAM",&
    & IOSTAT=FLAG,STATUS="REPLACE")
    IF(FLAG/=0) RETURN
    WRITE(OUT_UNIT) R,SHAPE(A)
    WRITE(OUT_UNIT,IOSTAT=FLAG) A
    CLOSE(OUT_UNIT)
END SUBROUTINE  WRITE_ARR7_F32
ENDMODULE FARR
