!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SHIMON PANFIL: INDUSTRIAL PHYSICS AND SIMULATIONS        !
! HTTP://INDUSTRIALPHYS.COM                                !
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE TSURF
USE FDEFS
USE VEC3D
IMPLICIT NONE
PUBLIC

TYPE :: TSURF_TYPE ! TRIANGULATED SURFACE
    INTEGER,ALLOCATABLE,DIMENSION(:,:)::TVEC   ! TRIANGLES (3,NT) 
    REAL(F32),ALLOCATABLE,DIMENSION(:,:)::VVEC ! VERTICES  (3,NV) 
    REAL(F32),ALLOCATABLE,DIMENSION(:,:)::NVEC ! NORMALS (3,NV) 
    INTEGER::NT ! NUMBER OF TRIANGLES
    INTEGER::NV ! NUMBER OF VERTICES
END TYPE TSURF_TYPE  

REAL(F32),PARAMETER :: MIN_NODE_DIST=1.0E-4
REAL(F32),PARAMETER :: AREA_ZERO=1.0E-12
REAL(F32),PARAMETER :: MIN_NORM=1.0E-12
INTEGER,PARAMETER::MAX_CON=8

CONTAINS

SUBROUTINE CLEAN_TSURF(S)
    TYPE(TSURF_TYPE),INTENT(INOUT) :: S
    IF(ALLOCATED(S%VVEC)) DEALLOCATE(S%VVEC)
    IF(ALLOCATED(S%TVEC)) DEALLOCATE(S%TVEC)
    IF(ALLOCATED(S%NVEC)) DEALLOCATE(S%NVEC)
    S%NV=0;S%NT=0
END SUBROUTINE CLEAN_TSURF

SUBROUTINE MAKE_TSURF(S,N,T,ERR)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    INTEGER,INTENT(IN):: N,T
    INTEGER,INTENT(OUT):: ERR
    IF((T<1).OR.(N<3)) THEN ! DO NOTHING
        ERR=117
        RETURN
    ENDIF
    CALL CLEAN_TSURF(S)
    S%NV=N;S%NT=T
    ALLOCATE(S%VVEC(3,S%NV),S%NVEC(3,S%NV),S%TVEC(3,S%NT),STAT=ERR)
END SUBROUTINE MAKE_TSURF

SUBROUTINE WRITE_TSURF_A(S,FNAME)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG,I
    INTEGER,PARAMETER::WRUNITA=11
    OPEN(UNIT=WRUNITA,FILE=FNAME,ACCESS="SEQUENTIAL",&
        STATUS="REPLACE",FORM="FORMATTED",IOSTAT=FLAG)
    IF(FLAG>0) RETURN
    WRITE(WRUNITA,*) S%NT,S%NV
    S%TVEC=S%TVEC-1 ! C INDICES
    DO I=1,S%NT
        WRITE(WRUNITA,*) S%TVEC(1,I),S%TVEC(2,I),S%TVEC(3,I)
    ENDDO
    S%TVEC=S%TVEC+1 ! FORTRAN INDICES  
    DO I=1,S%NV
        WRITE(WRUNITA,*) S%VVEC(1,I),S%VVEC(2,I),S%VVEC(3,I)
    ENDDO
    DO I=1,S%NV
        WRITE(WRUNITA,*) S%NVEC(1,I),S%NVEC(2,I),S%NVEC(3,I)
    ENDDO
    CLOSE(WRUNITA)
END SUBROUTINE WRITE_TSURF_A

SUBROUTINE READ_TSURF_A(S,FNAME)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG,I,N,T
    INTEGER,PARAMETER::RDUNITA=12
    OPEN(UNIT=RDUNITA,FILE=FNAME,ACCESS="SEQUENTIAL",&
        STATUS="OLD",IOSTAT=FLAG)
    IF(FLAG>0) RETURN
    READ(RDUNITA,*) T,N
    CALL MAKE_TSURF(S,N,T,FLAG)
    IF(FLAG/=0) RETURN
    DO I=1,S%NT
        READ(RDUNITA,*) S%TVEC(1,I),S%TVEC(2,I),S%TVEC(3,I)
    ENDDO
    S%TVEC=S%TVEC+1 ! FORTRAN INDICES
    DO I=1,S%NV
        READ(RDUNITA,*) S%VVEC(1,I),S%VVEC(2,I),S%VVEC(3,I)
    ENDDO
    DO I=1,S%NV
        READ(RDUNITA,*) S%NVEC(1,I),S%NVEC(2,I),S%NVEC(3,I)
    ENDDO
    CLOSE(RDUNITA)
END SUBROUTINE READ_TSURF_A

SUBROUTINE WRITE_TSURF_B(S,FNAME)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG
    INTEGER,PARAMETER::WRUNITB=12
    OPEN(UNIT=WRUNITB,FILE=FNAME,ACCESS="STREAM",&
        STATUS="REPLACE",IOSTAT=FLAG)
    IF(FLAG>0) RETURN
    WRITE(WRUNITB) S%NT,S%NV
    S%TVEC=S%TVEC-1 ! C INDICES
    WRITE(WRUNITB) S%TVEC
    S%TVEC=S%TVEC+1 ! FORTRAN INDICES
    WRITE(WRUNITB) S%VVEC
    WRITE(WRUNITB) S%NVEC    
    CLOSE(WRUNITB)
END SUBROUTINE WRITE_TSURF_B

SUBROUTINE READ_TSURF_B(S,FNAME,FLAG)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    CHARACTER(*) :: FNAME
    INTEGER :: FLAG,N,T
    INTEGER,PARAMETER::RDUNITB=11
    OPEN(UNIT=RDUNITB,FILE=FNAME,ACCESS="STREAM",&
        STATUS="OLD",IOSTAT=FLAG)
    IF(FLAG/=0) RETURN
    READ(RDUNITB) T,N
    CALL MAKE_TSURF(S,N,T,FLAG)
    IF(FLAG/=0) RETURN
    READ(RDUNITB) S%TVEC
    S%TVEC=S%TVEC+1 ! FORTRAN INDICES
    READ(RDUNITB) S%VVEC
    READ(RDUNITB) S%NVEC
    CLOSE(RDUNITB)    
    FLAG=0
END SUBROUTINE READ_TSURF_B


SUBROUTINE MAKE_CON_TRG(NCVEC,CTVEC,S,FLAG)
    INTEGER,DIMENSION(:),INTENT(INOUT)::NCVEC
    INTEGER,DIMENSION(MAX_CON,*),INTENT(INOUT)::CTVEC
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    INTEGER,INTENT(OUT)::FLAG
    INTEGER::J,K,N
    NCVEC(1:S%NV)=0
    CTVEC(:,1:S%NV)=0
    DO J=1,S%NT
        DO K=1,3
            N=S%TVEC(K,J) !NODE
            NCVEC(N)=NCVEC(N)+1 
            IF(NCVEC(N)>MAX_CON) THEN
                FLAG=J
                RETURN
            ENDIF
            CTVEC(NCVEC(N),N)=J
        ENDDO
    ENDDO
    FLAG=0
ENDSUBROUTINE MAKE_CON_TRG

SUBROUTINE MAKE_UNIT_SPHERE(S)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    REAL(F32)::X,Y,Z,R,X0,Y0,Z0
    INTEGER::J

    X0=0.0;Y0=0.0;Z0=0.0;
    DO J=1,S%NV
        X0=X0+S%VVEC(1,J)
        Y0=Y0+S%VVEC(2,J)
        Z0=Z0+S%VVEC(3,J)
    ENDDO
    X0=X0/S%NV;Y0=Y0/S%NV;Z0=Z0/S%NV !CENTER
    DO J=1,S%NV
        X=S%VVEC(1,J)-X0;Y=S%VVEC(2,J)-Y0;Z=S%VVEC(3,J)-Z0
        R=SQRT(X**2+Y**2+Z**2)
        S%VVEC(1,J)=X/R+X0
        S%VVEC(2,J)=Y/R+Y0
        S%VVEC(3,J)=Z/R+Z0
    ENDDO
ENDSUBROUTINE MAKE_UNIT_SPHERE

SUBROUTINE MAKE_UNIT_SEMI_SPHERE(S)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    REAL(F32)::X,Y,Z,R,X0,Y0,Z0
    INTEGER::J

    X0=0.0;Y0=0.0;Z0=0.0;
    DO J=1,S%NV
        X0=X0+S%VVEC(1,J)
        Y0=Y0+S%VVEC(2,J)
        Z0=Z0+S%VVEC(3,J)
    ENDDO
    X0=X0/S%NV;Y0=Y0/S%NV;Z0=Z0/S%NV !CENTER
    DO J=1,S%NV
        X=S%VVEC(1,J)-X0;Y=S%VVEC(2,J)-Y0;Z=S%VVEC(3,J)-Z0
        R=SQRT(X**2+Y**2+Z**2)
        S%VVEC(1,J)=X/R+X0
        S%VVEC(2,J)=Y/R+Y0
        IF(Z>=ZERO_F32) THEN
            S%VVEC(3,J)=Z/R+Z0
        ELSE
            S%VVEC(3,J)=Z0
        ENDIF
    ENDDO
ENDSUBROUTINE MAKE_UNIT_SEMI_SPHERE

SUBROUTINE MAKE_FAN(S,NA,FLAG)
    TYPE(TSURF_TYPE),INTENT(OUT):: S
    INTEGER::NA,FLAG
    INTEGER::N,T,J
    REAL(F32)::PHI,DPHI
    FLAG=0
    IF(NA<5) THEN !MAKES NO SENSE
        FLAG=1
        RETURN
    ENDIF
    N=NA+1;T=NA
    CALL MAKE_TSURF(S,N,T,FLAG)
    IF(FLAG/=0) RETURN
    S%VVEC=0
    S%NVEC=0
    S%NVEC(3,:)=ONE_F32
    DPHI=TWO_F32*PI_F32/NA
    DO J=1,NA
        PHI=DPHI*(J-1)
        S%VVEC(1,J)=COS(PHI)
        S%VVEC(2,J)=SIN(PHI)
        !S%VVEC(3,J) REMAINS ZERO
    ENDDO
    !S%VVEC(:,NA+1) REMAINS ZERO
    DO J=1,NA-1
        S%TVEC(1,J)=NA+1;S%TVEC(2,J)=J;S%TVEC(3,J)=J+1
    ENDDO
    S%TVEC(1,NA)=NA+1;S%TVEC(2,J)=NA;S%TVEC(3,J)=1
ENDSUBROUTINE MAKE_FAN
END MODULE TSURF
