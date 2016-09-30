!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations        !
! http://industrialphys.com                                !
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE TSURF
USE FDEFS
USE VEC3D
IMPLICIT NONE
PUBLIC

TYPE :: TSURF_TYPE ! TRIANGULATED SURFACE
    INTEGER,ALLOCATABLE,DIMENSION(:,:)::TVEC! (3,NT) TRIANGLES
    REAL(F32),ALLOCATABLE,DIMENSION(:,:)::VVEC !(3,NV) VERTICES
    INTEGER::NT ! NUMBER OF TRIANGLES
    INTEGER::NV ! NUMBER OF VERTICES
    REAL(F32),DIMENSION(3,3)::RPR !orientation
    REAL(F32),DIMENSION(3)::ORG !origin
END TYPE TSURF_TYPE  

REAL(F32),PARAMETER :: MIN_NODE_DIST=1.0E-8
REAL(F32),PARAMETER :: AREA_ZERO=1.0E-12
INTEGER,PARAMETER::MAX_CON=8

CONTAINS

SUBROUTINE CLEAN_TSURF(S)
    TYPE(TSURF_TYPE),INTENT(INOUT) :: S
    IF(ALLOCATED(S%VVEC)) DEALLOCATE(S%VVEC)
    IF(ALLOCATED(S%TVEC)) DEALLOCATE(S%TVEC)
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
    S%ORG=ZERO_F32
    S%RPR=ZERO_F32
    S%RPR(1,1)=ONE_F32;S%RPR(2,2)=ONE_F32;S%RPR(3,3)=ONE_F32
    S%NV=N;S%NT=T
    ALLOCATE(S%VVEC(3,S%NV),S%TVEC(3,S%NT),STAT=ERR)
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
    WRITE(WRUNITA,*) 'VERTICES'
    DO I=1,S%NV
        WRITE(WRUNITA,*) S%VVEC(1,I),S%VVEC(2,I),S%VVEC(3,I)
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
    CLOSE(RDUNITB)    
    FLAG=0
END SUBROUTINE READ_TSURF_B

! REGULAR POLYHEDRA INSCRIBED IN UNIT SPHERE
! CONSTRUCTED WITH MAPLE
SUBROUTINE TETRAHEDRON(S)
    TYPE(TSURF_TYPE),INTENT(OUT):: S
    INTEGER:: FLAG,N,T
    REAL(F32)::S3
    S3=1.0/SQRT(3.0)
    N=4
    T=4
    CALL MAKE_TSURF(S,N,T,FLAG)
    S%VVEC(1,1)=S3;S%VVEC(2,1)=S3;S%VVEC(3,1)=S3
    S%VVEC(1,2)=S3;S%VVEC(2,2)=-S3;S%VVEC(3,2)=-S3
    S%VVEC(1,3)=-S3;S%VVEC(2,3)=S3;S%VVEC(3,3)=-S3
    S%VVEC(1,4)=-S3;S%VVEC(2,4)=-S3;S%VVEC(3,4)=S3
    
    S%TVEC(1,1)=1;S%TVEC(2,1)=2;S%TVEC(3,1)=3
    S%TVEC(1,2)=1;S%TVEC(2,2)=4;S%TVEC(3,2)=2
    S%TVEC(1,3)=1;S%TVEC(2,3)=3;S%TVEC(3,3)=4
    S%TVEC(1,4)=2;S%TVEC(2,4)=4;S%TVEC(3,4)=3
ENDSUBROUTINE TETRAHEDRON

SUBROUTINE HEXAHEDRON(S)
    TYPE(TSURF_TYPE),INTENT(OUT):: S
    INTEGER:: FLAG,N,T
    REAL(F32)::S3
    INTEGER::P(4,6),J
    S3=1.0/SQRT(3.0)
    N=14
    T=24
    CALL MAKE_TSURF(S,N,T,FLAG)
    S%VVEC(1,1)=S3;S%VVEC(2,1)=S3;S%VVEC(3,1)=S3
    S%VVEC(1,2)=S3;S%VVEC(2,2)=S3;S%VVEC(3,2)=-S3
    S%VVEC(1,3)=S3;S%VVEC(2,3)=-S3;S%VVEC(3,3)=S3
    S%VVEC(1,4)=S3;S%VVEC(2,4)=-S3;S%VVEC(3,4)=-S3
    S%VVEC(1,5)=-S3;S%VVEC(2,5)=S3;S%VVEC(3,5)=S3
    S%VVEC(1,6)=-S3;S%VVEC(2,6)=S3;S%VVEC(3,6)=-S3
    S%VVEC(1,7)=-S3;S%VVEC(2,7)=-S3;S%VVEC(3,7)=S3
    S%VVEC(1,8)=-S3;S%VVEC(2,8)=-S3;S%VVEC(3,8)=-S3
    P(1,1)=1;P(2,1)=3;P(3,1)=4;P(4,1)=2
    P(1,2)=5;P(2,2)=6;P(3,2)=8;P(4,2)=7
    P(1,3)=1;P(2,3)=2;P(3,3)=6;P(4,3)=5
    P(1,4)=3;P(2,4)=7;P(3,4)=8;P(4,4)=4
    P(1,5)=1;P(2,5)=5;P(3,5)=7;P(4,5)=3
    P(1,6)=4;P(2,6)=8;P(3,6)=6;P(4,6)=2
    DO J=1,6
    S%VVEC(:,8+J)=0.25*(S%VVEC(:,P(1,J))+S%VVEC(:,P(2,J)) &
                 &   +S%VVEC(:,P(3,J))+S%VVEC(:,P(4,J)))
    ENDDO
    DO J=1,6
    S%TVEC(1,4*(J-1)+1)=8+J;S%TVEC(2,4*(J-1)+1)=P(1,J);S%TVEC(3,4*(J-1)+1)=P(2,J)
    S%TVEC(1,4*(J-1)+2)=8+J;S%TVEC(2,4*(J-1)+2)=P(2,J);S%TVEC(3,4*(J-1)+2)=P(3,J)
    S%TVEC(1,4*(J-1)+3)=8+J;S%TVEC(2,4*(J-1)+3)=P(3,J);S%TVEC(3,4*(J-1)+3)=P(4,J)
    S%TVEC(1,4*(J-1)+4)=8+J;S%TVEC(2,4*(J-1)+4)=P(4,J);S%TVEC(3,4*(J-1)+4)=P(1,J)
    ENDDO
ENDSUBROUTINE HEXAHEDRON

SUBROUTINE OCTAHEDRON(S)
    TYPE(TSURF_TYPE),INTENT(OUT):: S
    INTEGER:: FLAG,N,T
    N=6
    T=8
    CALL MAKE_TSURF(S,N,T,FLAG)
    S%VVEC(1,1)=1.0;S%VVEC(2,1)=0.0;S%VVEC(3,1)=0.0
    S%VVEC(1,2)=-1.0;S%VVEC(2,2)=0.0;S%VVEC(3,2)=0.0
    S%VVEC(1,3)=0.0;S%VVEC(2,3)=1.0;S%VVEC(3,3)=0.0
    S%VVEC(1,4)=0.0;S%VVEC(2,4)=-1.0;S%VVEC(3,4)=0.0
    S%VVEC(1,5)=0.0;S%VVEC(2,5)=0.0;S%VVEC(3,5)=1.0
    S%VVEC(1,6)=0.0;S%VVEC(2,6)=0.0;S%VVEC(3,6)=-1.0
    
    S%TVEC(1,1)=5;S%TVEC(2,1)=1;S%TVEC(3,1)=3
    S%TVEC(1,2)=6;S%TVEC(2,2)=3;S%TVEC(3,2)=1
    S%TVEC(1,3)=5;S%TVEC(2,3)=4;S%TVEC(3,3)=1
    S%TVEC(1,4)=6;S%TVEC(2,4)=1;S%TVEC(3,4)=4
    S%TVEC(1,5)=5;S%TVEC(2,5)=3;S%TVEC(3,5)=2
    S%TVEC(1,6)=6;S%TVEC(2,6)=2;S%TVEC(3,6)=3
    S%TVEC(1,7)=5;S%TVEC(2,7)=2;S%TVEC(3,7)=4
    S%TVEC(1,8)=6;S%TVEC(2,8)=4;S%TVEC(3,8)=2
ENDSUBROUTINE OCTAHEDRON

SUBROUTINE DODECAHEDRON(S)
    TYPE(TSURF_TYPE),INTENT(OUT):: S
    INTEGER:: FLAG,N,T
    REAL(F32)::x,y,z
    INTEGER::P(5,12),J
    x=sqrt(5.0)-1.0
    y=x+2.0
    z=2*sqrt(3.0)
    N=32
    T=60
    CALL MAKE_TSURF(S,N,T,FLAG)
    S%VVEC(1,1)=x;  S%VVEC(2,1)=-y; S%VVEC(3,1)=0
    S%VVEC(1,2)=-x; S%VVEC(2,2)=-y; S%VVEC(3,2)=0
    S%VVEC(1,3)=-2; S%VVEC(2,3)=-2; S%VVEC(3,3)=2
    S%VVEC(1,4)=0;  S%VVEC(2,4)=-x; S%VVEC(3,4)=y
    S%VVEC(1,5)=2;  S%VVEC(2,5)=-2; S%VVEC(3,5)=2
    S%VVEC(1,6)=0;  S%VVEC(2,6)=x;  S%VVEC(3,6)=-y
    S%VVEC(1,7)=-2; S%VVEC(2,7)=2;  S%VVEC(3,7)=-2
    S%VVEC(1,8)=-x; S%VVEC(2,8)=y;  S%VVEC(3,8)=0
    S%VVEC(1,9)=x;  S%VVEC(2,9)=y;   S%VVEC(3,9)=0
    S%VVEC(1,10)=2;  S%VVEC(2,10)=2; S%VVEC(3,10)=-2
    S%VVEC(1,11)=0;  S%VVEC(2,11)=-x; S%VVEC(3,11)=-y
    S%VVEC(1,12)=-2; S%VVEC(2,12)=-2; S%VVEC(3,12)=-2
    S%VVEC(1,13)=-y; S%VVEC(2,13)=0; S%VVEC(3,13)=-x
    S%VVEC(1,14)=-y;  S%VVEC(2,14)=0; S%VVEC(3,14)=x
    S%VVEC(1,15)=-2;  S%VVEC(2,15)=2; S%VVEC(3,15)=2
    S%VVEC(1,16)=0;  S%VVEC(2,16)=x;  S%VVEC(3,16)=y
    S%VVEC(1,17)=2; S%VVEC(2,17)=2;  S%VVEC(3,17)=2
    S%VVEC(1,18)=y; S%VVEC(2,18)=0;  S%VVEC(3,18)=x
    S%VVEC(1,19)=y;  S%VVEC(2,19)=0;   S%VVEC(3,19)=-x
    S%VVEC(1,20)=2;  S%VVEC(2,20)=-2; S%VVEC(3,20)=-2

!    P(1,1)=8;   P(2,1)=9;   P(3,1)=17;  P(4,1)=16;  P(5,1)=15
    P(1,1)=8;   P(2,1)=15;   P(3,1)=16;  P(4,1)=17;  P(5,1)=9
    P(1,2)=6;   P(2,2)=7;   P(3,2)=8;   P(4,2)=9;   P(5,2)=10
    P(1,3)=1;   P(2,3)=5;   P(3,3)=4;   P(4,3)=3;   P(5,3)=2
    P(1,4)=1;   P(2,4)=2;   P(3,4)=12;  P(4,4)=11;  P(5,4)=20
    P(1,5)=4;   P(2,5)=5;   P(3,5)=18;  P(4,5)=17;  P(5,5)=16
    P(1,6)=6;   P(2,6)=10;  P(3,6)=19;  P(4,6)=20;  P(5,6)=11
    P(1,7)=3;   P(2,7)=4;   P(3,7)=16;  P(4,7)=15;  P(5,7)=14
    P(1,8)=6;   P(2,8)=11;  P(3,8)=12;  P(4,8)=13;  P(5,8)=7
    P(1,9)=9;   P(2,9)=17;  P(3,9)=18;  P(4,9)=19;  P(5,9)=10
    P(1,10)=1;  P(2,10)=20; P(3,10)=19; P(4,10)=18; P(5,10)=5
    P(1,11)=7;  P(2,11)=13; P(3,11)=14; P(4,11)=15; P(5,11)=8
    P(1,12)=2;  P(2,12)=3;  P(3,12)=14; P(4,12)=13; P(5,12)=12

    DO J=1,12
    S%VVEC(:,20+J)=0.2*(S%VVEC(:,P(1,J))+S%VVEC(:,P(2,J)) &
        &   +S%VVEC(:,P(3,J))+S%VVEC(:,P(4,J))+S%VVEC(:,P(5,J)))
    ENDDO
    S%VVEC=S%VVEC/z
    DO J=1,12
    S%TVEC(1,5*(J-1)+1)=20+J;S%TVEC(2,5*(J-1)+1)=P(1,J);S%TVEC(3,5*(J-1)+1)=P(2,J)
    S%TVEC(1,5*(J-1)+2)=20+J;S%TVEC(2,5*(J-1)+2)=P(2,J);S%TVEC(3,5*(J-1)+2)=P(3,J)
    S%TVEC(1,5*(J-1)+3)=20+J;S%TVEC(2,5*(J-1)+3)=P(3,J);S%TVEC(3,5*(J-1)+3)=P(4,J)
    S%TVEC(1,5*(J-1)+4)=20+J;S%TVEC(2,5*(J-1)+4)=P(4,J);S%TVEC(3,5*(J-1)+4)=P(5,J)
    S%TVEC(1,5*(J-1)+5)=20+J;S%TVEC(2,5*(J-1)+5)=P(5,J);S%TVEC(3,5*(J-1)+5)=P(1,J)
    ENDDO
ENDSUBROUTINE DODECAHEDRON

SUBROUTINE ICOSAHEDRON(S)
    TYPE(TSURF_TYPE),INTENT(OUT):: S
    INTEGER:: FLAG,N,T
    REAL(F32)::x,z
    x=(sqrt(5.0)+1.0)*0.5
    z=sqrt(0.5*sqrt(5.0)+2.5)
    N=12
    T=20
    CALL MAKE_TSURF(S,N,T,FLAG)
    S%VVEC(1,1)=0;  S%VVEC(2,1)=x;  S%VVEC(3,1)=1
    S%VVEC(1,2)=0;  S%VVEC(2,2)=x;  S%VVEC(3,2)=-1
    S%VVEC(1,3)=0;  S%VVEC(2,3)=-x; S%VVEC(3,3)=1
    S%VVEC(1,4)=0;  S%VVEC(2,4)=-x; S%VVEC(3,4)=-1
    S%VVEC(1,5)=1;  S%VVEC(2,5)=0;  S%VVEC(3,5)=x
    S%VVEC(1,6)=1;  S%VVEC(2,6)=0;  S%VVEC(3,6)=-x
    S%VVEC(1,7)=-1; S%VVEC(2,7)=0;  S%VVEC(3,7)=x
    S%VVEC(1,8)=-1; S%VVEC(2,8)=0;  S%VVEC(3,8)=-x
    S%VVEC(1,9)=x;  S%VVEC(2,9)=1;  S%VVEC(3,9)=0
    S%VVEC(1,10)=x; S%VVEC(2,10)=-1;S%VVEC(3,10)=0
    S%VVEC(1,11)=-x;S%VVEC(2,11)=1; S%VVEC(3,11)=0
    S%VVEC(1,12)=-x;S%VVEC(2,12)=-1;S%VVEC(3,12)=0
    S%VVEC=S%VVEC/z
    S%TVEC(1,1)=1;  S%TVEC(2,1)=5;      S%TVEC(3,1)=9
    S%TVEC(1,2)=1;  S%TVEC(2,2)=9;      S%TVEC(3,2)=2
    S%TVEC(1,3)=1;  S%TVEC(2,3)=2;      S%TVEC(3,3)=11
    S%TVEC(1,4)=1;  S%TVEC(2,4)=11;     S%TVEC(3,4)=7
    S%TVEC(1,5)=1;  S%TVEC(2,5)=7;      S%TVEC(3,5)=5
    S%TVEC(1,6)=2;  S%TVEC(2,6)=9;      S%TVEC(3,6)=6
    S%TVEC(1,7)=2;  S%TVEC(2,7)=6;      S%TVEC(3,7)=8
    S%TVEC(1,8)=2;  S%TVEC(2,8)=8;      S%TVEC(3,8)=11
    S%TVEC(1,9)=3;  S%TVEC(2,9)=4;      S%TVEC(3,9)=10
    S%TVEC(1,10)=3; S%TVEC(2,10)=10;    S%TVEC(3,10)=5
    S%TVEC(1,11)=3; S%TVEC(2,11)=5;     S%TVEC(3,11)=7
    S%TVEC(1,12)=3; S%TVEC(2,12)=7;     S%TVEC(3,12)=12
    S%TVEC(1,13)=3; S%TVEC(2,13)=12;    S%TVEC(3,13)=4
    S%TVEC(1,14)=4; S%TVEC(2,14)=12;    S%TVEC(3,14)=8
    S%TVEC(1,15)=4; S%TVEC(2,15)=8;     S%TVEC(3,15)=6
    S%TVEC(1,16)=4; S%TVEC(2,16)=6;     S%TVEC(3,16)=10
    S%TVEC(1,17)=5; S%TVEC(2,17)=10;    S%TVEC(3,17)=9
    S%TVEC(1,18)=6; S%TVEC(2,18)=9;     S%TVEC(3,18)=10
    S%TVEC(1,19)=7; S%TVEC(2,19)=11;    S%TVEC(3,19)=12
    S%TVEC(1,20)=8; S%TVEC(2,20)=12;    S%TVEC(3,20)=11
ENDSUBROUTINE ICOSAHEDRON

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

SUBROUTINE DELETE_REDUNDANT(S,ND_START,TR_START,FLAG)
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    INTEGER::ND_START,TR_START,FLAG
    INTEGER,ALLOCATABLE,DIMENSION(:,:)::TVEC_TMP
    REAL(F32),ALLOCATABLE,DIMENSION(:,:)::VVEC_TMP
    INTEGER::ND_NEW,NT_NEW
    INTEGER::J,K,I,M,N
    FLAG=0
    ND_NEW=S%NV
    NT_NEW=S%NT
    DO J=ND_START,S%NV
        K=J+1
        DO WHILE(K<=ND_NEW)
            IF(DIST(J,K)>MIN_NODE_DIST) THEN 
                K=K+1
            ELSE ! COINCIDING NODES
                DO I=TR_START,S%NT !UPDATE TRIANGLES
                    DO M=1,3
                        N=S%TVEC(M,I)
                        IF(N==K) S%TVEC(M,I)=J
                        IF(N>K) S%TVEC(M,I)=S%TVEC(M,I)-1
                    ENDDO
                ENDDO
                ND_NEW=ND_NEW-1
                DO M=K+1,S%NV !UPDATE NODES
                    S%VVEC(:,M-1)=S%VVEC(:,M)
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    DO J=TR_START,S%NT
        K=J+1
        DO WHILE(K<=S%NT)
            IF(.NOT.TRG_EQ(J,K)) THEN 
                K=K+1
            ELSE ! COINCIDING TRIANGLES
                NT_NEW=NT_NEW-1
                DO M=K+1,S%NV 
                    S%TVEC(:,M-1)=S%TVEC(:,M)
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    IF((NT_NEW.EQ.S%NT).AND.(ND_NEW.EQ.S%NV)) RETURN
    ALLOCATE(VVEC_TMP(3,ND_NEW),TVEC_TMP(3,NT_NEW),STAT=FLAG)
    IF(FLAG/=0) RETURN
    VVEC_TMP(:,:)=S%VVEC(:,1:ND_NEW)
    TVEC_TMP(:,:)=S%TVEC(:,1:NT_NEW)
    CALL MAKE_TSURF(S,ND_NEW,NT_NEW,FLAG)
    IF(FLAG/=0) RETURN
    S%VVEC=VVEC_TMP
    S%TVEC=TVEC_TMP
    DEALLOCATE(VVEC_TMP,TVEC_TMP)
CONTAINS
    REAL(F32) FUNCTION DIST(N1,N2) ! DISTANCE BETWEEN NODES
    INTEGER::N1,N2
    REAL(F32):: X,Y,Z
    X=S%VVEC(1,N1)-S%VVEC(1,N2)
    Y=S%VVEC(2,N1)-S%VVEC(2,N2)
    Z=S%VVEC(3,N1)-S%VVEC(3,N2)
    DIST=SQRT(X**2+Y**2+Z**2)
    ENDFUNCTION DIST
    LOGICAL FUNCTION TRG_EQ(T1,T2) ! COINSIDING TRIANGLES
    INTEGER::T1,T2
    INTEGER::N1,N2,N3,M1,M2,M3
    N1=S%TVEC(1,T1);N2=S%TVEC(2,T1);N3=S%TVEC(3,T1)
    M1=S%TVEC(1,T2);M2=S%TVEC(2,T2);M3=S%TVEC(3,T2)
    TRG_EQ=((N1.EQ.M1).OR.(N1.EQ.M2).OR.(N1.EQ.M3)) &
            & .AND.((N2.EQ.M1).OR.(N2.EQ.M2).OR.(N2.EQ.M3))  &
            & .AND.((N3.EQ.M1).OR.(N3.EQ.M2).OR.(N3.EQ.M3))
    ENDFUNCTION TRG_EQ
ENDSUBROUTINE DELETE_REDUNDANT

SUBROUTINE REFINE_TRIANGULATION_2(SOLD,SNEW,FLAG)
    TYPE(TSURF_TYPE),INTENT(IN):: SOLD
    TYPE(TSURF_TYPE),INTENT(OUT):: SNEW
    INTEGER,INTENT(OUT)::FLAG
    INTEGER::TOLD,NOLD,J,N,N1,N2,N3,M1,M2,M3
    FLAG=0
    TOLD=SOLD%NT
    NOLD=SOLD%NV
    CALL MAKE_TSURF(SNEW,NOLD+3*TOLD,4*TOLD,FLAG)
    IF(FLAG/=0) RETURN
    SNEW%VVEC(:,1:NOLD)=SOLD%VVEC
    N=NOLD
    DO J=1,TOLD
        N1=SOLD%TVEC(1,J);N2=SOLD%TVEC(2,J);N3=SOLD%TVEC(3,J)
        M1=N+1;M2=N+2;M3=N+3
        SNEW%VVEC(:,M1)=0.5*(SNEW%VVEC(:,N1)+SNEW%VVEC(:,N2))
        SNEW%VVEC(:,M2)=0.5*(SNEW%VVEC(:,N2)+SNEW%VVEC(:,N3))
        SNEW%VVEC(:,M3)=0.5*(SNEW%VVEC(:,N3)+SNEW%VVEC(:,N1))
        N=N+3
        SNEW%TVEC(1,4*(J-1)+1)=N1
        SNEW%TVEC(2,4*(J-1)+1)=M1
        SNEW%TVEC(3,4*(J-1)+1)=M3
        SNEW%TVEC(1,4*(J-1)+2)=M1
        SNEW%TVEC(2,4*(J-1)+2)=N2
        SNEW%TVEC(3,4*(J-1)+2)=M2
        SNEW%TVEC(1,4*(J-1)+3)=M2
        SNEW%TVEC(2,4*(J-1)+3)=N3
        SNEW%TVEC(3,4*(J-1)+3)=M3
        SNEW%TVEC(1,4*(J-1)+4)=M1
        SNEW%TVEC(2,4*(J-1)+4)=M2
        SNEW%TVEC(3,4*(J-1)+4)=M3
    ENDDO
  !  CALL DELETE_REDUNDANT(SNEW,NOLD,1,FLAG)
ENDSUBROUTINE REFINE_TRIANGULATION_2
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
        if(z>=zero_f32) then
            S%VVEC(3,J)=Z/R+Z0
        else
            S%VVEC(3,J)=Z0
        endif
    ENDDO
ENDSUBROUTINE MAKE_UNIT_SEMI_SPHERE

SUBROUTINE MAKE_REPERS(RVEC,S,FLAG)
    REAL(F32),DIMENSION(3,3,*)::RVEC
    TYPE(TSURF_TYPE),INTENT(INOUT):: S
    INTEGER,INTENT(OUT)::FLAG
    INTEGER::J,N1,N2,N3
    REAL(F32)::R1(3),U(3),V(3),W(3)   
    DO J=1,S%NT
        N1=S%TVEC(1,J)
        N2=S%TVEC(2,J)
        N3=S%TVEC(3,J)
        R1=S%VVEC(:,N1)
        U=S%VVEC(:,N2)-R1
        U=U/LENGTH_F32(U)
        V=S%VVEC(:,N3)-R1 
        V=V-U*DOT_F32(V,U)
        V=V/LENGTH_F32(V)
        W=CROSS_F32(U,V)
        RVEC(:,1,J)=U
        RVEC(:,2,J)=V
        RVEC(:,3,J)=W
    ENDDO
    FLAG=0
END SUBROUTINE MAKE_REPERS
! not shifted
subroutine get_trg(trg,ts,t)
    real(f32)::trg(3,3)
    TYPE(TSURF_TYPE):: tS
    integer::t
    real(f32)::v1(3),v2(3),v3(3)
    v1=ts%vvec(:,ts%tvec(1,t))
    v2=ts%vvec(:,ts%tvec(2,t))
    v3=ts%vvec(:,ts%tvec(3,t))
    trg(:,1)=v1(1)*ts%rpr(:,1)+v1(2)*ts%rpr(:,2)+v1(3)*ts%rpr(:,3)
    trg(:,2)=v2(1)*ts%rpr(:,1)+v2(2)*ts%rpr(:,2)+v2(3)*ts%rpr(:,3)
    trg(:,3)=v3(1)*ts%rpr(:,1)+v3(2)*ts%rpr(:,2)+v3(3)*ts%rpr(:,3)
endsubroutine get_trg
subroutine get_trg_shifted(trg,ts,t)
    real(f32)::trg(3,3)
    TYPE(TSURF_TYPE):: tS
    integer::t
    real(f32)::v1(3),v2(3),v3(3)
    v1=ts%vvec(:,ts%tvec(1,t))
    v2=ts%vvec(:,ts%tvec(2,t))
    v3=ts%vvec(:,ts%tvec(3,t))
    trg(:,1)=v1(1)*ts%rpr(:,1)+v1(2)*ts%rpr(:,2)+v1(3)*ts%rpr(:,3)+ts%org
    trg(:,2)=v2(1)*ts%rpr(:,1)+v2(2)*ts%rpr(:,2)+v2(3)*ts%rpr(:,3)+ts%org
    trg(:,3)=v3(1)*ts%rpr(:,1)+v3(2)*ts%rpr(:,2)+v3(3)*ts%rpr(:,3)+ts%org
endsubroutine get_trg_shifted

END MODULE TSURF