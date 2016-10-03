!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE VEC3D
USE FDEFS
IMPLICIT NONE
PUBLIC
CONTAINS

! FLOAT

REAL(F32) FUNCTION LENGTH_F32(A) 
    REAL(F32)::A(3)
    LENGTH_F32=SQRT(A(1)**2+A(2)**2+A(3)**2)
ENDFUNCTION LENGTH_F32

REAL(F32) FUNCTION DOT_F32(A,B) 
    REAL(F32)::A(3),B(3)
    DOT_F32=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
ENDFUNCTION DOT_F32

COMPLEX(F32) FUNCTION CDOT_F32(A,B)
    COMPLEX(F32)::A(3)
    REAL(F32)::B(3)
    CDOT_F32=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
ENDFUNCTION CDOT_F32

FUNCTION CROSS_F32(A,B) RESULT(C)
    REAL(F32)::A(3),B(3),C(3)
    C(1)=A(2)*B(3)-A(3)*B(2)
    C(2)=A(3)*B(1)-A(1)*B(3)
    C(3)=A(1)*B(2)-A(2)*B(1)
ENDFUNCTION CROSS_F32

FUNCTION CCROSS_F32(A,B) RESULT(C)
    COMPLEX(F32)::A(3),C(3)
    REAL(F32)::B(3)
    C(1)=A(2)*B(3)-A(3)*B(2)
    C(2)=A(3)*B(1)-A(1)*B(3)
    C(3)=A(1)*B(2)-A(2)*B(1)
ENDFUNCTION CCROSS_F32

FUNCTION MXV_F32(A,B) RESULT(C)
    REAL(F32)::A(3,3),C(3),B(3)
        C(1)=A(1,1)*B(1)+A(1,2)*B(2)+A(1,3)*B(3)
        C(2)=A(2,1)*B(1)+A(2,2)*B(2)+A(2,3)*B(3)
        C(3)=A(3,1)*B(1)+A(3,2)*B(2)+A(3,3)*B(3)
ENDFUNCTION MXV_F32

FUNCTION CMXV_F32(A,B) RESULT(C)
    COMPLEX(F32)::A(3,3),C(3)
    REAL(F32)::B(3)    
        C(1)=A(1,1)*B(1)+A(1,2)*B(2)+A(1,3)*B(3)
        C(2)=A(2,1)*B(1)+A(2,2)*B(2)+A(2,3)*B(3)
        C(3)=A(3,1)*B(1)+A(3,2)*B(2)+A(3,3)*B(3)
ENDFUNCTION CMXV_F32

!making rotation matrix RxRyRz
FUNCTION ROT_MAT_F32 (AX,AY,AZ) RESULT(M)
    REAL(F32)::AX,AY,AZ,M(3,3)
    REAL(F32)::CZ,SZ,CY,SY,CX,SX
    CX=COS(AX*TORAD_F32);SX=SIN(AX*TORAD_F32)
    CY=COS(AY*TORAD_F32);SY=SIN(AY*TORAD_F32)
    CZ=COS(AZ*TORAD_F32);SZ=SIN(AZ*TORAD_F32)
    M(1,1)=CY*CZ;          M(1,2)=-CY*SZ;         M(1,3)=SY
    M(2,1)=SX*SY*CZ+CX*SZ; M(2,2)=-SX*SY*SZ+CX*CZ;M(2,3)=-SX*CY
    M(3,1)=-CX*SY*CZ+SX*SZ;M(3,2)=CX*SY*SZ+SX*CZ; M(3,3)=CX*CY
ENDFUNCTION ROT_MAT_F32

! GET PATCH PARAMETERS
FUNCTION GET_TRGV(TRG) RESULT (RES)
    REAL(F32)::TRG(3,3),RES(3)
    REAL(F32):: R(3),U(3),V(3),A,P,H
    R=TRG(:,1)
    U=TRG(:,2)-R
    V=TRG(:,3)-R
    A=LENGTH_F32(U)
    P=DOT_F32(V,U)/A
    H=SQRT(DOT_F32(V,V)-P**2)
    RES(1)=A;RES(2)=P;RES(3)=H
ENDFUNCTION GET_TRGV 

FUNCTION TRG_AREA(TRG) RESULT (RES)
    REAL(F32)::TRG(3,3),RES
    REAL(F32):: W(3),U(3),V(3)
    U=TRG(:,2)-TRG(:,1)
    V=TRG(:,3)-TRG(:,1)
    W=CROSS_F32(U,V)
    RES=HALF_F32*LENGTH_F32(W)
ENDFUNCTION  TRG_AREA

! DOUBLE

REAL(F64) FUNCTION LENGTH_F64(A) 
    REAL(F64)::A(3)
    LENGTH_F64=SQRT(A(1)**2+A(2)**2+A(3)**2)
ENDFUNCTION LENGTH_F64

REAL(F64) FUNCTION DOT_F64(A,B) 
    REAL(F64)::A(3),B(3)
    DOT_F64=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
ENDFUNCTION DOT_F64

COMPLEX(F64) FUNCTION CDOT_F64(A,B)
    COMPLEX(F64)::A(3)
    REAL(F64)::B(3)
    CDOT_F64=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
ENDFUNCTION CDOT_F64

FUNCTION CROSS_F64(A,B) RESULT(C)
    REAL(F64)::A(3),B(3),C(3)
    C(1)=A(2)*B(3)-A(3)*B(2)
    C(2)=A(3)*B(1)-A(1)*B(3)
    C(3)=A(1)*B(2)-A(2)*B(1)
ENDFUNCTION CROSS_F64

FUNCTION CCROSS_F64(A,B) RESULT(C)
    COMPLEX(F64)::A(3),C(3)
    REAL(F64)::B(3)
    C(1)=A(2)*B(3)-A(3)*B(2)
    C(2)=A(3)*B(1)-A(1)*B(3)
    C(3)=A(1)*B(2)-A(2)*B(1)
ENDFUNCTION CCROSS_F64

FUNCTION MXV_F64(A,B) RESULT(C)
    REAL(F64)::A(3,3),C(3),B(3)
        C(1)=A(1,1)*B(1)+A(1,2)*B(2)+A(1,3)*B(3)
        C(2)=A(2,1)*B(1)+A(2,2)*B(2)+A(2,3)*B(3)
        C(3)=A(3,1)*B(1)+A(3,2)*B(2)+A(3,3)*B(3)
ENDFUNCTION MXV_F64

FUNCTION CMXV_F64(A,B) RESULT(C)
    COMPLEX(F64)::A(3,3),C(3)
    REAL(F64)::B(3)    
        C(1)=A(1,1)*B(1)+A(1,2)*B(2)+A(1,3)*B(3)
        C(2)=A(2,1)*B(1)+A(2,2)*B(2)+A(2,3)*B(3)
        C(3)=A(3,1)*B(1)+A(3,2)*B(2)+A(3,3)*B(3)
ENDFUNCTION CMXV_F64

!making rotation matrix RxRyRz
FUNCTION ROT_MAT_F64 (AX,AY,AZ) RESULT(M)
    REAL(F64)::AX,AY,AZ,M(3,3)
    REAL(F64)::CZ,SZ,CY,SY,CX,SX
    CX=COS(AX*TORAD_F64);SX=SIN(AX*TORAD_F64)
    CY=COS(AY*TORAD_F64);SY=SIN(AY*TORAD_F64)
    CZ=COS(AZ*TORAD_F64);SZ=SIN(AZ*TORAD_F64)
    M(1,1)=CY*CZ;          M(1,2)=-CY*SZ;         M(1,3)=SY
    M(2,1)=SX*SY*CZ+CX*SZ; M(2,2)=-SX*SY*SZ+CX*CZ;M(2,3)=-SX*CY
    M(3,1)=-CX*SY*CZ+SX*SZ;M(3,2)=CX*SY*SZ+SX*CZ; M(3,3)=CX*CY
ENDFUNCTION ROT_MAT_F64

ENDMODULE VEC3D
