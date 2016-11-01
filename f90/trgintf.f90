!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Shimon Panfil                                           !
! Copyright (c) Shimon Panfil Industrial Physics and Simulations  !
! http://industrialphys.com                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE TRGINT_MODULE
    IMPLICIT NONE
    PRIVATE
    REAL(8),PARAMETER :: ZZERO=1.0e-12_8
    REAL(8),PARAMETER :: EPS=1.0e-12_8
    INTEGER,PUBLIC,PARAMETER :: CHARGE_TRG_LENGTH=12
    PUBLIC TRGINT
    PUBLIC POTCTV
CONTAINS
FUNCTION NEWG(U,V,Q,ZZ,A,B1,C) RESULT(RES)
    REAL(8),INTENT(IN) :: U,V,Q,ZZ,A,B1,C
    REAL(8) :: VV,Z,VZ,SVQ,S,T,WW,QQ,VQ,SS,TT,BQV
    REAL(8) :: G1,G2,G3,G4,G5,L1,L2,L3,G6,B
    REAL(8) :: VV2,VV1,BVVQQ,TS,S2T2,CVTS,RES
    VV=V*V
    VV2=(VV+1.0_8)
    VV1=VV2-2.0_8
    B=B1/VV2
    Z=SQRT(ZZ)
    QQ=Q*Q
    WW=ZZ+QQ
    S=SQRT(U*U+WW)-U
    T=S+2.0_8*U
    VZ=V*Z
    VQ=V*Q
    SVQ=S+VQ
    G1=2.0_8*C*Z*VV2*ATAN2(SVQ,VZ)
    SS=S*S
    TT=T*T
    TS=T-S
    CVTS=C*V*TS
    S2T2=0.5_8*(SS+TT)*VV
    BQV=B*Q*V
    G2=0.25_8*A*V*(SS-TT)+0.5_8*B*S2T2-BQV*VV*T-BQV*S+CVTS
    L1=LOG(SVQ*SVQ+VZ*VZ)
    L2=LOG(S)
    L3=LOG(VV2)
    BVVQQ=2.0_8*B*VV*QQ
    G3=L2*(A*V*WW+B*WW+2.0_8*Q*C)
    G4=(L3+L2-L1)*(B*S2T2-BQV*VV1*TS+CVTS)
    G5=L1*(-0.5_8*(VV*VV+1.0_8)*WW*B+C*VV1*Q)
    G6=(L1-L2)*BVVQQ
    RES=G1+G2+G3+G4+G5+G6
END FUNCTION NEWG
FUNCTION TRGINT(x,y,z,x0IN,y0IN,z0IN,q0,x1,y1,z1,q1,x2,y2,z2,q2) RESULT(RES)
    REAL(8),INTENT(IN) :: x,y,z,x0IN,y0IN,z0IN,q0,x1,y1,z1,q1,x2,y2,z2,q2
    REAL(8) :: x10,y10,z10,x20,y20,z20,X0,Y0,Z0
    REAL(8) :: nx,ny,nz,kx,ky,kz,n2x,n2y,n2z
    REAL(8) :: d1,d2,c2,s2,u2,u1,v2,zz,u3,v3,v4,a,b,c
    REAL(8) :: PT1,pt,pb,qt,qb,PB1,ppt,ppb,wtt,wbb,wtb,wbt,Vb,Vt,res
    x10=x1-x0IN
    y10=y1-y0IN
    z10=z1-z0IN
    d1=sqrt(x10*x10+y10*y10+z10*z10)
    if(d1<EPS) THEN
        RES=0.0_8
        GOTO 666
    END IF
    x20=x2-x0IN
    y20=y2-y0IN
    z20=z2-z0IN
    d2=sqrt(x20*x20+y20*y20+z20*z20)
    if(d2<EPS) THEN
        RES=0.0_8
        GOTO 666
    END IF
    nx=x10/d1
    ny=y10/d1
    nz=z10/d1
    n2x=x20/d2
    n2y=y20/d2
    n2z=z20/d2
    c2=nx*n2x+ny*n2y+nz*n2z
    s2=sqrt(1.0_8-c2*c2)
    u2=d2*c2
    u1=d1
    v2=d2*s2
    kx=(n2x-c2*nx)/s2
    ky=(n2y-c2*ny)/s2
    kz=(n2z-c2*nz)/s2
    x0=X0IN-x
    y0=Y0IN-y
    z0=Z0IN-z
    u3=nx*x0+ny*y0+nz*z0
    v3=kx*x0+ky*y0+kz*z0
    zz=x0*x0+y0*y0+z0*z0-u3*u3-v3*v3
    if(zz<ZZERO) zz=ZZERO
    v4=v2+v3
    a=(q1-q0)/d1
    b=(d1*(q2-q0)-d2*c2*(q1-q0))/(d1*v2)
    c=q0-a*u3-b*v3
    PT1=(u2-u1)/v2
    PB1=u2/v2
    ppt=sqrt(1.0_8+PT1*PT1)
    ppb=sqrt(1.0_8+PB1*PB1)
    qt=(u3+u1-PT1*v3)/ppt
    qb=(u3-PB1*v3)/ppb
    wtt=ppt*v4+PT1*qt
    wtb=ppb*v4+PB1*qb
    wbt=ppt*v3+PT1*qt
    wbb=ppb*v3+PB1*qb
    pt=PT1/ppt
    pb=PB1/ppb
    Vt=sqrt((1.0_8+pt)/(1.0_8-pt))
    Vb=sqrt((1.0_8+pb)/(1.0_8-pb))  
    res=ff(wtt,wbt,wtb,wbb,qt,qb,Vt,Vb,zz,a,b,c)
666 CONTINUE
contains
FUNCTION FF(WTT,WBT, WTB, WBB,QT, QB,VT, VB,ZZ, A, B,C) RESULT(RES)
    REAL(8),INTENT(IN) :: WTT,WBT,WTB,WBB,QT,QB,VT,VB,ZZ,A,B,C
    REAL(8) :: RES
    RES=(-NEWG(WTT,VT,QT,ZZ,A,B,C)+NEWG(WBT,VT,QT,ZZ,A,B,C))/(VT*VT+1.0_8) &
    +(NEWG(WTB,VB,QB,ZZ,A,B,C)-NEWG(WBB,VB,QB,ZZ,A,B,C))/(VB*VB+1.0_8)
END FUNCTION FF
END FUNCTION TRGINT
FUNCTION POTCTV(X,Y,Z,CTV,N) RESULT(RES)
    REAL(8),INTENT(IN) :: X,Y,Z
    REAL(8),DIMENSION(:,:),INTENT(IN)::CTV
    INTEGER,INTENT(IN) ::N
    REAL(8)::RES
    INTEGER :: J
    RES=0.0_8
    DO J=1,N
     RES=RES+TRGINT(X,Y,Z,CTV(1,J),CTV(2,J),CTV(3,J),CTV(4,J),&
         CTV(5,J),CTV(6,J),CTV(7,J),CTV(8,J),CTV(9,J),CTV(10,J),CTV(11,J),CTV(12,J))
    END DO
END FUNCTION POTCTV
END MODULE TRGINT_MODULE
