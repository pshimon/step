!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Shimon Panfil                                           !
! Copyright (c) Shimon Panfil Industrial Physics and Simulations  !
! http://industrialphys.com                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE int_trgf
    IMPLICIT NONE
    PRIVATE
    REAL(8),PARAMETER :: ZZERO=1.0e-12_8
    REAL(8),PARAMETER :: EPS=1.0e-12_8
    INTEGER,PUBLIC,PARAMETER :: CHARGE_TRG_LENGTH=12
    PUBLIC int_trg

CONTAINS
FUNCTION int_trg(x,y,z,x0,y0,z0,q0,xx1,yy1,zz1,q1,xx2,yy2,zz2,q2) RESULT(RES)
    REAL(8),INTENT(IN) :: x,y,z,x0,y0,z0,q0,xx1,yy1,zz1,q1,xx2,yy2,zz2,q2
    REAL(8) ::rr11,rr12,rr13,rr22,rr23,rr33,m,mm 
    REAL(8) ::nrm,w1,w2,r1,r2,v1,v2,u1,u2,s1u,s2u,s1d,s2d 
    REAL(8) ::a1,a2,b1,b2,b,c1,c2,c,r3m,r1m,r2m,res 
    REAL(8) ::x1,y1,z1,x2,y2,z2,x3,y3,z3
     x1=xx1-x0
     y1=yy1-y0
     z1=zz1-z0
     x2=xx2-x0
     y2=yy2-y0
     z2=zz2-z0
     x3=x-x0
     y3=y-y0
     z3=z-z0
    rr11=x1**2+y1**2+z1**2
    r1=sqrt(rr11)
    rr12=x1*x2+y1*y2+z1*z2
    rr13=x1*x3+y1*y3+z1*z3
    rr22=x2**2+y2**2+z2**2
    r2=sqrt(rr22)
    rr23=x2*x3+y2*y3+z2*z3
    rr33=x3**2+y3**2+z3**2
    mm=0.25d0*(rr11+rr22-2.0d0*rr12)
    m=sqrt(mm)
    nrm=0.25d0*sqrt(rr11*rr22-rr12*rr12)/mm
    s1d=-rr13/r1
    s2d=-rr23/r2
    s1u=s1d+r1
    s2u=s2d+r2
    w1=sqrt(rr33-rr13*rr13/rr11)
    w2=sqrt(rr33-rr23*rr23/rr22)
    a1=(q1-q2)/r1
    a2=(q1-q2)/r2
    c=(q1+q2-2.0d0*q0)*m-0.25d0*(q1-q2)*(rr11-rr22)/m
    b=2.0d0*q0*m+(q1-q2)*0.5d0*(rr13-rr23)/m
    b1=(b+rr13/rr11*c)/r1
    b2=(b+rr23/rr22*c)/r2
    c1=c/rr11
    c2=c/rr22
    r1m=0.5d0*(rr11-rr12)
    r2m=0.5d0*(rr12-rr22)
    r3m=0.5d0*(rr13-rr23)
    v1=r1m/(m*r1)
    v2=r2m/(m*r2)
    u1=(r1m*rr13/rr11-r3m)/m
    u2=(r2m*rr23/rr22-r3m)/m
    res=nrm*(intg(s1u,w1,v1,u1,a1,b1,c1) &
        -intg(s1d,w1,v1,u1,a1,b1,c1)     &
        -intg(s2u,w2,v2,u2,a2,b2,c2)     &
        +intg(s2d,w2,v2,u2,a2,b2,c2))
contains
FUNCTION intg( x, w, v, u, aa, bb, cc) RESULT(RES) 
    REAL(8),INTENT(IN) :: x,w,v,u,aa,bb,cc
    REAL(8) :: RES,v1p,v1m,vv,dd,d,s,sp,sm,b,a,c,ln1,ln2,ln3
    REAL(8) :: at,ii,jj,dduu,v1m2,v1p2,v4,kk
     v1p=1.0d0+v
     v1m=1.0d0-v
     vv=v1p*v1m
     dd=w**2*vv-u**2
     d=sqrt(dd)
     s=sqrt(x*x+w*w)
     sp=s+x
     sm=s-x
     b=s+v*x+u
     a=d*(s*v+x)
     c=dd+u*b
     ln1=log(sp)
     ln2=log(sm)
     ln3=log(b)
     ii= 0.5d0*(x*s+w*w*ln1)
     at=atan2(a,c)
     jj= -x+(x-v*u/vv)*ln3-u/vv*ln2+d/vv*at
     dduu=d**2-u**2
     v1m2=1.0d0/v1m**2
     v1p2=1.0d0/v1p**2
     v4=1.0d0/vv**2
     kk= -0.25d0*x**2                                  &
    +0.25d0*u*(sm/v1m+sp/v1p)                          &
    +0.25d0*(dduu*(1.0d0+v*v)*v4+2.0d0*x**2+w**2)*ln3  & 
    +0.125d0*dduu*(ln2*v1m2+ln1*v1p2)                  &
    +d*u*v*v4*at
    res= aa*ii+bb*jj+cc*kk
END FUNCTION intg
END FUNCTION int_trg
END MODULE int_trgf

