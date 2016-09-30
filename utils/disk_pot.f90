!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION  DISK_POT(R,Z,A) RESULT(RES)
    USE FDEFS
    REAL(F64)::R,Z,A,RES,res2
    REAL(F64),parameter::phi0=zero_f64,phi1=half_f64*PI_F64
    INTEGER,PARAMETER:: JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1
    REAL,PARAMETER::EPS=1.e-6
    call qromb(res2)
    res=two_f64*res2
    CONTAINS
    FUNCTION DP(PHI) RESULT(RES)
    REAL(F64)::PHI,RES
    REAL(F64)::AR1,AR2,W2,RC,R0,R1,R2
    RC=R*COS(PHI)
    AR1=A-RC
    AR2=A+RC
    W2=R**2+Z**2
    R0=SQRT(W2)
    R1=SQRT(A**2-TWO_F64*A*RC+W2)
    R2=SQRT(A**2+TWO_F64*A*RC+W2)
    RES=R1+R2-TWO_F64*R0+RC*LOG((R1+AR1)/(R0-RC)*(R0+RC)/(R2+AR2))
    ENDFUNCTION DP
    SUBROUTINE trapzd(s,n)
      INTEGER:: n
      REAL(F64):: s
      INTEGER it,j
      REAL(F64):: del,acc,tnm,x
      if (n.eq.1) then
        s=0.5*(phi1-phi0)*(dp(phi0)+dp(phi1))
      else
        it=2**(n-2)
        tnm=it
        del=(phi1-phi0)/tnm
        x=phi0+0.5*del
        acc=0.
        do  j=1,it
          acc=acc+dp(x)
          x=x+del
        enddo
        s=0.5*(s+(phi1-phi0)*acc/tnm)
      endif
      return
    ENDSUBROUTINE trapzd
    SUBROUTINE qromb(ss)
    REAL(F64):: ss
    INTEGER:: j
    REAL(F64):: dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do j=1,JMAX
        call trapzd(s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
    enddo
    ENDSUBROUTINE qromb
    ENDFUNCTION DISK_POT

    FUNCTION  DISK_POT0(Z,A) RESULT(RES)
    USE FDEFS
    REAL(F64)::Z,A,RES
    REAL(F64),PARAMETER::FCT=TWO_F64*PI_F64;
    RES=FCT*(SQRT(A**2+Z**2)-ABS(Z))
    ENDFUNCTION DISK_POT0

