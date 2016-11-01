!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! adopted from Numerical Recipies
SUBROUTINE polint(xa,ya,n,x,y,dy)
    use fdefs
    INTEGER:: n
!    INTEGER,PARAMETER::NMAX=10
    REAL(f64):: dy,x,y,xa(n),ya(n)
    INTEGER:: i,m,ns
    REAL(f64):: den,dif,dift,ho,hp,w,c(n),d(n)
    ns=1
    dif=abs(x-xa(1))
    do  i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
            ns=i
            dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
    enddo
    y=ya(ns)
    ns=ns-1
    do  m=1,n-1
        do  i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
            enddo
            if (2*ns.lt.n-m)then
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            endif
            y=y+dy
    enddo
    return
END
