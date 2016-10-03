!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Shimon Panfil                                           !
! Copyright (c) Shimon Panfil Industrial Physics and Simulations  !
! http://industrialphys.com                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE four1(dataarr,nn,isig)
    INTEGER isig,nn
    REAL dataarr(2*nn)
    INTEGER i,istep,j,m,mmax,n
    REAL tempi,tempr
    REAL theta,wi,wpi,wpr,wr,wtemp               
    n=2*nn                                                        
    j=1
    do i=1,n,2                  
        if(j.gt.i)then
            tempr=dataarr(j)        
            tempi=dataarr(j+1)
            dataarr(j)=dataarr(i)
            dataarr(j+1)=dataarr(i+1)
            dataarr(i)=tempr
            dataarr(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
            j=j-m
            m=m/2
            goto 1
        endif
        j=j+m
    enddo 
    mmax=2                         
2   if (n.gt.mmax) then            
        istep=2*mmax
        theta=6.28318530717959d0/REAL(isig*mmax)               
        wpr=-2.d0*sin(0.5d0*theta)**2                            
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do  m=1,mmax,2         
            do  i=m,n,istep
                j=i+mmax         
                tempr=wr*dataarr(j)-wi*dataarr(j+1)
                tempi=wr*dataarr(j+1)+wi*dataarr(j)
                dataarr(j)=dataarr(i)-tempr
                dataarr(j+1)=dataarr(i+1)-tempi
                dataarr(i)=dataarr(i)+tempr
                dataarr(i+1)=dataarr(i+1)+tempi
            enddo 
            wtemp=wr              
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
        enddo 
        mmax=istep
        goto 2                        
    endif                         
    return
END





