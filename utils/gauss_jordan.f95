
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Shimon Panfil                                           !
! Copyright (c) Shimon Panfil Industrial Physics and Simulations  !
! http://industrialphys.com                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gauss_jordan(a,b)
    REAL, DIMENSION(:,:), INTENT(INOUT) :: a,b
    INTEGER, DIMENSION(size(a,1)) :: ipiv,indxr,indxc
    INTEGER i,icol,irow,j,k,l,ll,m,n                            
    REAL big,dum,pivinv 
    n=size(a,1)
    m=size(b,2)
    ipiv=0
    icol=0
    irow=0
    do  i=1,n                                 
        big=0.0                                     
        do  j=1,n                            
            if(ipiv(j).ne.1)then                 
                do  k=1,n
                    if (ipiv(k).eq.0) then
                        if (abs(a(j,k)).ge.big)then
                            big=abs(a(j,k))
                            irow=j
                            icol=k
                        endif
                    endif
                enddo 
            endif
        enddo 
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
            do  l=1,n
                dum=a(irow,l)
                a(irow,l)=a(icol,l)
                a(icol,l)=dum
            enddo 
            do  l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
            enddo 
        endif
        indxr(i)=irow                         
        indxc(i)=icol                              
        pivinv=1.0/a(icol,icol)
        a(icol,icol)=1.0
        a(icol,:)=a(icol,:)*pivinv
        b(icol,:)=b(icol,:)*pivinv
        do  ll=1,n                          
            if(ll.ne.icol)then               
                dum=a(ll,icol)
                a(ll,icol)=0.0
                a(ll,:)=a(ll,:)-a(icol,:)*dum
                b(ll,:)=b(ll,:)-b(icol,:)*dum
            endif
        enddo 
    enddo                       
    do  l=n,1,-1                             
        if(indxr(l).ne.indxc(l))then               
            do  k=1,n                           
                dum=a(k,indxr(l))                
                a(k,indxr(l))=a(k,indxc(l))
                a(k,indxc(l))=dum
          enddo 
        endif
    enddo 
    return                                     
END SUBROUTINE gauss_jordan

