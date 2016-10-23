!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM farr_tst
    USE  farr
    IMPLICIT NONE
    REAL(F32),ALLOCATABLE,DIMENSION(:,:,:,:,:,:,:) ::Arr
    integer::flag=0,R,ar
    character(1)::numstr 

 !   integer::s(7)=[2,3,4,5,6,7,8]
 !   allocate(arr(s(1),s(2),s(3),s(4),s(5),s(6),s(7)),stat=flag)
 !   print *,flag
    do r=1,7
        write(numstr,'(i1)') r
        call READ_ARR7_F32(Arr,ar,'arr_'//numstr//'.bin',FLAG)

        if(flag/=0) then
            print*,'error in reading:','arr_'//numstr//'.bin'
            stop
        endif
        print *,ar,shape(arr)
        call write_ARR7_F32(Arr,ar,'farr_'//numstr//'.bin',FLAG)
        if(flag/=0) then
            print*,'error in writing:','farr_'//numstr//'.bin'
            stop
        endif
    enddo
    print*,Arr(2,2,2,2,2,2,2)

END PROGRAM farr_tst
    

