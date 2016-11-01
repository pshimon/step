!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM farr_tst
!    USE  farr
    USE CDATA
    IMPLICIT NONE
    !REAL(F32),ALLOCATABLE,DIMENSION(:,:,:,:,:,:,:) ::Arr
    integer::flag=0,R,ar
    character(1)::numstr 
    TYPE(CDATA_F32)::B
     do r=1,7
        write(numstr,'(i1)') r
      !  call READ_ARR7_F32(Arr,ar,'arr_'//numstr//'.bin',FLAG)
       call READ_CDATA_F32(B,'arr_'//numstr//'.bin',FLAG)
        if(flag/=0) then
            print*,'error in reading:','arr_'//numstr//'.bin'
            stop
        endif
        print *,B%SHP(1)-F32_LBL
       ! call write_ARR7_F32(Arr,ar,'farr_'//numstr//'.bin',FLAG)
        call WRITE_CDATA_F32(B,'farr_'//numstr//'.bin',FLAG)
        if(flag/=0) then
            print*,'error in writing:','farr_'//numstr//'.bin'
            stop
        endif
    enddo
 !   print*,Arr(2,2,2,2,2,2,2)

END PROGRAM farr_tst
    

