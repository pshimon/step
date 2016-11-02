!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM data_buf_tst
    USE data_buf
    IMPLICIT NONE
    REAL(F32) ::Arr7(2,3,4,5,6,7,8)
    character(1)::numstr
    integer::flag=0,R
    r=7
    write(numstr,'(i1)') r
    call READ_ARR_7_F32('DataBufFlt'//numstr//'.bin',FLAG,arr7,2,3,4,5,6,7,8)
    if(flag/=0) then
        print*,'error in reading:','DataBufFlt'//numstr//'.bin','flag',flag
    else
        call write_ARR_7_F32('farrFlt_'//numstr//'.bin',FLAG,arr7,2,3,4,5,6,7,8)
    endif
    print*,Arr7(2,2,2,2,2,2,2)

    r=6
    write(numstr,'(i1)') r
    call READ_ARR_6_F32('DataBufFlt'//numstr//'.bin',FLAG,arr7(:,:,:,:,:,:,1),2,3,4,5,6,7)
    if(flag/=0) then
        print*,'error in reading:','DataBufFlt'//numstr//'.bin','flag',flag
    else
        call write_ARR_6_F32('farrFlt_'//numstr//'.bin',FLAG,arr7(:,:,:,:,:,:,1),2,3,4,5,6,7)
    endif
    print*,Arr7(2,2,2,2,2,2,1)

    r=5
    write(numstr,'(i1)') r
    call READ_ARR_5_F32('DataBufFlt'//numstr//'.bin',FLAG,arr7(:,:,:,:,:,1,1),2,3,4,5,6)
    if(flag/=0) then
        print*,'error in reading:','DataBufFlt'//numstr//'.bin','flag',flag
    else
        call write_ARR_5_F32('farrFlt_'//numstr//'.bin',FLAG,arr7(:,:,:,:,:,1,1),2,3,4,5,6)
    endif
    print*,Arr7(2,2,2,2,2,1,1)

    r=4
    write(numstr,'(i1)') r
    call READ_ARR_4_F32('DataBufFlt'//numstr//'.bin',FLAG,arr7(:,:,:,:,1,1,1),2,3,4,5)
    if(flag/=0) then
        print*,'error in reading:','DataBufFlt'//numstr//'.bin','flag',flag
    else
        call write_ARR_4_F32('farrFlt_'//numstr//'.bin',FLAG,arr7(:,:,:,:,1,1,1),2,3,4,5)
    endif
    print*,Arr7(2,2,2,2,1,1,1)

    r=3
    write(numstr,'(i1)') r
    call READ_ARR_3_F32('DataBufFlt'//numstr//'.bin',FLAG,arr7(:,:,:,1,1,1,1),2,3,4)
    if(flag/=0) then
        print*,'error in reading:','DataBufFlt'//numstr//'.bin','flag',flag
    else
        call write_ARR_3_F32('farrFlt_'//numstr//'.bin',FLAG,arr7(:,:,:,1,1,1,1),2,3,4)
    endif
    print*,Arr7(2,2,2,1,1,1,1)
    r=2
    write(numstr,'(i1)') r
    call READ_ARR_2_F32('DataBufFlt'//numstr//'.bin',FLAG,arr7(:,:,1,1,1,1,1),2,3)
    if(flag/=0) then
        print*,'error in reading:','DataBufFlt'//numstr//'.bin','flag',flag
    else
        call write_ARR_2_F32('farrFlt_'//numstr//'.bin',FLAG,arr7(:,:,1,1,1,1,1),2,3)
    endif
    print*,Arr7(2,2,1,1,1,1,1)
     r=1
    write(numstr,'(i1)') r
    call READ_ARR_1_F32('DataBufFlt'//numstr//'.bin',FLAG,arr7(:,1,1,1,1,1,1),2)
    if(flag/=0) then
        print*,'error in reading:','DataBufFlt'//numstr//'.bin','flag',flag
    else
        call write_ARR_1_F32('farrFlt_'//numstr//'.bin',FLAG,arr7(:,1,1,1,1,1,1),2)
    endif
    print*,Arr7(2,1,1,1,1,1,1)

END PROGRAM data_buf_tst
    

