!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM CMP_TSURF
    USE TSURF
    IMPLICIT NONE
    TYPE(TSURF_TYPE):: S1,S2
    INTEGER::FLAG,cac,dt
    CHARACTER(LEN=100) :: FNAME1,fname2
    real(f32)::dv,dn
    CAC=COMMAND_ARGUMENT_COUNT()
    IF(CAC/=2) THEN
        PRINT *,"USAGE CMP_TSURF surf1 surf2"
        STOP
    ENDIF
    CALL GET_COMMAND_ARGUMENT(1, FNAME1)
    CALL GET_COMMAND_ARGUMENT(2, FNAME2)
    call READ_TSURF_B(S1,TRIM(FNAME1),FLAG)
    call READ_TSURF_B(S2,TRIM(FNAME2),FLAG)
    if((s1%nt/=s2%nt).or.(s1%nv/=s2%nv)) then
    print *,'surfaces are different'
        stop
    endif

    dt=sum(abs(s1%tvec-s2%tvec))
    if(dt/=0) then
        print *,'triangles are different'
    call WRITE_TSURF_A(S1,'CMP-'//TRIM(FNAME1))
    call WRITE_TSURF_A(S2,'CMP-'//TRIM(FNAME2))
       
    endif
    dv=maxval(abs(s1%vvec-s2%vvec))
    dn=maxval(abs(s1%nvec-s2%nvec))
    print*,'maxdiff in vvec:',dv
    print*,'maxdiff in nvec:',dn
ENDPROGRAM  CMP_TSURF





