!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MK_SPHERE
    USE TSURF
    IMPLICIT NONE
    TYPE(TSURF_TYPE):: S0,S1
    INTEGER::FLAG,j
    character(1)::num
    character(30)::str
    CALL ICOSAHEDRON(S0)
    do j=0,2
        write(num,'(i1)') 2*j
        str='USPH_'//num//'.tsb'
        CALL WRITE_TSURF_B(S0,trim(str))
        print *,str,s0%nt,s0%nv
        CALL REFINE_TRIANGULATION_2(S0,S1,FLAG)
        IF(FLAG/=0) THEN
            PRINT*,'REFINE_TRIANGULATION_2:',flag
            STOP
        ENDIF
        call MAKE_UNIT_SPHERE(S1)
        write(num,'(i1)') 2*j+1
        str='USPH_'//num//'.tsb'
         print *,str,s1%nt,s1%nv
        CALL WRITE_TSURF_B(S1,trim(str))
        CALL REFINE_TRIANGULATION_2(S1,S0,FLAG)
        IF(FLAG/=0) THEN
            PRINT*,'REFINE_TRIANGULATION_2:',flag
            STOP
        ENDIF
        call MAKE_UNIT_SPHERE(S0)
       
    enddo

ENDPROGRAM  MK_SPHERE





