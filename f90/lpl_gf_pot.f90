!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations        !
! http://industrialphys.com                                !
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE LPL_GF_POT
USE LPL_GF
USE T_SURF
IMPLICIT NONE
PUBLIC
CONTAINS

FUNCTION LPL_GF_POT_L1(DST,S,Q) RESULT (RES)
    REAL(F32)::DST(3)
    TYPE(TSURF_TYPE):: S
    REAL(F64)::Q(:),RES
    INTEGER::T
    INTEGER::N1,N2,N3

    RES=ZERO_F64
    DO T=1,S%NT
        N1=S%TVEC(1,T)
        N2=S%TVEC(2,T)
        N3=S%TVEC(3,T)
        RES=RES+LPL_GF_L1(DST,S%VVEC(:,N1),S%VVEC(:,N2),S%VVEC(:,N3),Q(N1),Q(N2),Q(N3))
    ENDDO
ENDFUNCTION LPL_GF_POT_L1

ENDMODULE LPL_GF_POT
