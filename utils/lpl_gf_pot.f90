!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations        !
! http://industrialphys.com                                !
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION LPL_GF_TRG(TRG,DST) RESULT (RES)
    USE  LPL_GF_PINT
    IMPLICIT NONE
    REAL(F32)::TRG(3,3),DST(3)
    REAL(F64)::RES
    RES=POT_TRG(GET_TRGV(TRG),DST-TRG(:,1))
ENDFUNCTION LPL_GF_TRG

FUNCTION LPL_GF_DISK(POS,DIR,RAD,DST) RESULT (RES)
    USE VEC3D
    IMPLICIT NONE
    REAL(F32)::POS(3),DIR(3),RAD,DST(3)    
    REAL(F64)::RES,DISK_POT
    EXTERNAL DISK_POT
    REAL(F64)::R,A,Z
    REAL(F32)::RV(3)
    RV=DST-POS
    Z=DOT_F32(RV,DIR)/LENGTH_F32(DIR)
    R=SQRT(DOT_F32(RV,RV)-Z**2)
    RES=DISK_POT(R,Z,RAD)
ENDFUNCTION LPL_GF_DISK

FUNCTION LPL_GF_tsurf(ts,DST,chrg) RESULT (RES)
    USE  tsurf
    IMPLICIT NONE
    TYPE(TSURF_TYPE):: tS
    REAL(F32)::DST(3)
    REAL(F64)::chrg(:),RES
    REAL(F32)::rv(3),trg(3,3) 
    REAL(F64)::LPL_GF_TRG
    external LPL_GF_TRG
    integer::t,n
    !rv=dst-ts%org
    res=0.0d0
    do t=1,ts%nt
        !call get_trg(trg,ts,t)
        res=res+LPL_GF_TRG(TRG,DST)*chrg(t)
    enddo
ENDFUNCTION LPL_GF_tsurf

