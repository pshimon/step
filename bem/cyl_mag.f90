!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CYL_MAG
USE FDEFS
IMPLICIT NONE
PUBLIC
REAL(F32)::MAG_R,MAG_H,MAG_POS(3),MAG_DIR(3) !MAGNET PARAMETERS
REAL(F64)::MAG_MZ
REAL(F64),PARAMETER::MAG_NORM=ONEOVER4PI_F64*10**4
contains
real(f64) function dip_pot(dst)
   USE VEC3D
    IMPLICIT NONE
    REAL(F32)::DST(3)   
    REAL(F64)::D
    REAL(F64)::R
    REAL(F32)::RV(3)
    RV=DST-MAG_POS
    d=PI_F32*MAG_R**2*MAG_H*MAG_NORM
    R=LENGTH_F32(RV)
    dip_pot=d*DOT_F32(RV,MAG_DIR)/LENGTH_F32(MAG_DIR)/r**3
endfunction dip_pot
ENDMODULE CYL_MAG
