!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE FDEFS
    USE,INTRINSIC:: ISO_C_BINDING
    IMPLICIT NONE
    PUBLIC
    INTEGER,PARAMETER::F32=KIND(1.0E0)
    INTEGER,PARAMETER::F64=KIND(1.0D0)
!        PI,E  40 DIGITS ENOUGH EVEN FOR QUAD
    REAL(F32),PARAMETER :: PI_F32=                          &
 &           3.141592653589793238462643383279502884197 
    REAL(F64),PARAMETER :: PI_F64=                          & 
 &           3.141592653589793238462643383279502884197D0
    REAL(F32),PARAMETER :: E_F32=                            & 
 &               2.718281828459045235360287471352662497757
    REAL(F64),PARAMETER :: E_F64=                            &
 &                 2.718281828459045235360287471352662497757D0
    REAL(F32),PARAMETER :: TINY_F32=1.0E-30 
    REAL(F32),PARAMETER :: HUGE_F32=1.0E30   
    REAL(F64),PARAMETER :: TINY_F64=1.0D-300
    REAL(F64),PARAMETER :: HUGE_F64=1.0D300 
    REAL(F32),PARAMETER :: SMALL_F32=1.0E-20 
    REAL(F32),PARAMETER :: LARGE_F32=1.0E20   
    REAL(F64),PARAMETER ::SMALL_F64=1.0D-200
    REAL(F64),PARAMETER ::LARGE_F64=1.0D200 

    REAL(F32),PARAMETER :: ZERO_F32=0.0E0
    REAL(F64),PARAMETER :: ZERO_F64=0.0D0
    REAL(F32),PARAMETER :: ONE_F32=1.0E0
    REAL(F64),PARAMETER :: ONE_F64=1.0D0
    REAL(F32),PARAMETER :: TWO_F32=2.0E0
    REAL(F64),PARAMETER :: TWO_F64=2.0D0
    REAL(F32),PARAMETER :: HALF_F32=0.5E0
    REAL(F64),PARAMETER :: HALF_F64=0.5D0
    REAL(F32),PARAMETER :: THIRD_F32=1.0E0/3.0E0
    REAL(F64),PARAMETER :: THIRD_F64=1.0D0/3.0D0
    REAL(F32),PARAMETER :: QUART_F32=0.25E0
    REAL(F64),PARAMETER :: QUART_F64=0.25D0

    COMPLEX(F32),PARAMETER::I_F32=CMPLX(0.0,1.0,F32)
    COMPLEX(F64),PARAMETER::I_F64=CMPLX(0.0,1.0,F64)

    REAL(F32),PARAMETER::ONEOVER4PI_F32=QUART_F32/PI_F32
    REAL(F64),PARAMETER::ONEOVER4PI_F64=QUART_F64/PI_F64
    REAL(F32),PARAMETER::TORAD_F32=PI_F32/180.0_F32
    REAL(F64),PARAMETER::TORAD_F64=PI_F64/180.0_F64
    
END MODULE FDEFS
