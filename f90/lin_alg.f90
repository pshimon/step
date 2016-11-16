!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE LIN_ALG
    USE DATA_BUF
    IMPLICIT NONE
    PUBLIC
CONTAINS
SUBROUTINE  GETRF_F64(N,A,IPIV,INFO)
    INTEGER::N,IPIV(N),INFO
    REAL(F64)::A(N,N)
    CALL DGETRF(N,N,A,N,IPIV,INFO)
ENDSUBROUTINE  GETRF_F64 
SUBROUTINE  GETRS_F64(N,A,K,B,IPIV,INFO)
    INTEGER::N,K,IPIV(N),INFO
    REAL(F64)::A(N,N),B(N,K)
    CALL DGETRS("N",N,K,A,N,IPIV,B,N,INFO)
ENDSUBROUTINE  GETRS_F64 
ENDMODULE LIN_ALG
