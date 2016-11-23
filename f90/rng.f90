!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shimon Panfil: Industrial Physics and Simulations                   !!
! http://industrialphys.com                                           !!
! THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RNG
    IMPLICIT NONE
    PUBLIC
CONTAINS
!       BASED ON GFORTRAN.PDF        
    SUBROUTINE INIT_RANDOM_SEED(CLOCK)
    INTEGER,INTENT(INOUT):: CLOCK
    INTEGER :: I, N 
    INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
    CALL RANDOM_SEED(SIZE = N)
    ALLOCATE(SEED(N))
    IF(CLOCK<0) CALL SYSTEM_CLOCK(COUNT=CLOCK)
    SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
    CALL RANDOM_SEED(PUT = SEED)
    DEALLOCATE(SEED)
    END SUBROUTINE INIT_RANDOM_SEED
END MODULE RNG
