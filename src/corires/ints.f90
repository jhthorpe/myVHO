!------------------------------------------------------------
! ints
!       - contains subroutines and functions dealing with 
!         integrals in the HO basis
!------------------------------------------------------------
MODULE ints

CONTAINS

!------------------------------------------------------------
! ints_qq
!       -calculates the value of <i|q^2|j> 
!       - Kindly provided by Devin Matthews
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_qq(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (i .EQ. j) THEN
    ints_qq = x+0.5D0
  ELSE IF (ABS(j-i) .EQ. 2) THEN
    ints_qq = 0.5D0*SQRT((x+1.0D0)*(x+2.0D0))
  ELSE
    ints_qq = 0.0D0
  END IF
END FUNCTION ints_qq

!------------------------------------------------------------

END MODULE ints
