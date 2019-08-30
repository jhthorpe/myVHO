!------------------------------------------------------------
! ints
!       - contains subroutines and functions dealing with 
!         integrals in the HO basis
!------------------------------------------------------------
MODULE ints

CONTAINS

!------------------------------------------------------------
! ints_p
!       -calculates the value of <i|p|j> 
!       - Kindly provided by Devin Matthews
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_p(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (j .EQ. i+1) THEN
    ints_p = -1.0D0*SQRT((x+1.0D0)/2.0D0)
  ELSE IF (ABS(j-i) .EQ. 1) THEN
    ints_p = SQRT((x+1.0D0)/2.0D0)
  ELSE
    ints_p = 0.0D0
  END IF
END FUNCTION ints_p

!------------------------------------------------------------
! ints_q
!       -calculates the value of <i|q|j> 
!       - Kindly provided by Devin Matthews
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_q(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (ABS(j-i) .EQ. 1) THEN
    ints_q = SQRT((x+1.0D0)/2.0D0)
  ELSE
    ints_q = 0.0D0
  END IF
END FUNCTION ints_q

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
! ints_qqq
!       -calculates the value of <i|q^3|j> 
!       - Kindly provided by Devin Matthews
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_qqq(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (ABS(j-i) .EQ. 1) THEN
    ints_qqq = 3.0D0*SQRT(0.125D0*(x+1.0D0)*(x+1.0D0)*(x+1.0D0))
  ELSE IF (ABS(j-i) .EQ. 3) THEN
    ints_qqq = 0.5D0*SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)/2.0D0)
  ELSE
    ints_qqq = 0.0D0
  END IF
END FUNCTION ints_qqq

!------------------------------------------------------------
! ints_pi
!       - calculates the value of <bra| π |ket>
!         given the respective bra, ket, and indices that 
!         differ at the i'th and j'th quantum numbers
!------------------------------------------------------------
! brai          : int, bra i'th quantum number
! braj          : int, bra j'th quantum number
! keti          : int, ket i'th quantum number
! ketj          : int, ket j'th quantum number
! zij           : real*8, ζ_ij
! zji           : real*8, ζ_ji
! wi            : real*8, ωi
! wj            : real*8, ωj

REAL(KIND=8) FUNCTION ints_pi(brai,braj,keti,ketj,zij,zji,wi,wj)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: zij,zji,wi,wj
  INTEGER, INTENT(IN) :: brai,braj,keti,ketj
  ints_pi = zij*SQRT(wj/wi)*ints_q(brai,keti)*ints_p(braj,ketj) + &
            zji*SQRT(wi/wj)*ints_q(braj,ketj)*ints_p(brai,keti)
END FUNCTION ints_pi
!------------------------------------------------------------

END MODULE ints
