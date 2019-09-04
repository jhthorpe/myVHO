!------------------------------------------------------------
! ints
!       - contains subroutines and functions dealing with 
!         integrals in the HO basis
!------------------------------------------------------------
MODULE ints
  USE sort

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
! ints_pq
!       -calculates the value of <i|pq|j> 
!       - Kindly provided by Devin Matthews
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_pq(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (i .EQ. j) THEN
    ints_pq = -0.5D0
  ELSE IF (j .EQ. i+2) THEN
    ints_pq = -0.5D0*SQRT((x+1.0D0)*(x+2.0D0))
  ELSE IF (j+2 .EQ. i) THEN
    ints_pq = 0.5D0*SQRT((x+1.0D0)*(x+2.0D0))
  ELSE
    ints_pq = 0.0D0
  END IF
END FUNCTION ints_pq


!------------------------------------------------------------
! ints_pp
!       -calculates the value of <i|p^2|j> 
!       - Kindly provided by Devin Matthews
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_pp(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (i .EQ. j) THEN
    ints_pp = -1.0D0*x-0.5D0
  ELSE IF (ABS(j-i) .EQ. 2) THEN
    ints_pp = 0.5D0*SQRT((x+1.0D0)*(x+2.0D0))
  ELSE
    ints_pp = 0.0D0
  END IF
END FUNCTION ints_pp

!------------------------------------------------------------
! ints_pqq
!       -calculates the value of <i|pqq|j> 
!       - James coded this one, so it should not be trusted
!         even half as much as Devin's :) 
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_pqq(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (j .EQ. i+1) THEN
    ints_pqq = -1.0D0*(x+3.0D0)*SQRT((x+1.0D0)/8.0D0)
  ELSE IF (j+1 .EQ. i) THEN
    ints_pqq = (x-1.0D0)*SQRT((x+1.0D0)/8.0D0)
  ELSE IF (j .EQ. i+3) THEN
    ints_pqq = -1.0D0*SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)/8.0D0)
  ELSE IF (j+3 .EQ. i) THEN
    ints_pqq = -1.0D0*SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)/8.0D0)
  ELSE
    ints_pqq = 0.0D0
  END IF
END FUNCTION ints_pqq


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
! ints_qp
!       -calculates the value of <i|qp|j> 
!       - Kindly provided by Devin Matthews
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_qp(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (i .EQ. j) THEN
    ints_qp = 0.5D0
  ELSE IF (j .EQ. i+2) THEN
    ints_qp = -0.5D0*SQRT((x+1.0D0)*(x+2.0D0))
  ELSE IF (j+2 .EQ. i) THEN
    ints_qp = 0.5D0*SQRT((x+1.0D0)*(x+2.0D0))
  ELSE
    ints_qp = 0.0D0
  END IF
END FUNCTION ints_qp

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
! ints_qpq
!       -calculates the value of <i|qpq|j> 
!       - Kindly provided by Devin Matthews
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_qpq(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (j .EQ. i+1) THEN
    ints_qpq = -1.0D0*(x+1.0D0)*SQRT((x+1.0D0)/8.0D0)
  ELSE IF (j+1 .EQ. i) THEN
    ints_qpq = (x+1.0D0)*SQRT((x+1.0D0)/8.0D0)
  ELSE IF (j .EQ. i+3) THEN
    ints_qpq = -1.0D0*SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)/8.0D0)
  ELSE IF (j+3 .EQ. i) THEN
    ints_qpq = SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)/8.0D0)
  ELSE
    ints_qpq = 0.0D0
  END IF
END FUNCTION ints_qpq

!------------------------------------------------------------
! ints_qqp
!       -calculates the value of <i|q^2p|j> 
!       - Kindly provided by Devin Matthews
!       - factors of 1/i must be accounted for
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_qqp(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (j .EQ. i+1) THEN
    ints_qqp = -1.0D0*(x-1.0D0)*SQRT((x+1.0D0)/8.0D0)
  ELSE IF (j+1 .EQ. i) THEN
    ints_qqp = (x+3.0D0)*SQRT((x+1.0D0)/8.0D0)
  ELSE IF (j .EQ. i+3) THEN
    ints_qqp = -1.0D0*SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)/8.0D0)
  ELSE IF (j+3 .EQ. i) THEN
    ints_qqp = SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)/8.0D0)
  ELSE
    ints_qqp = 0.0D0
  END IF
END FUNCTION ints_qqp

!------------------------------------------------------------
! ints_qqqq
!       -calculates the value of <i|q^4|j> 
!       - Kindly provided by Devin Matthews
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_qqqq(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (j .EQ. i) THEN
    ints_qqqq = 0.75D0*(2.0D0*x*x+2.0D0*x+1.0D0)
  ELSE IF (ABS(j-i) .EQ. 2) THEN
    ints_qqqq = 0.5D0*(2.0D0*x+3.0D0)*SQRT((x+1.0D0)*(x+2.0D0))
  ELSE IF (ABS(j-i) .EQ. 4) THEN
    ints_qqqq = 0.25D0*SQRT((x+1.0D0)*(x+2.0D0)*(x+3.0D0)&
                   *(x+4.0D0))
  ELSE
    ints_qqqq = 0.0D0
  END IF
END FUNCTION ints_qqqq


!------------------------------------------------------------
! ints_qqqqq
!       -calculates the value of <i|q^5|j> 
!       - Kindly provided by Devin Matthews
!------------------------------------------------------------
REAL(KIND=8) FUNCTION ints_qqqqq(i,j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  REAL(KIND=8) :: x
  x = 1.0D0*MIN(i,j)
  IF (ABS(j-i) .EQ. 1) THEN
    ints_qqqqq = 1.25D0*(2.0D0*x*x+4.0D0*x+3.0D0)*&
                    SQRT(0.5D0*(x+1.0D0))
  ELSE IF (ABS(j-i) .EQ. 3) THEN
    ints_qqqqq = 1.25D0*(x+2.0D0)*SQRT(0.5D0*&
                    (x+1.0D0)*(x+2.0D0)*(x+3.0D0))
  ELSE IF (ABS(j-i) .EQ. 5) THEN
    ints_qqqqq = 0.25D0*SQRT(0.5D0*(x+1.0D0)*&
                    (x+2.0D0)*(x+3.0D0)*(x+4.0D0)*&
                    (x+5.0D0))
  ELSE
    ints_qqqqq = 0.0D0
  END IF
END FUNCTION ints_qqqqq

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
! ints_phi3
!       - calculates the value of 
!       <bra |1/6 Σ_{i,j,k} Φ_3 q_i*q_j*q_k |ket>
!       - where i is fixed and bra,ket can only 
!         differ in the i'th index
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! i             : int, fixed vibrational number
! bra           : 1D int, bra quantum numbers
! ket           : 1D int, ket quantum numbers
! phi3          : 3D real*8, cubic force constants 
REAL(KIND=8) FUNCTION ints_phi3(nvib,i,bra,ket,phi3)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,i
  REAL(KIND=8) :: s1
  INTEGER :: j
  s1 = 0.0D0 
  !ijj
  DO j=0,i-1
    s1 = s1 + phi3(i,j,j)*ints_qq(bra(j),ket(j))
  END DO
  DO j=i+1,nvib-1
    s1 = s1 + phi3(i,j,j)*ints_qq(bra(j),ket(j))
  END DO
  s1 = 3.0D0*s1*ints_q(bra(i),ket(i))
  !iii
  s1 = s1 + phi3(i,i,i)*ints_qqq(bra(i),ket(i))
  ints_phi3 = 1.0D0/6.0D0*s1
END FUNCTION ints_phi3

!------------------------------------------------------------
! ints_mu0
!       Evaluates the matrix element
!       
!       <bra | μ0_{α,β} |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, α rotational axis
! b             : int, β rotaitonal axis
! mu0           : 2D real*8, zeroth order mu (rot,rot)
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_mu0(nvib,a,b,mu0,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: mu0
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,a,b
  IF (ALL(bra(0:nvib-1) .EQ. ket(0:nvib-1))) THEN 
    ints_mu0 = mu0(a,b) 
  ELSE 
    ints_mu0 = 0.0D0
  END IF
END FUNCTION ints_mu0

!------------------------------------------------------------
! ints_V0
!       Evalutes the matrix element
!       
!       <bra| 1/2 Σ_i q_i^2 |ket>       
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! phi2          : 1D real*8, quadratic force constants
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_V0(nvib,phi2,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib
  REAL(KIND=8) :: s1
  INTEGER :: i
  s1 = 0.0D0
  DO i=0,nvib-1
    IF (ALL(bra(0:i-1) .EQ. ket(0:i-1)) .AND. &
        ALL(bra(i+1:nvib-1) .EQ. ket(i+1:nvib-1))) THEN 
      s1 = s1 + phi2(i)*ints_qq(bra(i),ket(i)) 
    END IF
  END DO
  ints_V0 = 0.5D0*s1
END FUNCTION ints_V0

!------------------------------------------------------------
! ints_T0
!       Evalutes the matrix element
!       
!       <bra| 1/2 Σ_i p_i^2 |ket>       
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! phi2          : 1D real*8, quadratic force constants
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_T0(nvib,phi2,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib
  REAL(KIND=8) :: s1
  INTEGER :: i
  s1 = 0.0D0
  DO i=0,nvib-1
    IF (ALL(bra(0:i-1) .EQ. ket(0:i-1)) .AND. &
        ALL(bra(i+1:nvib-1) .EQ. ket(i+1:nvib-1))) & 
      s1 = s1 - phi2(i)*ints_pp(bra(i),ket(i)) 
  END DO
  ints_T0 = 0.5D0*s1
END FUNCTION ints_T0

!------------------------------------------------------------
! ints_mu1
!       Evalutes the matrix element
!
!       <bra| Σ_i μ_{α,β}^i q_i |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, α rotational axis
! b             : int, β rotaitonal axis
! mu1           : 3D real*8, first order mu (vib,rot,rot)
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_mu1(nvib,a,b,mu1,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: mu1
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,a,b
  REAL(KIND=8) :: s1
  INTEGER :: i
  s1 = 0.0D0
  DO i=0,nvib-1
    IF (ALL(bra(0:i-1) .EQ. ket(0:i-1)) .AND. &
        ALL(bra(i+1:nvib-1) .EQ. ket(i+1:nvib-1))) &
      s1 = s1 + mu1(i,a,b)*ints_q(bra(i),ket(i))
  END DO
  ints_mu1 = s1
END FUNCTION ints_mu1

!------------------------------------------------------------
! ints_mu0pi
!       Evalutes the matrix element
!
!       <bra| μ_{α,β}.π_{β} |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, α rotational axis
! b             : int, β rotaitonal axis
! mu0           : 2D real*8, zeroth order mu (rot,rot)
! phi2          : 1D real*8, quadratic force constants (vib)
! zeta          : 3D real*8, Coriolis constants (vib,vib,rot)
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_mu0pi(nvib,a,b,mu0,phi2,zeta,&
                                 bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: zeta
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: mu0
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,a,b
  INTEGER, DIMENSION(0:nvib-1) :: diff,loc
  REAL(KIND=8) :: s1
  INTEGER :: i,j,ndiff
  s1 = 0.0D0
  diff = 0
  ndiff = 0
  !get differences and locations
  DO i=0,nvib-1 
    IF (bra(i) .NE. ket(i)) THEN
      loc(ndiff) = i
      diff(ndiff) = ABS(bra(i) - ket(i))
      ndiff = ndiff + 1
    END IF
  END DO
  ! i != j b/c of coriolis zetas
  IF (ndiff .NE. 2 .OR. ANY(diff(0:ndiff-1) .GT. 1)) THEN
    ints_mu0pi = 0.0D0
    RETURN
  ELSE
    i = loc(0)
    j = loc(1) 
    s1 = zeta(i,j,b)*SQRT(phi2(j)/phi2(i))*&
           ints_q(bra(i),ket(i))*ints_p(bra(j),ket(j))
    s1 = s1 + zeta(j,i,b)*SQRT(phi2(i)/phi2(j))*&
           ints_q(bra(j),ket(j))*ints_p(bra(i),ket(i))
  END IF
  ints_mu0pi = mu0(a,b)*s1
END FUNCTION ints_mu0pi

!------------------------------------------------------------
! ints_V1
!       Evalutes the matrix element
!
!       <bra| 1/6 Σ_i,j,k Φ_ijk q_i q_j q_k |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! phi3          : 3D real*8, cubic force constants
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 

REAL(KIND=8) FUNCTION ints_V1(nvib,phi3,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib
  INTEGER, DIMENSION(0:nvib-1) :: diff,loc
  INTEGER, DIMENSION(0:2) :: v
  REAL(KIND=8) :: s1
  INTEGER :: i,j,k,ndiff
  s1 = 0.0D0
  diff = 0
  ndiff = 0
  !check the obvious ones
  DO i=0,nvib-1 
    IF (bra(i) .NE. ket(i)) THEN
      loc(ndiff) = i
      diff(ndiff) = ABS(bra(i) - ket(i))
      ndiff = ndiff + 1
    END IF
  END DO
  !too large of a difference
  IF (ndiff .GT. 3 .OR. MAXVAL(diff(0:ndiff-1)) .GT. 3) THEN
    ints_V1 = 0.0D0
    RETURN
  !go through all phi
  ELSE
    DO k=0,nvib-1
      DO j=0,nvib-1
        DO i=0,nvib-1
          v = [i,j,k] 
          CALL sort_int_ijk(v)
          IF (ALL(bra(0:v(0)-1) .EQ. ket(0:v(0)-1)) .AND.&
              ALL(bra(v(0)+1:v(1)-1) .EQ. ket(v(0)+1:v(1)-1)) .AND. &
              ALL(bra(v(1)+1:v(2)-1) .EQ. ket(v(1)+1:v(2)-1)) .AND. &
              ALL(bra(v(2)+1:nvib-1) .EQ. ket(v(2)+1:nvib-1))) THEN
            s1 = s1 + phi3(i,j,k)*ints_q3(nvib,[i,j,k],bra,ket)
          END IF
        END DO
      END DO
    END DO
  END IF
  ints_V1 = 1.0D0/6.0D0*s1
END FUNCTION ints_V1

!------------------------------------------------------------
! ints_q3
!       Evalutes the matrix element 
!
!       <bra| q_i q_j q_k |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! ids           : 1D int, [i,j,k]
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_q3(nvib,ids,bra,ket)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket,ids
  INTEGER, INTENT(IN) :: nvib
  LOGICAL, DIMENSION(0:3-1,0:3-1) :: isq
  INTEGER, DIMENSION(0:3-1,0:3-1) :: positions
  INTEGER, DIMENSION(0:3-1) :: values,counts
  REAL(KIND=8) :: val,temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 3
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1 .OR. i .EQ. 2) THEN
          isq(0,j) = .TRUE.
        ELSE
          isq(0,j) = .FALSE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1 .OR. i .EQ. 2) THEN
          isq(counts(j),j) = .TRUE.
        ELSE
          isq(counts(j),j) = .FALSE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

 !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> Q
  ! 2 -> Q
  DO i=0,nids-1
    !no we reached the end 
    IF (counts(i) .EQ. 0) THEN
      EXIT
    !q
    ELSE IF (counts(i) .EQ. 1) THEN
      IF (isq(0,i)) THEN
        temp = ints_q(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_q3  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    !qq
    ELSE IF (counts(i) .EQ. 2) THEN
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_qq(bra(values(i)),ket(values(i)))
      ELSE
          WRITE(*,*) "ints_q3  : ERROR"
          WRITE(*,*) "There is a 2 index case that has been missed"
          STOP
      END IF

      !qqq
    ELSE IF (counts(i) .EQ. 3) THEN
      IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i)) THEN
        temp = ints_qqq(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_q3  : ERROR"
        WRITE(*,*) "There is a 3 index case that has been missed"
        STOP
      END IF
    ELSE
      WRITE(*,*) "ints_q3  : ERROR"
      WRITE(*,*) "There is a bad case here"
      STOP
    END IF
    val = val*temp
  END DO
  ints_q3 = val
END FUNCTION ints_q3

!------------------------------------------------------------
! ints_mu2
!       Evaluates the matrix element
!
!       <bra| 1/2 Σ_i,j μ_{α,β}^{i,j} qi qj |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, α rotational axis
! b             : int, β rotaitonal axis
! mu2           : 2D real*8, zeroth order mu (rot,rot)
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_mu2(nvib,a,b,mu2,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,a,b
  INTEGER, DIMENSION(0:nvib-1) :: diff,loc
  REAL(KIND=8) :: s1
  INTEGER :: i,j,ndiff
  s1 = 0.0D0
  diff = 0
  ndiff = 0
  !get differences and locations
  DO i=0,nvib-1 
    IF (bra(i) .NE. ket(i)) THEN
      loc(ndiff) = i
      diff(ndiff) = ABS(bra(i) - ket(i))
      ndiff = ndiff + 1
    END IF
  END DO
  !If bra = ket, we get qi^2
  IF (ndiff .EQ. 0) THEN
    DO i=0,nvib-1
      s1 = s1 + mu2(i,i,a,b)*ints_qq(bra(i),ket(i)) 
    END DO
  !otherwise, two must differ
  ELSE IF (ndiff .EQ. 2) THEN
    i = loc(0)
    j = loc(1)
    s1 = s1 + mu2(i,j,a,b)*ints_q(bra(i),ket(i))*&
              ints_q(bra(j),ket(j))
    s1 = s1 + mu2(j,i,a,b)*ints_q(bra(j),ket(j))*&
              ints_q(bra(i),ket(i))
  ELSE
    ints_mu2 = 0.0D0
    RETURN
  END IF
  ints_mu2 = 0.5D0*s1
END FUNCTION ints_mu2

!------------------------------------------------------------
! ints_mu1pi
!       Evaluates the matrix element
!
!       <bra| Σ_i μ_{α,β}^i qi π_β |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, α rotational axis
! b             : int, β rotational axis
! mu1           : 3D real*8, first order μ matrix
! phi2          : 1D real*8, quadratic force constants
! zeta          : 3D real*8, Coriolis zeta matrix
! bra           : 1D int, LHS vibrational quantum numbers
! ket           : 1D int, RHS vibrational quantum numbers

REAL(KIND=8) FUNCTION ints_mu1pi(nvib,a,b,mu1,phi2,zeta,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: mu1,zeta
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,a,b
  INTEGER, DIMENSION(0:2) :: v
  REAL(KIND=8) :: s1,tol
  INTEGER :: i,j,k
  tol = 1.0D-15
  s1 = 0.0D0
  DO k=0,nvib-1
    DO j=0,nvib-1
      IF (ABS(zeta(j,k,b)) .LT. tol) CYCLE
      DO i=0,nvib-1
        IF (ABS(mu1(i,a,b)) .LT. tol) CYCLE
        v = [i,j,k]
        CALL sort_int_ijk(v) 
        IF (ALL(bra(0:v(0)-1) .EQ. ket(0:v(0)-1)) .AND. &
            ALL(bra(v(0)+1:v(1)-1) .EQ. ket(v(0)+1:v(1)-1)) .AND. &
            ALL(bra(v(1)+1:v(2)-1) .EQ. ket(v(1)+1:v(2)-1)) .AND. &
            ALL(bra(v(2)+1:nvib-1) .EQ. ket(v(2)+1:nvib-1)) ) THEN
          s1 = s1 + mu1(i,a,b)*zeta(j,k,b)*SQRT(phi2(k)/phi2(j))*&
                    ints_q2p(nvib,[i,j,k],bra,ket)
        END IF
      END DO
    END DO
  END DO
  ints_mu1pi = s1
END FUNCTION ints_mu1pi

!------------------------------------------------------------
! ints_q2p
!       Evaluates the matrix element
!
!       <bra| q_i q_j p_k |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! ids           : 1D int, [i,j,k]
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_q2p(nvib,ids,bra,ket)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket,ids
  INTEGER, INTENT(IN) :: nvib
  LOGICAL, DIMENSION(0:3-1,0:3-1) :: isq,isp
  INTEGER, DIMENSION(0:3-1,0:3-1) :: positions
  INTEGER, DIMENSION(0:3-1) :: values,counts
  REAL(KIND=8) :: val,temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 3
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1) THEN
          isq(0,j) = .TRUE.
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE.
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1) THEN
          isq(counts(j),j) = .TRUE.
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE.
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

 !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> Q
  ! 2 -> P
  DO i=0,nids-1
    !no we reached the end 
    IF (counts(i) .EQ. 0) THEN
      EXIT
    !q,p
    ELSE IF (counts(i) .EQ. 1) THEN
      IF (isq(0,i)) THEN
        temp = ints_q(bra(values(i)),ket(values(i)))
      ELSE IF (isp(0,i)) THEN
        temp = ints_p(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_q2p  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    !qq,qp
    ELSE IF (counts(i) .EQ. 2) THEN
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_qq(bra(values(i)),ket(values(i)))
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_qp(bra(values(i)),ket(values(i)))
      ELSE
          WRITE(*,*) "ints_q2p  : ERROR"
          WRITE(*,*) "There is a 2 index case that has been missed"
          STOP
      END IF

      !qqp
    ELSE IF (counts(i) .EQ. 3) THEN
      IF (isq(0,i) .AND. isq(1,i) .AND. isp(2,i)) THEN
        temp = ints_qqp(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_q2p  : ERROR"
        WRITE(*,*) "There is a 3 index case that has been missed"
        STOP
      END IF
    ELSE
      WRITE(*,*) "ints_q2p  : ERROR"
      WRITE(*,*) "There is a bad case here"
    END IF
    val = val*temp
  END DO
  ints_q2p = val 
END FUNCTION ints_q2p

!------------------------------------------------------------
! ints_pimu1
!       Evaluates the matrix element
!
!       <bra| π_α Σ_i μ_{α,β}^i qi |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, α rotational axis
! b             : int, β rotational axis
! mu1           : 3D real*8, first order μ matrix
! phi2          : 1D real*8, quadratic force constants
! zeta          : 3D real*8, Coriolis zeta matrix
! bra           : 1D int, LHS vibrational quantum numbers
! ket           : 1D int, RHS vibrational quantum numbers

REAL(KIND=8) FUNCTION ints_pimu1(nvib,a,b,mu1,phi2,zeta,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: mu1,zeta
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,a,b
  INTEGER, DIMENSION(0:2) :: v
  REAL(KIND=8) :: s1,tol
  INTEGER :: i,j,k
  tol = 1.0D-15
  s1 = 0.0D0
  DO k=0,nvib-1
    DO j=0,nvib-1
      IF (ABS(zeta(j,k,b)) .LT. tol) CYCLE
      DO i=0,nvib-1
        IF (ABS(mu1(i,a,b)) .LT. tol) CYCLE
        v = [i,j,k]
        CALL sort_int_ijk(v) 
        IF (ALL(bra(0:v(0)-1) .EQ. ket(0:v(0)-1)) .AND. &
            ALL(bra(v(0)+1:v(1)-1) .EQ. ket(v(0)+1:v(1)-1)) .AND. &
            ALL(bra(v(1)+1:v(2)-1) .EQ. ket(v(1)+1:v(2)-1)) .AND. &
            ALL(bra(v(2)+1:nvib-1) .EQ. ket(v(2)+1:nvib-1)) ) THEN
          s1 = s1 + mu1(i,a,b)*zeta(j,k,a)*SQRT(phi2(k)/phi2(j))*&
                    ints_eval_qpq(nvib,[i,j,k],bra,ket)
        END IF
      END DO
    END DO
  END DO
  ints_pimu1 = s1
END FUNCTION ints_pimu1

!------------------------------------------------------------
! ints_eval_qpq
!       Evaluates the matrix element
!
!       <bra| p_i q_j q_k |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! ids           : 1D int, [i,j,k]
! bra           : 1D int, LHS vibrational QN 
! ket           : 1D int, RHS vibrational QN 
REAL(KIND=8) FUNCTION ints_eval_qpq(nvib,ids,bra,ket)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket,ids
  INTEGER, INTENT(IN) :: nvib
  LOGICAL, DIMENSION(0:3-1,0:3-1) :: isq,isp
  INTEGER, DIMENSION(0:3-1,0:3-1) :: positions
  INTEGER, DIMENSION(0:3-1) :: values,counts
  REAL(KIND=8) :: val,temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 3
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2) THEN
          isq(0,j) = .TRUE.
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE.
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2) THEN
          isq(counts(j),j) = .TRUE.
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE.
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

 !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> P
  ! 2 -> Q
  DO i=0,nids-1
    !no we reached the end 
    IF (counts(i) .EQ. 0) THEN
      EXIT
    !q,p
    ELSE IF (counts(i) .EQ. 1) THEN
      IF (isq(0,i)) THEN
        temp = ints_q(bra(values(i)),ket(values(i)))
      ELSE IF (isp(0,i)) THEN
        temp = ints_p(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_eval_qpq : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    !pq,qq,qp
    ELSE IF (counts(i) .EQ. 2) THEN
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_qq(bra(values(i)),ket(values(i)))
      ELSE IF (isp(0,i) .AND. isq(1,i)) THEN
        temp = ints_pq(bra(values(i)),ket(values(i)))
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_qp(bra(values(i)),ket(values(i)))
      ELSE
          WRITE(*,*) "ints_eval_qpq  : ERROR"
          WRITE(*,*) "There is a 2 index case that has been missed"
          STOP
      END IF

      !qpq
    ELSE IF (counts(i) .EQ. 3) THEN
      IF (isq(0,i) .AND. isp(1,i) .AND. isq(2,i)) THEN
        temp = ints_qpq(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_eval_qpq  : ERROR"
        WRITE(*,*) "There is a 3 index case that has been missed"
        STOP
      END IF
    ELSE
      WRITE(*,*) "ints_eval_qpq : ERROR"
      WRITE(*,*) "There is a bad case here"
    END IF
    val = val*temp
  END DO
  ints_eval_qpq = val 
END FUNCTION ints_eval_qpq

!------------------------------------------------------------
! ints_pimu0pi
!       Evaluates the matrix element
!
!       <bra| π_α μ_{α,β} π_β |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, α rotational axis
! b             : int, β rotational axis
! mu0           : 2D real*8, zeroth order μ 
! phi2          : 1D real*8, quadratic force constants
! zeta          : 3D real*8, Coriolis zeta 
! bra           : 1D int, LHS vibrational quantum numbers
! ket           : 1D int, RHS vibrational quantum numbers

REAL(KIND=8) FUNCTION ints_pimu0pi(nvib,a,b,mu0,phi2,zeta,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: zeta
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: mu0
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib,a,b
  INTEGER, DIMENSION(0:3) :: v
  REAL(KIND=8) :: s1,tol
  INTEGER :: i,j,k,l
  s1 = 0.0D0
  tol = 1.0D-15
  IF (ABS(mu0(a,b)) .LT. tol) THEN
    ints_pimu0pi = 0.0D0
    RETURN
  END IF 
  DO j=0,nvib-1
    DO i=0,nvib-1
      IF (ABS(zeta(i,j,a)) .LT. tol) CYCLE
      DO l=0,nvib-1
        DO k=0,nvib-1
          IF (ABS(zeta(k,l,b)) .LT. tol) CYCLE
          v = [i,j,k,l] 
          CALL sort_int_ijkl(v)
          IF (ALL(bra(0:v(0)-1) .EQ. ket(0:v(0)-1)) .AND. &
              ALL(bra(v(0)+1:v(1)-1) .EQ. ket(v(0)+1:v(1)-1)) .AND. &
              ALL(bra(v(1)+1:v(2)-1) .EQ. ket(v(1)+1:v(2)-1)) .AND. &
              ALL(bra(v(2)+1:v(3)-1) .EQ. ket(v(2)+1:v(3)-1)) .AND. &
              ALL(bra(v(3)+1:nvib-1) .EQ. ket(v(3)+1:nvib-1)) ) THEN
            s1 = s1 + ints_eval_qpqp(nvib,[i,j,k,l],bra,ket)
          END IF
        END DO
      END DO
    END DO
  END DO
  ints_pimu0pi = s1
END FUNCTION ints_pimu0pi

!------------------------------------------------------------
! ints_eval_qpqp
!       - evaluate 2nd order coriolis integrals in the 
!         harmonic oscillator basis
!       - this assumes the orthogonality checks have 
!         been performed
!------------------------------------------------------------
! ndim          : int, number of dimensions
! bra          : 1D int, LHS quantum numbers
! ket          : 1D int, RHF quantum numbers
! ids           : 1D int, [i,j,k,l] ids
! val           : real*8, value of integrals

REAL(KIND=8) FUNCTION ints_eval_qpqp(ndim,ids,bra,ket)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket,ids
  INTEGER, INTENT(IN) :: ndim
  LOGICAL, DIMENSION(0:4-1,0:4-1) :: isq,isp
  INTEGER, DIMENSION(0:4-1,0:4-1) :: positions
  INTEGER, DIMENSION(0:4-1) :: values,counts
  REAL(KIND=8) :: val,temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 4
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

  !WRITE(*,*) "int ids",ids+7

  !Gather information on positions
  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2) THEN
          isq(0,j) = .TRUE.
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE.
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2) THEN
          isq(counts(j),j) = .TRUE.
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE.
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO
  !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> P
  ! 2 -> Q
  ! 3 -> P
  DO i=0,nids-1
    IF (counts(i) .EQ. 0) THEN
      EXIT

    ELSE IF (counts(i) .EQ. 1) THEN
      !q
      IF (isq(0,i)) THEN
        temp = ints_q(bra(values(i)),ket(values(i)))
      !p
      ELSE IF (isp(0,i)) THEN
        temp = ints_p(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_eval_qpqp  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 2) THEN
      !qq
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_qq(bra(values(i)),ket(values(i)))
      !qp
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_qp(bra(values(i)),ket(values(i)))
      !pq
      ELSE IF (isp(0,i) .AND. isq(1,i)) THEN
        temp = ints_pq(bra(values(i)),ket(values(i)))
      !pp
      ELSE IF (isp(0,i) .AND. isp(1,i)) THEN
        temp = ints_pp(bra(values(i)),ket(values(i)))
      ELSE
        WRITE(*,*) "ints_eval_qpqp  : ERROR"
        WRITE(*,*) "There is a 2 index case that has been missed"
        STOP
      END IF

    ELSE
      WRITE(*,*) "ints_eval_qpqp  : ERROR"
      WRITE(*,*) "There is a bad case here"
      STOP
    END IF

    val = val*temp
  END DO

  !accont for 1/i factor of the two p's
  ints_eval_qpqp = -1.0D0*val
END FUNCTION ints_eval_qpqp

!------------------------------------------------------------
! ints_V2
!       Evalutes the matrix element
!
!       <bra| 1/24 Σ_ijkl Φ_ijkl qi qj qk ql |ket>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! phi4          : 4D real*8, quartic force constants
! bra           : 1D int, LHS vibrational quatum numbers
! ket           : 1D int, RHS vibrational quatum numbers

REAL(KIND=8) FUNCTION ints_V2(nvib,phi4,bra,ket)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: phi4
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket
  INTEGER, INTENT(IN) :: nvib
  INTEGER, DIMENSION(0:3) :: v
  REAL(KIND=8) :: s1,tol
  INTEGER :: i,j,k,l
  s1 = 0.0D0
  tol = 1.0D-15
  DO l=0,nvib-1
    DO k=0,nvib-1
      DO j=0,nvib-1
        DO i=0,nvib-1
          IF (ABS(phi4(i,j,k,l)) .LT. tol) CYCLE
          v = [i,j,k,l] 
          CALL sort_int_ijkl(v)
          IF (ALL(bra(0:v(0)-1) .EQ. ket(0:v(0)-1)) .AND. &
              ALL(bra(v(0)+1:v(1)-1) .EQ. ket(v(0)+1:v(1)-1)) .AND. &
              ALL(bra(v(1)+1:v(2)-1) .EQ. ket(v(1)+1:v(2)-1)) .AND. &
              ALL(bra(v(2)+1:v(3)-1) .EQ. ket(v(2)+1:v(3)-1)) .AND. &
              ALL(bra(v(3)+1:nvib-1) .EQ. ket(v(3)+1:nvib-1)) ) THEN
            s1 = s1 + ints_eval_qqqq(nvib,[i,j,k,l],bra,ket)
          END IF
        END DO
      END DO
    END DO
  END DO
  ints_V2 = 1.0D0/24.0D0*s1
END FUNCTION ints_V2

!------------------------------------------------------------
! ints_eval_qqqq
!       Evalutes the matrix element
!   
!       <bra| q q q q |ket>
!------------------------------------------------------------
! ndim          : int, number of dimensions
! bra          : 1D int, LHS quantum numbers
! ket          : 1D int, RHF quantum numbers
! ids           : 1D int, [i,j,k,l] ids
! val           : real*8, value of integrals

REAL(KIND=8) FUNCTION ints_eval_qqqq(ndim,ids,bra,ket)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: bra,ket,ids
  INTEGER, INTENT(IN) :: ndim
  INTEGER, DIMENSION(0:4-1,0:4-1) :: positions
  INTEGER, DIMENSION(0:4-1) :: values,counts
  REAL(KIND=8) :: val,temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 4
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

  !WRITE(*,*) "int ids",ids+7

  !Gather information on positions
  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

  !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> Q
  ! 2 -> Q
  ! 3 -> Q
  DO i=0,nids-1
    IF (counts(i) .EQ. 0) THEN
      EXIT

    !q
    ELSE IF (counts(i) .EQ. 1) THEN
      temp = ints_q(bra(values(i)),ket(values(i)))

    !qq
    ELSE IF (counts(i) .EQ. 2) THEN
      temp = ints_qq(bra(values(i)),ket(values(i)))

    !qqq
    ELSE IF (counts(i) .EQ. 3) THEN
      temp = ints_qqq(bra(values(i)),ket(values(i)))

    !qqqq
    ELSE IF (counts(i) .EQ. 4) THEN
      temp = ints_qqqq(bra(values(i)),ket(values(i)))

    ELSE
      WRITE(*,*) "ints_eval_qqqq  : ERROR"
      WRITE(*,*) "There is a bad case here"
      STOP
    END IF

    val = val*temp
  END DO
  ints_eval_qqqq = val
END FUNCTION ints_eval_qqqq

!------------------------------------------------------------

END MODULE ints
