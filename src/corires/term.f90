MODULE term
  USE ints

CONTAINS

!------------------------------------------------------------
! term_1
!       - calculates contribution of term1 (below) to Beff
!         for a given bra,ket, and rotational state
!
!       term1 = 1/2 <i|1/2 Σ_{r,s} μ_{α,α}^{r,s}*q_r*q_s |i>
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, rotational mode
! ket           : 1D int, ket QN's
! mu2           : 4D real*8, second order mu matrix

REAL(KIND=8) FUNCTION term_1(nvib,a,ket,mu2)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  INTEGER, DIMENSION(0:), INTENT(IN) :: ket
  INTEGER, INTENT(IN) :: nvib,a
  INTEGER, DIMENSION(0:nvib-1) :: bra
  REAL(KIND=8) :: temp
  INTEGER :: j
  temp = 0.0D0
  bra = ket
  DO j=0,nvib-1
    temp = temp + mu2(j,j,a,a)*ints_qq(bra(j),ket(j))
  END DO
  term_1 = 0.25D0*temp

END FUNCTION term_1

!------------------------------------------------------------
! term_2 
!       - calculates contributions of term2 (below) to Beff
!         of a given bra,ket, and rotational state
!
!     term2 = 1/I_α^2 * Σ_i!=k <i|π_α|k><k|π_a|i>/(E_i - E_k)
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! a             : int, rotational axis
! ket           : 1D int, bra quantum numbers
! Be            : 1D real*8, Be's
! phi2          : 1D real*8, quadratic force constants 
! zeta          : 3D real*8, Coriolis zetas (vib,vib,rot)

REAL(KIND=8) FUNCTION term_2(nvib,a,ket,Be,phi2,zeta)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: zeta
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Be,phi2
  INTEGER, DIMENSION(0:), INTENT(IN) :: ket
  INTEGER, INTENT(IN) :: nvib,a
  INTEGER, DIMENSION(0:nvib-1) :: bra
  REAL(KIND=8) :: temp,tol,Ei,Ek
  INTEGER :: i,j
  tol = 1.0D-15
  temp = 0.0D0 

  Ei = 0.0D0
  DO i=0,nvib-1
    Ei = Ei + phi2(i)*(1.0D0*ket(i)+0.5D0) 
  END DO

  !WRITE(*,*) "|",ket(0:nvib-1),">",Ei

  !Each vibrational qn can differ by only +/- 1 
  !(B/C no diagonal coriolis zetas to give 0,+/-2)
  DO j=0,nvib-2
    DO i=j+1,nvib-1
      IF (ABS(zeta(i,j,a)) .LT. tol) CYCLE

      !-1,-1
      bra = ket      
      bra(i) = ket(i) - 1
      bra(j) = ket(j) - 1
      IF (bra(i) .GE. 0 .AND. bra(j) .GE. 0) THEN
        Ek = Ei - phi2(i) - phi2(j)
        !WRITE(*,*) "<",bra(0:nvib-1),"|",Ek
        temp = temp - &
               ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))*&
               ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))/&
               (Ei - Ek)
      END IF

      !-1,+1
      bra = ket      
      bra(i) = ket(i) - 1
      bra(j) = ket(j) + 1
      IF (bra(i) .GE. 0) THEN
        Ek = Ei - phi2(i) + phi2(j)
        !WRITE(*,*) "<",bra(0:nvib-1),"|",Ek
        temp = temp - &
               ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))*&
               ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))/&
               (Ei - Ek)
      END IF

      !+1,-1
      bra = ket      
      bra(i) = ket(i) + 1
      bra(j) = ket(j) - 1
      IF (bra(j) .GE. 0) THEN
        Ek = Ei + phi2(i) - phi2(j)
        !WRITE(*,*) "<",bra(0:nvib-1),"|",Ek
        temp = temp - &
               ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))*&
               ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                       zeta(j,i,a),phi2(i),phi2(j))/&
               (Ei - Ek)
      END IF

      !+1,+1
      bra = ket      
      bra(i) = ket(i) + 1
      bra(j) = ket(j) + 1
      Ek = Ei + phi2(i) + phi2(j)
      !WRITE(*,*) "<",bra(0:nvib-1),"|",Ek
      temp = temp - &
             ints_pi(ket(i),ket(j),bra(i),bra(j),zeta(i,j,a),&
                     zeta(j,i,a),phi2(i),phi2(j))*&
             ints_pi(bra(i),bra(j),ket(i),ket(j),zeta(i,j,a),&
                     zeta(j,i,a),phi2(i),phi2(j))/&
             (Ei - Ek)
       
    END DO 
  END DO

  term_2 = temp*4.0D0*Be(a)*Be(a)
  
END FUNCTION term_2

!------------------------------------------------------------

END MODULE term
