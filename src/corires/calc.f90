!------------------------------------------------------------
! calc
!       - module containing subroutines for calculating
!         Beff
!------------------------------------------------------------
MODULE calc
  USE conv
  USE ints
  USE term

CONTAINS

!------------------------------------------------------------
! calc_states
!       - calculates Beff for given states
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational mode numbering offset
! nstates       : int, number of states
! l2h           : 1D int, labeling -> harmonic numbering
! states        : 2D int, states to calculated (vQns, state)
! phi2          : 1D real*8, quadratic force constants
! phi3          : 3D real*8, cubic force constants
! Be            : 1D real*8, rotational constants
! zeta          : 3D real*8, coriolis zetas
! mu1           : 3D real*8, order 1 mu terms (vib,rot,rot)
! mu2           : 4D real*8, order 2 mu terms (vib,vib,rot,rot)
SUBROUTINE calc_states(nvib,voff,nstates,l2h,states,phi2,&
                                    phi3,Be,zeta,mu1,mu2,didq)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3,zeta,mu1,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2,Be
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: states
  INTEGER, DIMENSION(0:), INTENT(IN) :: l2h
  INTEGER, INTENT(IN) :: nvib,voff,nstates
  REAL(KIND=8), DIMENSION(0:2,0:nstates-1) :: Beff
  INTEGER, DIMENSION(0:nvib-1) :: bra,ket
  REAL(KIND=8), DIMENSION(0:2) :: Beff_0
  REAL(KIND=8) :: val
  INTEGER :: i,j,a,b,n

  !Be
  bra = 0
  ket = 0
  Beff_0 = Be
  WRITE(*,*) "                                                 X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  WRITE(*,'(2x,A2)',ADVANCE='no') "Be"
  WRITE(*,'(8x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff_0(2))

  !TESTING TESTING
  WRITE(*,*) "TESTING TESTING TESTING TESTING"
  !term 1 test
  DO a=0,2
    DO b=0,2
      DO i=0,nvib-1
        Beff_0(a) = Beff_0(a) + Be(a)*Be(a)*0.75D0*didq(i,a,b)**2.0D0*2.0D0*Be(b)
        Beff_0(a) = Beff_0(a) + Be(a)*Be(a)*0.75D0*didq(i,a,b)**2.0D0*Be(b)
      END DO
    END DO 
  END DO

  DO a=0,2
    DO i=0,nvib-1
      DO j=0,nvib-1
        Beff_0(a) = Beff_0(a) - Be(a)*Be(a)*0.5D0*zeta(i,j,a)*zeta(i,j,a)*&
                    (phi2(i) - phi2(j))**2.0D0/(phi2(i)*phi2(j)*(phi2(i) + phi2(j)))
      END DO
    END DO
  END DO

  DO a=0,2
    DO i=0,nvib-1
      DO j=0,nvib-1
        Beff_0(a) = Beff_0(a) + Be(a)*Be(a)*0.5D0*phi3(i,i,j)*didq(j,a,a)/phi2(j)**1.5D0
      END DO
    END DO
  END DO

  WRITE(*,'(2x,A2)',ADVANCE='no') "B0"
  WRITE(*,'(8x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff_0(2))

  !Be - B0
  WRITE(*,'(2x,A7)',ADVANCE='no') "Be - B0"
  WRITE(*,'(3x,F18.7)',ADVANCE='no') conv_cm2MHz(Be(0) - Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Be(1) - Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Be(2) - Beff_0(2))

  STOP
  !END TESTING TESTING

  !B0
  !term1
  DO a=0,2
    DO b=0,2
      Beff_0(a) = Beff_0(a) + term_1(nvib,a,b,ket,mu2)*Be(a)*Be(a)
    END DO
  END DO

  !H(1)H(1) term
  DO a=0,2
    Beff_0(a) = Beff_0(a) + term_2(nvib,a,ket,Be,phi2,zeta)*Be(a)*Be(a)
  END DO

  !H(1)H(3) term
  DO a=0,2
    DO b=0,2
    Beff_0(a) = Beff_0(a) + term_3(nvib,a,b,ket,phi2,phi3,mu1)*Be(a)*Be(a)
    END DO
  END DO


  WRITE(*,'(2x,A2)',ADVANCE='no') "B0"
  WRITE(*,'(8x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff_0(2))

  !Be - B0
  WRITE(*,'(2x,A7)',ADVANCE='no') "Be - B0"
  WRITE(*,'(3x,F18.7)',ADVANCE='no') conv_cm2MHz(Be(0) - Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Be(1) - Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Be(2) - Beff_0(2))

  WRITE(*,*)
  WRITE(*,*)

  IF (nstates .EQ. 0) RETURN

  !calculate Be's per state
  WRITE(*,*) " State                                           X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  DO n=0,nstates-1
 
    !get harmonic numbering
    DO j=0,nvib-1
      ket(l2h(j)-1) = states(j,n)
    END DO

    !order 0 term
    Beff(0:2,n) = Be

    !term1
    bra = ket
    DO a=0,2
      DO b=0,2
        Beff(a,n) = Beff(a,n) + term_1(nvib,a,b,ket,mu2)*Be(a)*Be(a)
      END DO
    END DO

    !term 2
    DO a=0,2
      Beff(a,n) = Beff(a,n) + term_2(nvib,a,ket,Be,phi2,zeta)*Be(a)*Be(a)
    END DO

    !term 3
    DO a=0,2
      DO b=0,2
        Beff(a,n) = Beff(a,n) + term_3(nvib,a,b,ket,phi2,phi3,mu1)*Be(a)*Be(a)
      END DO
    END DO
    
    !print
    WRITE(*,'(1x,999(I2,1x))',ADVANCE='no') states(0:nvib-1,n)
    WRITE(*,'(4x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(0,n))
    WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(1,n))
    WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff(2,n))
    
  END DO
  
  WRITE(*,*)
  WRITE(*,*)
  
  !calculate Beff transitions
  WRITE(*,*) " B_eff(v0) --> B_eff(vi)                         X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  DO n=0,nstates-1

    !get harmonic numbering
    DO j=0,nvib-1
      ket(l2h(j)-1) = states(j,n)
    END DO

    !Print
    WRITE(*,'(1x,999(I2,1x))',ADVANCE='no') states(0:nvib-1,n)
    WRITE(*,'(4x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(0,n) - Beff_0(0))
    WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(1,n) - Beff_0(1))
    WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff(2,n) - Beff_0(2))

  END DO

  WRITE(*,*)
  WRITE(*,*)

END SUBROUTINE calc_states
!------------------------------------------------------------

END MODULE calc
!------------------------------------------------------------
