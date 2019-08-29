!------------------------------------------------------------
! calc
!       - module containing subroutines for calculating
!         Beff
!------------------------------------------------------------
MODULE calc
  USE conv

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
                                    phi3,Be,zeta,mu1,mu2)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3,zeta,mu1
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2,Be
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: states
  INTEGER, DIMENSION(0:), INTENT(IN) :: l2h
  INTEGER, INTENT(IN) :: nvib,voff,nstates
  REAL(KIND=8), DIMENSION(0:2,0:nstates-1) :: Beff
  INTEGER, DIMENSION(0:nvib-1) :: bra,ket
  REAL(KIND=8), DIMENSION(0:2) :: Beff_0
  INTEGER :: i,j

  WRITE(*,*) "       State               X        Y        Z"
  WRITE(*,*) "-----------------------------------------------------"

  !Calculate Beff for 0 state
  Beff_0 = Be

  !calculate Be's per state
  DO i=0,nstates-1

    !get harmonic numbering
    DO j=0,nvib-1
      ket(l2h(j)-1) = states(j,i)
    END DO

    !order 0 term
    Beff(0:2,i) = Be

    !linear term

    !coriolis term

    !last term 
    
    WRITE(*,'(1x,999(I2,1x))',ADVANCE='no') ket(0:nvib-1)
    WRITE(*,'(4x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(0,i))
    WRITE(*,'(4x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(1,i))
    WRITE(*,'(4x,F18.7)') conv_cm2MHz(Beff(2,i))
    
  END DO
  
  
  !calculate Beff changes

END SUBROUTINE
!------------------------------------------------------------

END MODULE calc
!------------------------------------------------------------
