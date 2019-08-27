!------------------------------------------------------------
! quad
!       - module for dealing with quadratic force constants
!------------------------------------------------------------
MODULE quad
  USE sort

CONTAINS

!------------------------------------------------------------
! quad_HO_eval
!       - evaluates contributions of quadratic force 
!         constants to the Hamiltonian in the Harmonic 
!         oscillator basis
!
!       Integrals are stored like:
!       Q2(2*i,k) -> <i|Qk^2|i>    , Q2(2*i+1,k) -> <i+2|Qk^2|i>
!
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! qPhi          : 1D int, quadratic FC quantum numbers
! Phi           : real*8, quadratic FC value
! Q2int         : 2D real*8, <i|q^2|i'> type integrals 
! quadval       : real*8, value to add to 

SUBROUTINE quad_HO_eval(ndim,PsiL,PsiR,qPhi,Phi,&
                            Q2int,quadval)
  IMPLICIT NONE 
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q2int
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: quadval
  REAL(KIND=8), INTENT(IN) :: Phi
  INTEGER, INTENT(IN) :: ndim,qPhi
  INTEGER, DIMENSION(0:1) :: v
  INTEGER :: i,j
  i = qPhi 
  !delta function for noninvolved dimensions
  v = [i,j]
  CALL sort_int_ij(v)
  IF (ALL(PsiL(0:v(0)-1) .EQ. PsiR(0:v(0)-1)) .AND.&
      ALL(PsiL(v(0)+1:v(1)-1) .EQ. PsiR(v(0)+1:v(1)-1)) .AND. &
      ALL(PsiL(v(1)+1:ndim-1) .EQ. PsiR(v(1)+1:ndim-1)) ) THEN
    ! q^2 can be v,v and v+2,v
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      j = PsiR(i)
      quadval = quadval + 0.5D0*phi*Q2int(2*j,i)
    ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
      j = PsiR(i)
      quadval = quadval + 0.5D0*phi*Q2int(2*j+1,i)
    END IF
  END IF
END SUBROUTINE quad_HO_eval

!------------------------------------------------------------

END MODULE quad
!------------------------------------------------------------
