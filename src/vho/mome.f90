!------------------------------------------------------------
! mome  
!       - module containing subroutines for dealing with 
!         momentum terms 
!------------------------------------------------------------
MODULE mome

CONTAINS
!------------------------------------------------------------
! mome_HO_eval
!       - evaluates contribution of momentum
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! qPhi          : int, quantum number of FC
! Phi           : real*8, quadratic force constant
! P2int         : 2D real*8, <i|p^2|i'> type integrals
! momeval       : real*8, value 
SUBROUTINE mome_HO_eval(ndim,PsiL,PsiR,qPhi,Phi,P2int,momeval)
  IMPLICIT NONE         
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: P2int
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: momeval
  REAL(KIND=8), INTENT(IN) :: Phi
  INTEGER, INTENT(IN) :: ndim,qPhi
  INTEGER :: i,j          
  i = qPhi
  !delta function for noninvolved dimensions
  IF (ALL(PsiL(0:i-1) .EQ. PsiR(0:i-1)) .AND.&
      ALL(PsiL(i+1:ndim-1) .EQ. PsiR(i+1:ndim-1))) THEN
    ! p^2 can be v,v and v+2,v
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      j = PsiR(i)
      momeval = momeval + 0.5*Phi*P2int(2*j,i)
    ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
      j = PsiR(i)
      momeval = momeval + 0.5*Phi*P2int(2*j+1,i) 
    END IF              
  END IF
END SUBROUTINE mome_HO_eval
  
!------------------------------------------------------------

END MODULE mome
!------------------------------------------------------------
