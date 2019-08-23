!------------------------------------------------------------
! cubi
!       - module containg subroutines dealing with 
!         cubic force constants
!------------------------------------------------------------
MODULE cubi

CONTAINS

!------------------------------------------------------------
! cubi_HO_eval
!       - evaluates contribution of a cubic force constant
!         to the potential energy in harmonic oscillator 
!         basis
!
!       Integrals are stored like:
!       Q3(2*i,k) -> <i+1|Qk^3|i>  , Q3(2*i+1,k) -> <i+3|Qk^3|i>
!------------------------------------------------------------
! diag          : bool, if true, compute diagonal terms
! ndim          : int, nubmer of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! qPhi          : 1D int, cubic fc dimension ids
! Phi           : real*8, cubic fc number
! Q1int         : 2D real*8, <i|q|i'> integrals   
! Q2int         : 2D real*8, <i|q^2|i'> integrals   
! Q3int         : 2D real*8, <i|q^3|i'> integrals   
! cubival       : real*8, value to iterate

SUBROUTINE cubi_HO_eval(diag,ndim,PsiL,PsiR,qPhi,&
                          Phi,Q1int,Q2int,Q3int,cubival)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q1int,Q2int,Q3int
  REAL(KIND=8), INTENT(INOUT) :: cubival
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,qPhi
  REAL(KIND=8), INTENT(IN) :: Phi
  LOGICAL, INTENT(IN) :: diag
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j,k,n,m,o
  i = qPhi(0)
  j = qPhi(1)
  k = qPhi(2)

  !check orthogonality of uninvolved dimensions
  IF (ANY(PsiL(0:i-1) .NE. PsiR(0:i-1)) .OR. &
      ANY(PsiL(i+1:j-1) .NE. PsiR(i+1:j-1)) .OR. &
      ANY(PsiL(j+1:k-1) .NE. PsiR(j+1:k-1)) .OR. &
      ANY(PsiL(k+1:ndim-1) .NE. PsiR(k+1:ndim-1)) &
  ) RETURN

  !type 1, q^3
  IF (i .EQ. j .AND. j .EQ. k) THEN
    IF (.NOT. diag) RETURN
    !v+1, v
    IF (PsiL(i) .EQ. PsiR(i)+1) THEN
      n = PsiR(i)
      cubival = cubival + Phi/6.0D0*Q3int(2*n,i)
    !v+3, v
    ELSE IF (PsiL(i) .EQ. PsiR(i)+3) THEN
      n = PsiR(i)
      cubival = cubival + Phi/6.0D0*Q3int(2*n+1,i)
    END IF

  !type 2, qi,qj^2
  ELSE IF (i .NE. j .AND. j .EQ. k) THEN
    ! i+1,i
    IF (PsiL(i) .EQ. PsiR(i)+1) THEN
      !j,j
      IF (PsiL(j) .EQ. PsiR(j)) THEN
        n = PsiR(i)
        m = PsiR(j)
        cubival = cubival + Phi/2.0D0*Q1int(n,i)*Q2int(2*m,j)
      !j+2,j
      ELSE IF (PsiL(j) .EQ. PsiR(j)+2) THEN
        n = PsiR(i)
        m = PsiR(j)
        cubival = cubival + Phi/2.0D0*Q1int(n,i)*Q2int(2*m+1,j)
      END IF
    END IF

   !type 3, qi^2,qj
  ELSE IF (i .EQ. j .AND. j .NE. k) THEN
    !k+1,k
    IF (PsiL(k) .EQ. PsiR(k)+1) THEN
      !i,i
      IF (PsiL(i) .EQ. PsiR(i)) THEN
        n = PsiR(i)
        m = PsiR(k)
        cubival = cubival + Phi/2.0D0*Q2int(2*n,i)*Q1int(m,k)
      !i+2,i
      ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
        n = PsiR(i)
        m = PsiR(k)
        cubival = cubival + Phi/2.0D0*Q2int(2*n+1,i)*Q1int(m,k)
      END IF
    END IF

  !type 4, qi,qj,qk
  ELSE IF (i .NE. j .AND. j .NE. k) THEN
    !i+1,i  j+1,j    k+1,k
    IF (PsiL(i) .EQ. PsiR(i)+1 .AND. &
        PsiL(j) .EQ. PsiR(j)+1 .AND. &
        PsiL(k) .EQ. PsiR(k)+1 ) THEN
      n = PsiR(i)
      m = PsiR(j)
      o = PsiR(k)
      cubival = cubival + Phi*Q1int(n,i)*Q1int(m,j)*Q1int(o,k)
    END IF
  ELSE
    WRITE(*,*) "cubi_HO_eval  : ERROR"
    WRITE(*,*) "There is a cubic force constant type "
    WRITE(*,*) "  that has not been coded correctly"
  END IF

END SUBROUTINE cubi_HO_eval
!------------------------------------------------------------

END MODULE cubi
!------------------------------------------------------------
