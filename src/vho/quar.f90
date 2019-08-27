!------------------------------------------------------------
! quar
!       - module containing subroutines dealing with 
!         quartic force constants
!------------------------------------------------------------
MODULE quar
  USE sort

CONTAINS
!------------------------------------------------------------
! quad_HO_eval
!       - evaluates contribution of a quadratic force constant
!         to the potential energy in harmonic oscillator 
!         basis
!
!       Integrals are stored like:
!        Q4(3*i,k) -> <i|Qk^4|i>    , Q4(3*i+1,k) -> <i+2|Qk^4|i>
!                                   , Q4(3*i+2,k) -> <i+4|Qk^4|i>
!------------------------------------------------------------
! diag          : bool, if true, compute diagonal terms
! ndim          : int, nubmer of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! qPhi          : 1D int, quartic fc dimension ids
! Phi           : real*8, quartic fc 
! Q1int         : 2D real*8, <i|q|i'> integrals   
! Q2int         : 2D real*8, <i|q^2|i'> integrals   
! Q3int         : 2D real*8, <i|q^3|i'> integrals   
! Q4int         : 2D real*8, <i|q^4|i'> integrals   
! quarval       : real*8, value to iterate

SUBROUTINE quar_HO_eval(diag,ndim,PsiL,PsiR,qPhi,&
                          Phi,Q1int,Q2int,Q3int,Q4int,quarval)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q1int,Q2int,Q3int,Q4int
  REAL(KIND=8), INTENT(INOUT) :: quarval
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,qPhi
  REAL(KIND=8), INTENT(IN) :: Phi
  LOGICAL, INTENT(IN) :: diag
  INTEGER, INTENT(IN) :: ndim
  INTEGER, DIMENSION(0:3) :: v
  INTEGER :: i,j,k,l,m,n,o,p
  i = qPhi(0)
  j = qPhi(1)
  k = qPhi(2)
  l = qPhi(3)
  !Check orthogonality of noninvolved terms
  v = [i,j,k,l]
  CALL sort_int_ijkl(v)
  IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR. &
      ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR. &
      ANY(PsiL(v(1)+1:v(2)-1) .NE. PsiR(v(1)+1:v(2)-1)) .OR. &
      ANY(PsiL(v(2)+1:v(3)-1) .NE. PsiR(v(2)+1:v(3)-1)) .OR. &
      ANY(PsiL(v(3)+1:ndim-1) .NE. PsiR(v(3)+1:ndim-1))) RETURN

  !type 1 : i i i i
  IF (i .EQ. j .AND. j .EQ. k .AND. k .EQ. l) THEN
    IF (.NOT. diag) RETURN
    !i, i
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      m = PsiR(i)
      quarval = quarval + Phi/24.0D0*Q4int(3*m,i)
    !i+2, i
    ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
      m = PsiR(i)
      quarval = quarval + Phi/24.0D0*Q4int(3*m+1,i)
    !i+4,i
    ELSE IF (PsiL(i) .EQ. PsiR(i)+4) THEN
      m = PsiR(i)
      quarval = quarval + Phi/24.0D0*Q4int(3*m+2,i)
    END IF

  !type 2 : i i i l 
  ELSE IF (i .EQ. j .AND. j .EQ. k .AND. k .NE. l) THEN
    ! l+1,l
    IF (PsiL(l) .EQ. PsiR(l)+1) THEN
      !i+1,i
      IF (PsiL(i) .EQ. PsiR(i)+1) THEN
        m = PsiR(i)
        p = PsiR(l)
        quarval = quarval + Phi/6.0D0*Q3int(2*m,i)*Q1int(p,l)
      !i+3,i
      ELSE IF (PsiL(i) .EQ. PsiR(i)+3) THEN
        m = PsiR(i)
        p = PsiR(l)
        quarval = quarval + Phi/6.0D0*Q3int(2*m+1,i)*Q1int(p,l)
      END IF
    END IF

  !type 3 : i j j j
  ELSE IF (i .NE. j .AND. j .EQ. k .AND. k .EQ. l) THEN
    !i+1,i
    IF (PsiL(i) .EQ. PsiR(i)+1) THEN
      !j+1,j
      IF (PsiL(j) .EQ. PsiR(j)+1) THEN
        m = PsiR(i)
        n = PsiR(j)
        quarval = quarval + Phi/6.0D0*Q1int(m,i)*Q3int(2*n,j)
      !j+3,j
      ELSE IF (PsiL(j) .EQ. PsiR(j)+3) THEN
        m = PsiR(i)
        n = PsiR(j)
        quarval = quarval + Phi/6.0D0*Q1int(m,i)*Q3int(2*n+1,j)
      END IF
    END IF

  !type 4 : i i k k
  ELSE IF (i .EQ. j .AND. j .NE. k .AND. k .EQ. l) THEN
    !i,i
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      !k,k
      IF (PsiL(k) .EQ. PsiR(k)) THEN
        m = PsiR(i)
        o = PsiR(k)
        quarval = quarval + Phi/4.0D0*Q2int(2*m,i)*Q2int(2*o,k)
      !k+2,k
      ELSE IF (PsiL(k) .EQ. PsiR(k)+2) THEN
        m = PsiR(i)
        o = PsiR(k)
        quarval = quarval + Phi/4.0D0*Q2int(2*m,i)*Q2int(2*o+1,k)
      END IF
    !i+2,i
    ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
      !k,k
      IF (PsiL(k) .EQ. PsiR(k)) THEN
        m = PsiR(i)
        o = PsiR(k)
        quarval = quarval + Phi/4.0D0*Q2int(2*m+1,i)*Q2int(2*o,k)
      !k+2,k
      ELSE IF (PsiL(k) .EQ. PsiR(k)+2) THEN
        m = PsiR(i)
        o = PsiR(k)
        quarval = quarval + Phi/4.0D0*Q2int(2*m+1,i)*Q2int(2*o+1,k)
      END IF
    END IF

  !type 5 : i i k l
  ELSE IF (i .EQ. j .AND. j .NE. k .AND. k .NE. l) THEN
    !k+1,k l+1,l
    IF (PsiL(k) .EQ. PsiR(k)+1 .AND. PsiL(l) .EQ. PsiR(l)+1) THEN
      !i,i
      IF (PsiL(i) .EQ. PsiR(i)) THEN
        m = PsiR(i)
        o = PsiR(k)
        p = PsiR(l)
        quarval = quarval + Phi/2.0D0*Q2int(2*m,i)*&
                            Q1int(o,k)*Q1int(p,l)
      !i+2,i
      ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
        m = PsiR(i)
        o = PsiR(k)
        p = PsiR(l)
        quarval = quarval + Phi/2.0D0*Q2int(2*m+1,i)*&
                            Q1int(o,k)*Q1int(p,l)
      END IF
    END IF

  !type 6 : i j k k
  ELSE IF (i .NE. j .AND. j .NE. k .AND. k .EQ. l) THEN
    !i,i+1  j,j+1
    IF (PsiL(i) .EQ. PsiR(i)+1 .AND. PsiL(j) .EQ. PsiR(j)+1) THEN
      !k,k
      IF (PsiL(k) .EQ. PsiR(k)) THEN
        m = PsiR(i)
        n = PsiR(j)
        o = PsiR(k)
        quarval = quarval + Phi/2.0D0*Q1int(m,i)*&
                            Q1int(n,j)*Q2int(2*o,k)
      !k+2,k
      ELSE IF (PsiL(k) .EQ. PsiR(k)+2) THEN
        m = PsiR(i)
        n = PsiR(j)
        o = PsiR(k)
        quarval = quarval + Phi/2.0D0*Q1int(m,i)*&
                            Q1int(n,j)*Q2int(2*o+1,k)
       END IF
    END IF

  !type 7 : i i k l 
  ELSE IF (i .EQ. j .AND. j .NE. k .AND. k .NE. l) THEN
    !k,k+1  l,l+1
    IF (PsiL(k) .EQ. PsiR(k)+1 .AND. PsiL(l) .EQ. PsiR(l)+1) THEN
      !i,i
      IF (PsiL(i) .EQ. PsiR(i)) THEN
        m = PsiR(i)
        o = PsiR(k)
        p = PsiR(l)
        quarval = quarval + Phi/2.0D0*Q2int(2*m,i)*Q1int(o,k)*Q1int(p,l)
      !i,i+2
      ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
        m = PsiR(i)
        o = PsiR(k)
        p = PsiR(l)
        quarval = quarval + Phi/2.0D0*Q2int(2*m+1,i)*Q1int(o,k)*Q1int(p,l)
      END IF
    END IF

  !type 8 : i j j l
  ELSE IF (i .NE. j .AND. j .EQ. k .AND. k .NE. l) THEN
    !i,i+1  l,l+1
    IF (PsiL(i) .EQ. PsiR(i)+1 .AND. PsiL(l) .EQ. PsiR(l)+1) THEN
      !j,j
      IF (PsiL(j) .EQ. PsiR(j)) THEN
        m = PsiR(i)
        n = PsiR(j)
        p = PsiR(l)
        quarval = quarval + Phi/2.0D0*Q1int(m,i)*Q2int(2*n,j)*Q1int(p,l)
      !j,j+2
      ELSE IF (PsiL(j) .EQ. PsiR(j)+2) THEN
        m = PsiR(i)
        n = PsiR(j)
        p = PsiR(l)
        quarval = quarval + Phi/2.0D0*Q1int(m,i)*Q2int(2*n+1,j)*Q1int(p,l)
      END IF
    END IF

  !type 9 : i j k k 
  ELSE IF (i .NE. j .AND. j .NE. k .AND. k .EQ. l) THEN
    !i,i+1 j,j+1
    IF (PsiL(i) .EQ. PsiR(i)+1 .AND. PsiL(j) .EQ. PsiR(j)+1) THEN
      !k,k
      IF (PsiL(k) .EQ. PsiR(k)) THEN
        m = PsiR(i)
        n = PsiR(j)
        o = PsiR(k)
        quarval = quarval + Phi/2.0D0*Q1int(m,i)*Q1int(n,j)*Q2int(2*o,k)
      ELSE IF (PsiL(k) .EQ. PsiR(k)+2) THEN
        m = PsiR(i)
        n = PsiR(j)
        o = PsiR(k)
        quarval = quarval + Phi/2.0D0*Q1int(m,i)*Q1int(n,j)*Q2int(2*o+1,k)
      END IF
    END IF

  !type 10 : i j k l
  ELSE IF (i .NE. j .AND. j .NE. k .AND. k .NE. l) THEN
    !i+1,i  j+1,j  k+1,k  l+1,l
    IF (PsiL(i) .EQ. PsiR(i)+1 .AND. PsiL(j) .EQ. PsiR(j)+1 &
        .AND. PsiL(k) .EQ. PsiR(k)+1 .AND. PsiL(l) .EQ. PsiR(l)+1) THEN
      m = PsiR(i)
      n = PsiR(j)
      o = PsiR(k)
      p = PsiR(l)
      quarval = quarval + Phi*Q1int(m,i)*Q1int(n,j)*&
                          Q1int(o,k)*Q1int(p,l)
    END IF
  ELSE
    WRITE(*,*) "quar_HO_eval  : ERROR"
    WRITE(*,*) "James, you missed a case!"
    WRITE(*,*) qPhi
    WRITE(*,*) PsiL,PsiR
  END IF
END SUBROUTINE quar_HO_eval

!------------------------------------------------------------
END MODULE quar
!------------------------------------------------------------
