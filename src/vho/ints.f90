!------------------------------------------------------------
! ints
!       - module supporting integral calculations
!------------------------------------------------------------
MODULE ints
  USE val

CONTAINS
!------------------------------------------------------------
! ints_key
!       - generate the proper key for determining quantum
!         numbers
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, nubmer of basis functions per d
! key           : 1D int, key
! error         : int, exit code

SUBROUTINE ints_key(ndim,nbas,key,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: key
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i
  error = 0
  key(ndim-1) = 1
  DO i=ndim-2,0,-1
    key(i) = key(i+1)*nbas(i+1)
  END DO
END SUBROUTINE ints_key

!------------------------------------------------------------
! ints_qnum
!       - returns quantum number array given the proper key 
!       - this assumes the product array has the final index
!         as the inner most loop
!------------------------------------------------------------
! ndim          : int, number of dimensions
! it            : int, location in 1D product array
! nbas          : 1D int, basis functons per dimension
! key           : 1D int, key vector
! qnum          : 1D int, quantum number vector
! error         : int, exti code

SUBROUTINE ints_qnum(ndim,it,nbas,key,qnum,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: qnum
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,key
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,it
  INTEGER :: i
  error = 0
  DO i=0,ndim-1
    qnum(i) = MOD((it)/key(i),nbas(i))
  END DO 
END SUBROUTINE ints_qnum

!------------------------------------------------------------
! ints_Norm
!       - calculates normalization integrals
!------------------------------------------------------------
! nabs          : int, number of abscissa
! q             : 1D real*8, abscissa
! W             : 1D real*8, weights
! HR            : 1D real*8, hermite polynomials
! norm          : real*8, normalization 
! error         : int, exit code

SUBROUTINE ints_norm(nabs,W,HR,norm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: W,HR
  REAL(KIND=8), INTENT(INOUT) :: norm
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nabs
  INTEGER :: i
  error = 0
  norm = 0.0D0
  DO i=0,nabs-1
    norm = norm + W(i)*HR(i)**2.0D0
  END DO 
  norm = 1.0D0/SQRT(norm)
END SUBROUTINE ints_norm

!------------------------------------------------------------
! ints_normalize
!       - normalize an integral 
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, left hand quantum numbers
! PsiR          : 1D int, right hand quantum numbers
! norm          : 1D real*8, list of normalization constants
! Hij           : real*8, integral to normalize
! error         : int, exit code

SUBROUTINE ints_normalize(ndim,PsiL,PsiR,Hij,norm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: norm
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL, PsiR
  REAL(KIND=8), INTENT(INOUT) :: Hij
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j
  error = 0
  DO i=0,ndim-1
    Hij = Hij*norm(PsiL(i))*norm(PsiR(i))
  END DO
END SUBROUTINE ints_normalize

!------------------------------------------------------------
! ints_normcalc
!       - precalculates all needed normalization constants
!------------------------------------------------------------
! nabs          : int, nubmer of abscissa
! nbas          : 1D int, number of basis functions
! W             : 1D real*8, weights
! Herm          : 2D real*8, hermite polynomials
! norm          : 1D real*8, normalization constants
! error         : int, exit code

SUBROUTINE ints_normcalc(nabs,nbas,W,Herm,norm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Herm
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: norm
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nabs
  REAL(KIND=8), DIMENSION(0:nabs-1) :: HR
  INTEGER :: i
  error = 0
  DO i=0,MAXVAL(nbas)-1
    HR = Herm(i,0:nabs-1)
    CALL ints_norm(nabs,W,HR,norm(i),error)
    IF (error .NE. 0) RETURN
    CALL val_check(norm(i),error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "ints_normcalc  : ERROR"
      WRITE(*,*) "This normalization constant had a bad value", i
      RETURN
    END IF
  END DO
END SUBROUTINE ints_normcalc
!------------------------------------------------------------
! ints_VTcalc
!       - precalculates all needed potential/kinetic 
!         integrals from the V.in files
!       - integrals are stored as:
!           VTint(i,j,k), where k is the dimension, and 
!                         i,j are the LHS and RHS quantum num
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nabs          : int, number of abscissa
! nbas          : 1D int, number of basis functions
! q             : 1D real*8, abscissa
! W             : 1D real*8, weights
! basK          : 1D real*8, basis function force constants
! norm          : 1D real*8, normalization constants
! Herm          : 2D real*8, hermite polynomials [order,abscissa] 
! Vij           : 2D real*8, potential energy [abscissa,dimension]
! VTint         : 3D real*8, precalcualted integrals
! error         : int, exit code

SUBROUTINE ints_VTcalc(ndim,nabs,nbas,q,W,basK,norm,Herm,Vij,&
                       VTint,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: VTint
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Herm,Vij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q,W,basK,norm
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nabs

  REAL(KIND=8), DIMENSION(0:nabs-1) :: HR,HL
  INTEGER :: i,j,k
  error = 0
  VTint = 0.0D0
  DO k=0,ndim-1
    DO j=0,nbas(k)-1
      HR = Herm(j,0:nabs-1)
      DO i=j,nbas(k)-1
        HL = Herm(i,0:nabs-1)

        !potential - HO integral
        CALL ints_VTint(nabs,q,W,HL,HR,basK(k),Vij(:,k),VTint(i,j,k),error)
        IF (error .NE. 0) RETURN
        VTint(i,j,k) = VTint(i,j,k)*norm(i)*norm(j)

        ! + HO value (this is the kinetic term)  
        IF (i .EQ. j) VTint(i,j,k) = VTint(i,j,k) &
                      + basK(k)*(1.0D0*i+0.5D0)

        CALL val_check(VTint(i,j,k),error)
        IF (error .NE. 0) THEN
          WRITE(*,*) "ints_VTcalc  : ERROR"
          WRITE(*,*) "Bad potential at i,j,k",i,j,k
          RETURN
        END IF

      END DO
    END DO
  END DO

END SUBROUTINE ints_VTcalc
!------------------------------------------------------------
! ints_VTint
!       - calculates integrals of the potential 
!------------------------------------------------------------
! nabs          : int, nubmer of abscissa
! q             : 1D real*8, abscissa
! W             : 1D real*8, weights
! HL            : 1D real*8, left hermite polynomials
! HR            : 1D real*8, right hermite polynomials
! basK          : 1D real*8, basis set force constants
! Vij           : 1D real*8, potential energy of dimension
! VTint          : real*8, integral to return
! error         : int, exit code

SUBROUTINE ints_VTint(nabs,q,W,HL,HR,basK,Vij,VTint,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q,W,HL,HR,Vij
  REAL(KIND=8), INTENT(INOUT) :: VTint 
  REAL(KIND=8), INTENT(IN) :: basK
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nabs
  INTEGER :: i
  error = 0
  VTint = 0
  DO i=0,nabs-1
    VTint = VTint + W(i)*HL(i)*HR(i)*&
           (Vij(i) - 0.5D0*basK*q(i)**2.0D0)
  END DO
END SUBROUTINE ints_VTint

!------------------------------------------------------------
! ints_VTput
!       - calculates VT contribution to hamiltonian
!       - assumes that the potential and kinetic terms between
!         normal coordinates are orthogonal
!         ie. <0,0,1|V1|0,0,0> is zero
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, left hand QN
! PsiR          : 1D int, right hand QN
! Hij           : real*8, value to add to
! VTint         : 3D real*8, array of VT integrals 
! error         : int, exit code

SUBROUTINE ints_VTput(ndim,PsiL,PsiR,VTint,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: VTint
  REAL(KIND=8), INTENT(INOUT) :: Hij
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: PsiL,PsiR
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: k
  error = 0
  DO k=0,ndim-1
    IF (ALL(PsiL(0:k-1) .EQ. PsiR(0:k-1)) .AND. &
        ALL(PsiL(k+1:ndim-1) .EQ. PsiR(k+1:ndim-1))) THEN
      Hij = Hij + VTint(PsiL(k),PsiR(k),k)
    END IF
  END DO
END SUBROUTINE ints_VTput

!------------------------------------------------------------
! ints_Q1put
!       - puts Q1 integrals into hamiltonian
!------------------------------------------------------------
! nabs          : int, number of abscissa
! PsiL          : 1D int, left hand QN
! PsiR          : 1D int, right hand QN
! Q1int         : 1D real*8, Q1 terms
! Hij           : real*8, value to change
! error         : int, exit code

!SUBROUTINE ints_Q1put(nabs,
!------------------------------------------------------------

END MODULE ints
!------------------------------------------------------------
