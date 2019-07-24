!------------------------------------------------------------
! ints
!       - module supporting integral calculations
!------------------------------------------------------------
MODULE ints

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
! ints_Vint
!       - calculates integrals of the potential 
!------------------------------------------------------------
! nabs          : int, nubmer of abscissa
! W             : 1D real*8, weights
! HL            : 1D real*8, left hermite polynomials
! HR            : 1D real*8, right hermite polynomials
! Vij           : 1D real*8, potential energy of dimension
! Vint          : real*8, integral to return
! error         : int, exit code

SUBROUTINE ints_Vint(nabs,W,HL,HR,Vij,Vint,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: W,HL,HR,Vij
  REAL(KIND=8), INTENT(INOUT) :: Vint 
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nabs
  INTEGER :: i
  error = 0
  Vint = 0
  DO i=0,nabs-1
    Vint = Vint + W(i)*HL(i)*HR(i)*Vij(i)
  END DO
END SUBROUTINE ints_Vint
!------------------------------------------------------------
! ints_VT
!       - calculates potential and kinetic integrals from 
!         the V.in files
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nabs          : int, number of abscissa
! nbas          : 1D int, number of basis functions
! PsiL          : 1D int, left hand quantum numbers
! PsiR          : 1D int, right hand quantum numbers
! basK          : 1D real*8, basis function force constants
! q             : 1D real*8, abscissa
! W             : 1D real*8, weights
! Vij           : 2D real*8, potential energies [abscissa,dim] 
! Herm          : 2D real*8, hermitie poly      [abscissa,bas]
! Hij           : real*8, integral to evaluate
! error         : int, exit code

SUBROUTINE ints_VT(ndim,nabs,nbas,PsiL,PsiR,basK,q,W,Vij,Herm,Hij,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Vij,Herm
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basK,q,W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: Hij
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nabs
  INTEGER :: i,j
  error = 0
  Hij = 0.0D0
   

END SUBROUTINE ints_VT
!------------------------------------------------------------
END MODULE ints
!------------------------------------------------------------
