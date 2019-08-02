!------------------------------------------------------------
! ints_HO
!       - module supporting integral calculations in the
!         HO basis
!------------------------------------------------------------
MODULE ints_HO
  USE val
  USE gauss
  USE key

CONTAINS

!------------------------------------------------------------
! ints_HO_Norm
!       - calculates normalization integrals
!       - outdated
!------------------------------------------------------------
! nabs          : int, number of abscissa
! q             : 1D real*8, abscissa
! W             : 1D real*8, weights
! HR            : 1D real*8, hermite polynomials
! norm          : real*8, normalization 
! error         : int, exit code

SUBROUTINE ints_HO_norm(nabs,W,HR,norm,error)
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
END SUBROUTINE ints_HO_norm

!------------------------------------------------------------
! ints_HO_normalize
!       - normalize an integral 
!       - outdated
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, left hand quantum numbers
! PsiR          : 1D int, right hand quantum numbers
! norm          : 1D real*8, list of normalization constants
! Hij           : real*8, integral to normalize
! error         : int, exit code

SUBROUTINE ints_HO_normalize(ndim,PsiL,PsiR,Hij,norm,error)
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
END SUBROUTINE ints_HO_normalize

!------------------------------------------------------------
! ints_HO_normcalc
!       - precalculates all needed normalization constants
!       - this function uses either MAXVAL(100,nbas+10) 
!         abscissa to calculate the normalization constants 
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, nubmer of basis funcitons
! nabs          : 1D int, number of abscissa
! W             : 2D real*8, weights
! Herm          : 3D real*8, hermite polynomials [abcs,bas qn,dim]
! norm          : 2D real*8, normalization constants [bas qn,dim]
! error         : int, exit code

SUBROUTINE ints_HO_normcalc(ndim,nbas,nabs,W,Herm,norm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: Herm
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: norm
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hin
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: qin,Win
  CHARACTER(LEN=20) :: fname
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,n,mbas
  error = 0
  mbas = MAXVAL(nbas)
  n = MAX(mbas+10,100) 
  ALLOCATE(qin(0:n-1))
  ALLOCATE(Win(0:n-1))
  ALLOCATE(Hin(0:n-1,0:mbas-1))

  !generate internal gaussian quadrature
  CALL gauss_hermite(n,qin,Win,error)
  IF (error .NE. 0) RETURN

  ! generate internal hermite polynomials
  DO j=0,mbas-1
    IF (j .EQ. 0) THEN
      Hin(0:n-1,0) = 1.0D0
    ELSE IF (j .EQ. 1) THEN
      Hin(0:n-1,1) = 2.0D0*qin(0:n-1)
    ELSE 
      Hin(0:n-1,j) = 2.0D0*qin(0:n-1)*Hin(0:n-1,j-1) &
                     - 2.0D0*(j-1)*Hin(0:n-1,j-2)
    END IF
  END DO

  ! Generate internal normalization
  DO k=0,ndim-1
    DO j=0,nbas(k)-1
      val = 0.0D0
      DO i=0,n-1
        val = val + Win(i)*Hin(i,j)*Hin(i,j) 
      END DO 
      norm(j,k) = SQRT(1.0D0/val)
    END DO
  END DO
  
  DEALLOCATE(qin)
  DEALLOCATE(Win)
  DEALLOCATE(Hin)
  
 ! OLD CODE
 ! DO k=0,ndim-1
 !   DO j=0,nbas(k)-1
 !     val = 0.0D0
 !     DO i=0,nabs(k)-1
 !       val = val + W(i,k)*Herm(i,j,k)*Herm(i,j,k) 
 !     END DO 
 !     norm(j,k) = SQRT(1.0D0/val)
 !   END DO
 ! END DO
END SUBROUTINE ints_HO_normcalc

!------------------------------------------------------------
! ints_HO_VTcalc
!       - precalculates all needed potential/kinetic 
!         integrals from the V.in files
!       - integrals are stored as:
!       VTint(j-i+keyI(i,k),k) -> <j|Vi + Ti|i>
!
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, number of basis functions
! nabs          : 1D int, number of abscissa
! q             : 2D real*8, abscissa
! W             : 2D real*8, weights
! basK          : 1D real*8, basis function force constants
! norm          : 1D real*8, normalization constants
! Herm          : 3D real*8, Herm poly [abscissa,basis qn,dimension]
! keyI          : 2D real*8, intermediate array [#1
! Vij           : 2D real*8, potential energy [abscissa,dimension]
! VTint         : 2D real*8, precalcualted integrals
! error         : int, exit code

SUBROUTINE ints_HO_VTcalc(ndim,nbas,nabs,q,W,basK,norm,Herm,keyI,Vij,&
                       VTint,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: Herm
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: VTint
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Vij,q,W,norm
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basK
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: keyI
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,a
  error = 0
  VTint = 0.0D0
  DO k=0,ndim-1
    !WRITE(*,*) "dimension :",k
    DO j=0,nbas(k)-1
  !    WRITE(*,*) "j=",j
      DO i=j,nbas(k)-1
 !       WRITE(*,*) "i=",i
        val = 0.0D0
        !potential - harmonic 
        DO a=0,nabs(k)-1
          val = val + W(a,k)*Herm(a,i,k)*Herm(a,j,k)*(Vij(a,k) &
                - (0.5D0*basK(k)*q(a,k)**2.0D0))
   !       WRITE(*,*) "a,val",a,val
        END DO
        val = val*norm(i,k)*norm(j,k)
   !     WRITE(*,*) "Normalied val:",val
      
        !kinetic energy part
        IF (i .EQ. j) val = val + bask(k)*(1.0D0*i+0.5D0)
        !WRITE(*,*) i,j,i-j,keyI(j,k),i-j+keyI(j,k) 
        VTint(i-j+keyI(j,k),k) = val

        CALL val_check(val,error)
        IF (error .NE. 0) THEN
          WRITE(*,*) "ints_HO_VTcalc  : ERROR"
          WRITE(*,*) "Bad VT integral at i,j,k",i,j,k
          RETURN
        END IF 
      END DO
    END DO
!    WRITE(*,*) '------------'
  END DO
 

  !DO k=0,ndim-1
  !  DO j=0,nbas(k)-1
  !    HR = Herm(j,0:nabs-1)
  !    DO i=j,nbas(k)-1
  !      HL = Herm(i,0:nabs-1)
  !
  !      !potential - HO integral
  !      CALL ints_HO_VTint(nabs,q,W,HL,HR,basK(k),Vij(:,k),VTint(i,j,k),error)
  !      IF (error .NE. 0) RETURN
  !      VTint(i,j,k) = VTint(i,j,k)*norm(i)*norm(j)
  !
  !      ! + HO value (this is the kinetic term)  
  !      IF (i .EQ. j) VTint(i,j,k) = VTint(i,j,k) &
  !                    + basK(k)*(1.0D0*i+0.5D0)
  !
  !      CALL val_check(VTint(i,j,k),error)
  !      IF (error .NE. 0) THEN
  !        WRITE(*,*) "ints_HO_VTcalc  : ERROR"
  !        WRITE(*,*) "Bad potential at i,j,k",i,j,k
  !        RETURN
  !      END IF

  !    END DO
  !  END DO
  !END DO

END SUBROUTINE ints_HO_VTcalc
!------------------------------------------------------------
! ints_HO_VTint
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

SUBROUTINE ints_HO_VTint(nabs,q,W,HL,HR,basK,Vij,VTint,error)
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
END SUBROUTINE ints_HO_VTint

!------------------------------------------------------------
! ints_HO_VTput
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

SUBROUTINE ints_HO_VTput(ndim,PsiL,PsiR,VTint,Hij,error)
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
END SUBROUTINE ints_HO_VTput

!------------------------------------------------------------
! ints_HO_Q1calc
!       - calculates integrals of the kind <i|qi|i'> 
!
!       Integrals are stored like:
!       Q1(i,k) -> <i+1|Qk|i>
!      
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! Q1int         : 2D real*8, Q1 terms
! error         : int, exit code

SUBROUTINE ints_HO_Q1calc(ndim,nbas,Q1int,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Q1int
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j,N
  error = 0 
  WRITE(*,*) "<i|q|i'> type integrals"
  DO j=0,ndim-1
    DO i=0,nbas(j)-1
      Q1int(i,j) = SQRT(0.5D0*(1.0D0*i+1.0D0))
    END DO
  END DO
END SUBROUTINE ints_HO_Q1calc

!------------------------------------------------------------
! ints_HO_Q2calc
!       - calculates integrals of the kind <i|qi^2|i'> 
!
!       Integrals are stored like:
!       Q2(2*i,k) -> <i|Qk^2|i>    , Q2(2*i+1,k) -> <i+2|Qk^2|i>
!      
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! Q2int         : 2D real*8, Q2 terms
! error         : int, exit code

SUBROUTINE ints_HO_Q2calc(ndim,nbas,Q2int,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Q2int
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j,N
  error = 0 
  WRITE(*,*) "<i|q^2|i'> type integrals"
  DO j=0,ndim-1
    DO i=0,nbas(j)-1
      Q2int(2*i,j) = 1.0D0*i + 0.5D0 
      Q2int(2*i+1,j) = 0.5D0*SQRT((1.0D0*i+1.0D0)*(1.0D0*i+2.0D0))
    END DO
  END DO
END SUBROUTINE ints_HO_Q2calc

!------------------------------------------------------------
! ints_HO_Q3calc
!       - calculates integrals of the kind <i|qi^3|i'> 
!
!       Integrals are stored like:
!       Q3(2*i,k) -> <i+1|Qk^3|i>  , Q3(2*i+1,k) -> <i+3|Qk^3|i>
!      
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! Q3int         : 2D real*8, Q3 terms
! error         : int, exit code

SUBROUTINE ints_HO_Q3calc(ndim,nbas,Q3int,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Q3int
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j,N
  error = 0 
  WRITE(*,*) "<i|q^3|i'> type integrals"
  DO j=0,ndim-1
    DO i=0,nbas(j)-1
      Q3int(2*i,j) = 3.0D0*SQRT(0.125D0*(1.0D0*i+1.0D0)**3.0D0) 
      Q3int(2*i+1,j) = SQRT(0.125D0*(1.0D0*i+1.0D0)*(1.0D0*i+2.0D0)*&
                          (1.0D0*i+3.0D0)) 
    END DO
  END DO
END SUBROUTINE ints_HO_Q3calc

!------------------------------------------------------------
! ints_HO_Q4calc
!       - calculates integrals of the kind <i|qi^4|i'> 
!
!       Integrals are stored like:
!        Q4(3*i,k) -> <i|Qk^4|i>    , Q4(3*i+1,k) -> <i+2|Qk^4|i>
!                                   , Q4(3*i+2,k) -> <i+4|Qk^4|i>
!      
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! Q4int         : 2D real*8, Q4 terms
! error         : int, exit code

SUBROUTINE ints_HO_Q4calc(ndim,nbas,Q4int,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Q4int
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j,N
  error = 0 
  WRITE(*,*) "<i|q^4|i'> type integrals"
  DO j=0,ndim-1
    DO i=0,nbas(j)-1
      Q4int(3*i,j) = 0.75D0*(2.0D0*i**2.0D0 + 2.0D0*i + 1.0D0) 
      Q4int(3*i+1,j) = 0.5D0*(2.0D0*i+3.0D0)*&
                     SQRT((1.0D0*i+1.0D0)*(1.0D0*i+2.0D0)) 
      Q4int(3*i+2,j) = 0.25D0*SQRT((1.0D0*i+1.0D0)*(1.0D0*i+2.0D0)&
                     *(1.0D0*i+3.0D0)*(1.0D0*i+4.0D0)) 
    END DO
  END DO
END SUBROUTINE ints_HO_Q4calc

!------------------------------------------------------------
! ints_HO_P2calc
!       - calculates integrals of the kind <i|pi^2|i'> 
!
!       Integrals are stored like:
!       P2(i,k) -> <2*i|Pk^2|i>    , P2(i,k)   -> <2*i+2|Pk^2|i>
!      
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! P2int         : 2D real*8, P2 terms
! error         : int, exit code

SUBROUTINE ints_HO_P2calc(ndim,nbas,P2int,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: P2int
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j,N
  error = 0 
  WRITE(*,*) "<i|p^2|i'> type integrals"
  DO j=0,ndim-1
    DO i=0,nbas(j)-1
      P2int(2*i,j) = 1.0D0*i + 0.5D0 
      P2int(2*i+1,j) = -0.5D0*SQRT((1.0D0*i+1.0D0)*(1.0D0*i+2.0D0)) 
    END DO
  END DO
END SUBROUTINE ints_HO_P2calc

!------------------------------------------------------------
! ints_HO_polyput
!       - calculates the contributions of the various 
!         force constants to the Hamiltonian in the HO basis
!
!       Integrals are stored like:
!       Q1(i,k) -> <i+1|Qk|i>
!       Q2(2*i,k) -> <i|Qk^2|i>    , Q2(2*i+1,k) -> <i+2|Qk^2|i>
!       Q3(2*i,k) -> <i+1|Qk^3|i>  , Q3(2*i+1,k) -> <i+3|Qk^3|i>
!       Q4(3*i,k) -> <i|Qk^4|i>    , Q4(3*i+1,k) -> <i+2|Qk^4|i>
!                                  , Q4(3*i+2,k) -> <i+4|Qk^4|i>
!       P2(2*i,k) -> <i|Pk^2|i>    , P2(2*i,k)   -> <i+2|Pk^2|i>
!       
!       Where i indicates the i'th quantum number of the 
!         k'th dimension
!
!       The integral
!          < i0 j0 k0 | Φ_122 | i1 j0 k0 >
!       Is then
!           1/6 * Φ_122 * <i1|Qi|i0> * <j0|Qj^2|j0> * <k0|k0>

!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! nQ2           : int, number of quadratic FC
! qQ2           : 1D int, quadratic FC quantum numbers
! Q2            : 1D real*8, quadratic FC values
! nQ3           : int, number of cubic FC
! qQ3           : 1D int, cubic FC quantum numbers
! Q3            : 1D real*8, cubic FC values
! nQ4           : int, number of quartic FC
! qQ4           : 1D int, quartic FC quantum numbers
! Q4            : 1D real*8, quartic FC values
! Q1int         : 2D real*8, <i|q|i'> type integrals 
! Q2int         : 2D real*8, <i|q^2|i'> type integrals 
! Q3int         : 2D real*8, <i|q^3|i'> type integrals 
! Q4int         : 2D real*8, <i|q^4|i'> type integrals 
! P2int         : 2D real*8, <i|p^2|i'> type integrals 
! Hij           : real*8, value of this hamiltonian element
! error         : int, error code

SUBROUTINE ints_HO_polyput(ndim,PsiL,PsiR,nQ2,qQ2,Q2,&
                           nQ3,qQ3,Q3,nQ4,qQ4,Q4,Q1int,Q2int,&
                           Q3int,Q4int,P2int,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q1int,Q2int,Q3int,&
                                                Q4int,P2int
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Q2,Q3,Q4
  INTEGER, DIMENSION(0:), INTENT(IN) :: qQ2,qQ3,qQ4
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: Hij
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nQ2,nQ3,nQ4
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8) :: quadval,cubival,quarval,momeval
  LOGICAL :: diag
  INTEGER :: i,j
  error = 0
  quadval = 0.0D0
  cubival = 0.0D0
  quarval = 0.0D0
  momeval = 0.0D0
  diag = .TRUE.

  !phi_ij (quadratic) terms
  DO i=0,nQ2-1
    CALL ints_HO_quadeval(ndim,PsiL,PsiR,qQ2(i),Q2(i),&
                        Q2int,quadval)     
  END DO

  !phi_ijk (cubic terms)
  DO i=0,nQ3-1
    CALL ints_HO_cubieval(diag,ndim,PsiL,PsiR,qQ3(3*i:3*i+2),&
                          Q3(i),Q1int,Q2int,Q3int,cubival)
  END DO

  !phi_ijkl (quartic terms)
  DO i=0,nQ4-1
    CALL ints_HO_quareval(diag,ndim,PsiL,PsiR,qQ4(4*i:4*i+3),&
                          Q4(i),Q1int,Q2int,Q3int,Q4int,quarval)
  END DO

  !p^2 (momentum) terms
  DO i=0,nQ2-1
    CALL ints_HO_momeval(ndim,PsiL,PsiR,qQ2(i),Q2(i),&
                        P2int,momeval)     
  END DO
  
  Hij = quadval + cubival + quarval + momeval

END SUBROUTINE ints_HO_polyput

!------------------------------------------------------------
! ints_HO_diagput
!       - calculates the contributions of the various 
!         force constants to the Hamiltonian in the HO basis
!
!       Integrals are stored like:
!       Q1(i,k) -> <i+1|Qk|i>
!       Q2(2*i,k) -> <i|Qk^2|i>    , Q2(2*i+1,k) -> <i+2|Qk^2|i>
!       Q3(2*i,k) -> <i+1|Qk^3|i>  , Q3(2*i+1,k) -> <i+3|Qk^3|i>
!       Q4(3*i,k) -> <i|Qk^4|i>    , Q4(3*i+1,k) -> <i+2|Qk^4|i>
!                                  , Q4(3*i+2,k) -> <i+4|Qk^4|i>
!       Where i indicates the i'th quantum number of the 
!         k'th dimension
!
!       The integral
!          < i0 j0 k0 | Φ_122 | i1 j0 k0 >
!       Is then
!           1/6 * Φ_122 * <i1|Qi|i0> * <j0|Qj^2|j0> * <k0|k0>

!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! nQ2           : int, number of cubic FC
! qQ2           : 1D int, cubic FC quantum numbers
! Q2            : 1D real*8, cubic FC values
! nQ3           : int, number of cubic FC
! qQ3           : 1D int, cubic FC quantum numbers
! Q3            : 1D real*8, cubic FC values
! nQ4           : int, number of quartic FC
! qQ4           : 1D int, quartic FC quantum numbers
! Q4            : 1D real*8, quartic FC values
! Q1int         : 2D real*8, <i|q|i'> type integrals 
! Q2int         : 2D real*8, <i|q^2|i'> type integrals 
! Q3int         : 2D real*8, <i|q^3|i'> type integrals 
! Q4int         : 2D real*8, <i|q^4|i'> type integrals 
! Hij           : real*8, value of this hamiltonian element
! error         : int, error code

SUBROUTINE ints_HO_diagput(ndim,PsiL,PsiR,nQ2,qQ2,Q2,&
                           nQ3,qQ3,Q3,nQ4,qQ4,Q4,Q1int,Q2int,&
                           Q3int,Q4int,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q1int,Q2int,Q3int,&
                                                Q4int
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Q2,Q3,Q4
  INTEGER, DIMENSION(0:), INTENT(IN) :: qQ2,qQ3,qQ4
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: Hij
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nQ2,nQ3,nQ4
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8) :: cubival,quarval
  LOGICAL :: diag
  INTEGER :: i,j
  error = 0
  cubival = 0.0D0
  quarval = 0.0D0
  diag = .FALSE.

  !phi_ijk (cubic terms)
  DO i=0,nQ3-1
    CALL ints_HO_cubieval(diag,ndim,PsiL,PsiR,qQ3(3*i:3*i+2),&
                          Q3(i),Q1int,Q2int,Q3int,cubival)
  END DO

  !phi_ijkl (quartic terms)
  DO i=0,nQ4-1
    CALL ints_HO_quareval(diag,ndim,PsiL,PsiR,qQ4(4*i:4*i+3),&
                          Q4(i),Q1int,Q2int,Q3int,Q4int,quarval)
  END DO

  
  Hij = cubival + quarval 

END SUBROUTINE ints_HO_diagput

!------------------------------------------------------------
! ints_HO_quadput
!       - calculates the value of a Hamiltonian element 
!         using gaussian quadrature for all the dimensions
!       
!       - note that the hermite polynomials and norm const are 
!         assumed to be already arranged so that we don't
!         need to access by quantum number. The LHS
!         values are accessed (2*i) and the RHS are (2*i+1)
!
!       - Vq is already in order       
!
!       - the code is pretty complicated, in order to maximize
!         memory locality
!
! #times abscissa is repeated : M/nabs(k)
! #repeats in a row           : key(k)
! #repeat cycles              : M/nabs(k)/key(k) 
! #length of repeat cycle     : nabs(k)*key(k) 
!
! Example
! nabs = 3,2
!  0  0
!  0  1  <- this is the first repeat cycle for 2nd dim
!  1  0
!  1  1  <- this is first repeat cycle for 1st dim, and 2nd for 2nd dim
!  2  0
!  2  1
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nabs          : 1D int, number of dimensions
! q             : 2D real*8, abscissa    [abscissa,dimension]
! W             : 2D real*8, weights     [weight,dimension]
! basK          : 1D real*8, basis frequencies
! Norm          : 1D real*8, normalization constants
! Heff          : 2D real*8, hermite poly
! PsiL          : 1D int, LHS quantum numbers 
! PsiR          : 1D int, RHS quantum numbers 
! Vq            : 1D int, potentials at abscissa
! Hij           : real*8, value to fill
! error         : int, error code

SUBROUTINE ints_HO_quadput(ndim,nabs,q,W,basK,Norm,Heff,PsiL,PsiR,Vq,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: q,W,Heff
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basK,Norm,Vq
  INTEGER, DIMENSION(0:), INTENT(IN) :: nabs,PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: Hij
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8), DIMENSION(0:PRODUCT(nabs)-1) :: val
  INTEGER, DIMENSION(0:ndim-1) :: key,ids
  INTEGER :: lo,hi,l0,h0
  INTEGER :: i,j,k,M
  error = 0
  M = PRODUCT(nabs)
  Hij = 0.0D0 

  CALL key_generate(ndim,nabs,key)
  val = Vq

  !Harmonic potential part
  DO k=0,ndim-1
    DO j=0,nabs(k)-1  !loop through unique abscissa
      l0 = j*key(k)
      h0 = l0 + key(k) - 1 
      DO i=0,(M/nabs(k))/key(k)-1 !loop through repeat cycles
        lo = i*key(k)*nabs(k)  + l0 
        hi = i*key(k)*nabs(k)  + h0 
        val(lo:hi) = val(lo:hi) - 0.5D0*basK(k)*q(j,k)**2.0D0
      END DO
    END DO
  END DO

  !Real Potential energy part
  DO k=0,ndim-1
    DO j=0,nabs(k)-1
      l0 = j*key(k)
      h0 = l0 + key(k) - 1 
      DO i=0,(M/nabs(k))/key(k)-1 !loop through repeat cycles
        lo = i*key(k)*nabs(k)  + l0 
        hi = i*key(k)*nabs(k)  + h0 
        val(lo:hi) = val(lo:hi)*Heff(j,2*k)*Heff(j,2*k+1)*&
                                     W(j,k)
      END DO
    END DO
  END DO

  !sum it up and normalize
  Hij = SUM(val(0:M-1))
  DO j=0,ndim-1
    Hij = Hij * Norm(2*j) 
    Hij = Hij * Norm(2*j+1)
  END DO

  !Kinetic energy part
  IF (ALL(PsiL .EQ. PsiR)) THEN
    DO k=0,ndim-1
      Hij = Hij + basK(k)*(1.0D0*PsiL(k)+0.5D0)
    END DO
  END IF

  CALL val_check(Hij,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "ints_HO_quadput  : ERROR"
    WRITE(*,*) "There is a bad value at this Matrix Element:"
    WRITE(*,*) "PsiL", PsiL
    WRITE(*,*) "PsiR", PsiR
  END IF

END SUBROUTINE ints_HO_quadput
!------------------------------------------------------------
! ints_HO_Q2eval
!       - evaluates quadratic force constant's contribution
!         to hamiltonian element in HO basis
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! qPhi          : 1D int, quadratic FC quantum numbers
! Phi           : real*8, quadratic FC value
! Q2int         : 2D real*8, <i|q^2|i'> type integrals 
! quadval       : real*8, value to add to 

SUBROUTINE ints_HO_quadeval(ndim,PsiL,PsiR,qPhi,Phi,&
                            Q2int,quadval)     
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q2int
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: quadval 
  REAL(KIND=8), INTENT(IN) :: Phi
  INTEGER, INTENT(IN) :: ndim,qPhi
  INTEGER :: i,j
  i = qPhi 
  !delta function for noninvolved dimensions
  IF (ALL(PsiL(0:i-1) .EQ. PsiR(0:i-1)) .AND.&
      ALL(PsiL(i+1:ndim-1) .EQ. PsiR(i+1:ndim-1)) ) THEN 
    ! q^2 can be v,v and v+2,v
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      j = PsiR(i)
      quadval = quadval + 0.5D0*phi*Q2int(2*j,i)
    ELSE IF (PsiL(i) .EQ. PsiR(i)+2) THEN
      j = PsiR(i)
      quadval = quadval + 0.5D0*phi*Q2int(2*j+1,i)
    END IF
  END IF
END SUBROUTINE ints_HO_quadeval

!------------------------------------------------------------
! ints_HO_momeval
!       - evaluates contribution of momentum
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! qPhi          : int, quantum number of FC
! Phi           : real*8, quadratic force constant
! P2int         : 2D real*8, <i|p^2|i'> type integrals
! momeval       : real*8, value 
SUBROUTINE ints_HO_momeval(ndim,PsiL,PsiR,qPhi,Phi,P2int,momeval)     
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
END SUBROUTINE ints_HO_momeval

!------------------------------------------------------------
! ints_HO_cubieval
!       - evaluates contribution of a cubic force constant
!         to the potential energy in harmonic oscillator 
!         basis
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

SUBROUTINE ints_HO_cubieval(diag,ndim,PsiL,PsiR,qPhi,&
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

  !check orthogonality of univolved dimensions
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
    WRITE(*,*) "ints_HO_cubieval  : ERROR"
    WRITE(*,*) "There is a cubic force constant type "
    WRITE(*,*) "  that has not been coded correctly"
  END IF

END SUBROUTINE ints_HO_cubieval

!------------------------------------------------------------
! ints_HO_quareval
!       - evaluates contribution of a quadratic force constant
!         to the potential energy in harmonic oscillator 
!         basis
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

SUBROUTINE ints_HO_quareval(diag,ndim,PsiL,PsiR,qPhi,&
                          Phi,Q1int,Q2int,Q3int,Q4int,quarval)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q1int,Q2int,Q3int,Q4int
  REAL(KIND=8), INTENT(INOUT) :: quarval
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,qPhi
  REAL(KIND=8), INTENT(IN) :: Phi
  LOGICAL, INTENT(IN) :: diag
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,j,k,l,m,n,o,p
  i = qPhi(0)
  j = qPhi(1)
  k = qPhi(2)
  l = qPhi(3)
  !Check orthogonality of noninvolved terms
  IF (ANY(PsiL(0:i-1) .NE. PsiR(0:i-1)) .OR. &
      ANY(PsiL(i+1:j-1) .NE. PsiR(i+1:j-1)) .OR. &
      ANY(PsiL(j+1:k-1) .NE. PsiR(j+1:k-1)) .OR. &
      ANY(PsiL(k+1:l-1) .NE. PsiR(k+1:l-1)) .OR. &
      ANY(PsiL(l+1:ndim-1) .NE. PsiR(l+1:ndim-1))) RETURN

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
    WRITE(*,*) "ints_HO_quareval  : ERROR"
    WRITE(*,*) "James, you missed a case!"
    WRITE(*,*) qPhi
    WRITE(*,*) PsiL,PsiR
  END IF
END SUBROUTINE ints_HO_quareval

!------------------------------------------------------------
END MODULE ints_HO
!------------------------------------------------------------
