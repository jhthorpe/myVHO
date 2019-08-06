!------------------------------------------------------------
! ints_HO
!       - module supporting integral calculations in the
!         HO basis
!------------------------------------------------------------
MODULE ints_HO
  USE val
  USE gauss
  USE key
  USE quad
  USE cubi
  USE quar
  USE mome

CONTAINS

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
    DO j=0,nbas(k)-1
      DO i=j,nbas(k)-1
        val = 0.0D0
        !potential - harmonic 
        DO a=0,nabs(k)-1
          val = val + W(a,k)*Herm(a,i,k)*Herm(a,j,k)*(Vij(a,k) &
                - (0.5D0*basK(k)*q(a,k)**2.0D0))
        END DO
        val = val*norm(i,k)*norm(j,k)
      
        !kinetic energy part
        IF (i .EQ. j) val = val + bask(k)*(1.0D0*i+0.5D0)
        VTint(i-j+keyI(j,k),k) = val

        CALL val_check(val,error)
        IF (error .NE. 0) THEN
          WRITE(*,*) "ints_HO_VTcalc  : ERROR"
          WRITE(*,*) "Bad VT integral at i,j,k",i,j,k
          RETURN
        END IF 
      END DO
    END DO
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
    CALL quad_HO_eval(ndim,PsiL,PsiR,qQ2(i),Q2(i),&
                        Q2int,quadval)     
  END DO

  !phi_ijk (cubic terms)
  DO i=0,nQ3-1
    CALL cubi_HO_eval(diag,ndim,PsiL,PsiR,qQ3(3*i:3*i+2),&
                          Q3(i),Q1int,Q2int,Q3int,cubival)
  END DO

  !phi_ijkl (quartic terms)
  DO i=0,nQ4-1
    CALL quar_HO_eval(diag,ndim,PsiL,PsiR,qQ4(4*i:4*i+3),&
                          Q4(i),Q1int,Q2int,Q3int,Q4int,quarval)
  END DO

  !p^2 (momentum) terms
  DO i=0,nQ2-1
    CALL mome_HO_eval(ndim,PsiL,PsiR,qQ2(i),Q2(i),&
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
    CALL cubi_HO_eval(diag,ndim,PsiL,PsiR,qQ3(3*i:3*i+2),&
                          Q3(i),Q1int,Q2int,Q3int,cubival)
  END DO

  !phi_ijkl (quartic terms)
  DO i=0,nQ4-1
    CALL quar_HO_eval(diag,ndim,PsiL,PsiR,qQ4(4*i:4*i+3),&
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
!  1  1  
!  2  0
!  2  1 <- this is the first repeat cycle for the 1st dim
!
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
  REAL(KIND=8) :: temp
  INTEGER :: lo,hi,l0,h0
  INTEGER :: i,j,k,M
  error = 0
  M = PRODUCT(nabs)
  Hij = 0.0D0 

  !IF (.FALSE.) THEN
  !slow code for testing
  !which I think is probably wrong...
  !CALL key_generate(ndim,nabs,key)
  !DO k=0,M-1
  !  CALL key_idx2ids(ndim,k,nabs,key,ids)
  !  temp = Vq(k) 
  !  DO j=0,ndim-1
  !    temp = temp - 0.5*bask(k)*q(ids(j),j)**2.0D0
  !  END DO 
  !  DO j=0,ndim-1
  !    temp = temp*W(ids(j),j)*Heff(ids(j),2*j)*Heff(ids(j),2*j+1)
  !  END DO
  !  Hij = Hij + temp 
  !END DO

  !ELSE
  
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
  !END IF
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
END MODULE ints_HO
!------------------------------------------------------------
