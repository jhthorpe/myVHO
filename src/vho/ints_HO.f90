!------------------------------------------------------------
! ints_HO
!       - module supporting integral calculations in the
!         HO basis
!------------------------------------------------------------
MODULE ints_HO
  USE val

CONTAINS
!------------------------------------------------------------
! ints_HO_key
!       - generate the proper key for determining quantum
!         numbers
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, nubmer of basis functions per d
! key           : 1D int, key
! error         : int, exit code

SUBROUTINE ints_HO_key(ndim,nbas,key,error)
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
END SUBROUTINE ints_HO_key

!------------------------------------------------------------
! ints_HO_qnum
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

SUBROUTINE ints_HO_qnum(ndim,it,nbas,key,qnum,error)
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
END SUBROUTINE ints_HO_qnum

!------------------------------------------------------------
! ints_HO_Norm
!       - calculates normalization integrals
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
!------------------------------------------------------------
! nabs          : int, nubmer of abscissa
! nbas          : 1D int, number of basis functions
! W             : 1D real*8, weights
! Herm          : 2D real*8, hermite polynomials
! norm          : 1D real*8, normalization constants
! error         : int, exit code

SUBROUTINE ints_HO_normcalc(nabs,nbas,W,Herm,norm,error)
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
    CALL ints_HO_norm(nabs,W,HR,norm(i),error)
    IF (error .NE. 0) RETURN
    CALL val_check(norm(i),error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "ints_HO_normcalc  : ERROR"
      WRITE(*,*) "This normalization constant had a bad value", i
      RETURN
    END IF
  END DO
END SUBROUTINE ints_HO_normcalc
!------------------------------------------------------------
! ints_HO_VTcalc
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

SUBROUTINE ints_HO_VTcalc(ndim,nabs,nbas,q,W,basK,norm,Herm,Vij,&
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
        CALL ints_HO_VTint(nabs,q,W,HL,HR,basK(k),Vij(:,k),VTint(i,j,k),error)
        IF (error .NE. 0) RETURN
        VTint(i,j,k) = VTint(i,j,k)*norm(i)*norm(j)

        ! + HO value (this is the kinetic term)  
        IF (i .EQ. j) VTint(i,j,k) = VTint(i,j,k) &
                      + basK(k)*(1.0D0*i+0.5D0)

        CALL val_check(VTint(i,j,k),error)
        IF (error .NE. 0) THEN
          WRITE(*,*) "ints_HO_VTcalc  : ERROR"
          WRITE(*,*) "Bad potential at i,j,k",i,j,k
          RETURN
        END IF

      END DO
    END DO
  END DO

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
      Q1int(i,j) = SQRT(0.5D0*(i+1.0D0))
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
  INTEGER :: i,j
  error = 0
  quadval = 0.0D0
  cubival = 0.0D0
  quarval = 0.0D0
  momeval = 0.0D0
  !phi_ij (quadratic) terms
  DO i=0,nQ2-1
    CALL ints_HO_quadeval(ndim,PsiL,PsiR,qQ2(i),Q2(i),&
                        Q2int,quadval)     
  END DO

  !p^2 (momentum) terms
  DO i=0,nQ2-1
    CALL ints_HO_momeeval(ndim,PsiL,PsiR,qQ2(i),Q2(i),&
                        P2int,momeval)     
  END DO
  
  Hij = quadval + cubival + quarval + momeval

END SUBROUTINE ints_HO_polyput

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
! ints_HO_momeeval
!       - evaluates contribution of momentum
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! qPhi          : int, quantum number of FC
! Phi           : real*8, quadratic force constant
! P2int         : 2D real*8, <i|p^2|i'> type integrals
! momeval       : real*8, value 
SUBROUTINE ints_HO_momeeval(ndim,PsiL,PsiR,qPhi,Phi,P2int,momeval)     
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
END SUBROUTINE ints_HO_momeeval
!------------------------------------------------------------
END MODULE ints_HO
!------------------------------------------------------------
