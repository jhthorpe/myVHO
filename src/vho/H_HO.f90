!------------------------------------------------------------
! H_HO
!       - module containing subroutines dealing with the 
!         Hamiltonian in harmonic oscillator basis
!------------------------------------------------------------
MODULE H_HO
  USE V
  USE fcon
  USE basis
  USE ints_HO
  USE memory
  USE evec

CONTAINS
!------------------------------------------------------------
! H_HO_build
!       - builds the Hamiltonian in the harmonic osciallator 
!         basis
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! nabs          : int, number of abscissa
! q             : 1D real*8, list of abscissa
! W             : 1D real*8, list of weights
! Hij           : 2D real*8, hamiltonian
! error         : int, error code

SUBROUTINE H_HO_build(job,ndim,nbas,nabs,mem,q,W,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: q,W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(IN) :: job,ndim,nabs
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Vij, Herm
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Q1,Q2,Q3,Q4,P1,P2,QP,PQ,basK
  INTEGER, DIMENSION(:), ALLOCATABLE :: qQ1,qQ2,qQ3,qQ4,qP1,qP2,qQP,qPQ
  REAL(KIND=8) :: qmem
  REAL(KIND=8) :: ti,tf
  INTEGER :: nQ1,nQ2,nQ3,nQ4,nP1,nP2,nQP,nPQ
  INTEGER :: memstat
  INTEGER :: i,j,N

  CALL CPU_TIME(ti)
  error = 0
  N = PRODUCT(nbas)
  WRITE(*,*) "Beginning Construction of HO Hamiltonian"
  WRITE(*,*) "Dimension of Hamiltonian",N
  WRITE(*,*)

  !analyze the memory situation 
  CALL memory_HObuild(job,mem,N,ndim,nbas,nabs,memstat,error) 
  IF (error .NE. 0) RETURN
  IF (memstat .NE. 2) THEN 
    WRITE(*,*) "Sorry, only incore memory coded"
    error = 1
    RETURN
  END IF

  !If we're calculating V with gauss quad
  IF (job .EQ. 2) THEN
    !Read in basis set information
    CALL basis_get(ndim,basK,error)
    IF (error .NE. 0) RETURN

    !Get potential energies at abscissa 
    ALLOCATE(Vij(0:nabs-1,0:ndim-1))
    CALL V_get(job,ndim,nabs,q,Vij,error)
    IF (error .NE. 0) RETURN

    !Generate Hermite Polynomials
    IF (memstat .EQ. 2) THEN
      ALLOCATE(Herm(0:MAXVAL(nbas)-1,0:nabs-1))
      CALL H_HO_Hermite_incore(nabs,nbas,q,Herm,error)
    END IF
    IF (error .NE. 0) RETURN
  END IF

  !Read in force constant data 
  !CALL fcon_get(ndim,nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
  !                 nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,nQP,qQP,QP,&
  !                 nPQ,qPQ,PQ,error)
  CALL fcon_get(job,ndim,nQ2,qQ2,Q2,nQ3,qQ3,Q3,nQ4,qQ4,Q4,error)
  IF (error .NE. 0) RETURN

  IF (job .EQ. 1) THEN
    CALL H_HO_build_poly_incore(ndim,nbas,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                                nQ4,qQ4,Q4,Hij,error)
  END IF

  STOP

  !Evaluate integrals
  IF (job .NE. 2 .AND. job .NE. 1 ) THEN 
    WRITE(*,*) "H_HO_build  : ERROR"
    WRITE(*,*) "Sorry, only jobtype 0,1 are supported" 
    error = 1
    RETURN
  END IF
  IF (memstat .EQ. 2) THEN
    ALLOCATE(Hij(0:N-1,0:N-1))
    CALL H_HO_build_incore(ndim,nbas,nabs,q,W,Vij,basK,&
                    nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                    nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,&
                    nQP,qQP,QP,nPQ,qPQ,PQ,Herm,Hij,error)
  END IF

  IF (ALLOCATED(Vij)) DEALLOCATE(Vij)
  IF (ALLOCATED(Herm)) DEALLOCATE(Herm)

  CALL CPU_TIME(tf)
  WRITE(*,'(1x,A29,I2,A4,F6.1,1x,A1)') "H_HO_build  : finished with status ", error," in ", tf-ti, "s"
  WRITE(*,*)

END SUBROUTINE H_HO_build

!------------------------------------------------------------
! H_HO_Hermite_incore
!       - preconstruct all the Hermite polynomials needed
!------------------------------------------------------------
! nabs          : int, nubmer of abscissa
! nbas          : 1D int, nubmer of basis functions
! q             : 1D real*8, abscissa
! Herm          : 2D real*8, hermite poly [order,abscissa]
! error         : int, exit code

SUBROUTINE H_HO_Hermite_incore(nabs,nbas,q,Herm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Herm
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nabs
  
  REAL(KIND=8) :: qval
  INTEGER :: i,j,imax

  imax = MAXVAL(nbas)
  
  Herm = 0.0D0
  DO j=0,nabs-1
    qval = q(j)
    IF (imax .EQ. 1) THEN
      Herm(0,j) = 1.0D0
      CYCLE
    ELSE IF (imax .EQ. 2) THEN
      Herm(0,j) = 1.0D0
      Herm(1,j) = 2.0*qval
    ELSE
      Herm(0,j) = 1.0D0
      Herm(1,j) = 2.0*qval
      DO i=1,imax-2
        Herm(i+1,j) = 2.0D0*qval*Herm(i,j) - 2.0D0*i*Herm(i-1,j)
      END DO
    END IF
  END DO

END SUBROUTINE H_HO_Hermite_incore

!------------------------------------------------------------
! H_HO_build_incore
!       - HO basis
!       - constructs hamiltonian holding everything in memory
!       - naming scheme is as follows: aXY
!       - a = n,q,nothing. n -> number, q -> list of quant.
!             number, nothing -> actual values
!       - X = Q,P,QP,PQ
!       - Y = 1,2,3,4, order. 
!       - example: nQ4 -> number of Q^4 type force constants 
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, number of basis functions
! nabs          : int, number of abscissa
! q             : 1D real*8, abscissa [abscissa,dimension]
! W             : 2D real*8, weights  [weight,dimension]
! Vij           : 2D real*8, potential [potential,dimension]
! basK          : 1D real*8, basis function force constant
! nXY           : int, force constant X of order Y 
! qXY           : 1D int,  QN of FC X of order Y
! XY            : 1D real*8, FC X of order Y
! Herm          : 2D real*8, hermite polynomials [order,abscissa]
! Hij           : 2D real*8, product basis hamiltonian
! error         : int, exit code
!------------------------------------------------------------
SUBROUTINE H_HO_build_incore(ndim,nbas,nabs,q,W,Vij,basK,&
                     nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                     nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,&
                     nQP,qQP,QP,nPQ,qPQ,PQ,Herm,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Herm,Vij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q,W,Q1,Q2,Q3,Q4,P1,&
                                             P2,QP,PQ,basK 
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,qQ1,qQ2,qQ3,&
                                        qQ4,qP1,qP2,qQP,qPQ
  INTEGER, INTENT(IN) :: ndim,nabs,nQ1,nQ2,nQ3,nQ4,nP1,nP2,nQP,nPQ
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: VTint
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: norm
  REAL(KIND=8), DIMENSION(0:ndim-1) :: Vk
  REAL(KIND=8), DIMENSION(0:nabs-1) :: HL,HR
  INTEGER, DIMENSION(0:ndim-1) :: PsiL,PsiR,key
  INTEGER :: i,j,k,N
  
  error = 0
  Hij = 0.0D0
  N = PRODUCT(nbas)

  !generate keys
  CALL ints_HO_key(ndim,nbas,key,error)
  IF (error .NE. 0) RETURN

  !Generate normalization constants
  WRITE(*,*) "Calculating Normalization Constants..."
  ALLOCATE(norm(0:MAXVAL(nbas)-1))
  CALL ints_HO_normcalc(nabs,nbas,W,Herm,norm,error)

  !Generate V integrals
  WRITE(*,*) "Calculating Potential Integrals..."
  ALLOCATE(VTint(0:MAXVAL(nbas)-1,0:MAXVAL(nbas)-1,0:ndim-1))
  CALL ints_HO_VTcalc(ndim,nabs,nbas,q,W,basK,norm,Herm,Vij,VTint,error)

  !Generate Q1 integrals
  !WRITE(*,*) "Calculating Q1 integrals"
  
  !Generate Q2 integrals

  !Generate Q3 integrals

  !generate Q4 integrals

  !Perhaps this can be reversed? Loop over dimensions, and then
  ! over i,j ?
  !Go through vector by vector
  WRITE(*,*) "Filling in Hamiltonian..."
  DO j=0,N-1

    !generate quantum number of R
    CALL ints_HO_qnum(ndim,j,nbas,key,PsiR,error)
    IF (error .NE. 0) RETURN

    !lower triangular
    DO i=j,N-1
      CALL ints_HO_qnum(ndim,i,nbas,key,PsiL,error)
      IF (error .NE. 0) RETURN
     
      !Construct each part of the integral
      !Potential + Kinetic precalculated
      CALL ints_HO_VTput(ndim,PsiL,PsiR,VTint,Hij(i,j),error)
   
      !force constants
      !CALL ints_Q1put(nabs,PsiL,PsiR,Q1int,Hij(i,j),error)

      !check the value
      CALL val_check(Hij(i,j),error) 
      IF (error .NE. 0) THEN
        WRITE(*,*) "H_HO_build_incore  : ERROR"
        WRITE(*,*) "This integral had a bad value",N
        WRITE(*,*) "PsiL",PsiL
        WRITE(*,*) "PsiR",PsiR
      END IF
  
    END DO
  
  END DO

  WRITE(*,*)
  DEALLOCATE(norm)
  DEALLOCATE(VTint)

END SUBROUTINE H_HO_build_incore

!------------------------------------------------------------
! H_HO_build_poly_incore
!       - build hamiltonian where the potential is 
!         approximated as a 4th order polynomial in the 
!         harmonic oscillator basis set
!       - all terms are held incore
!       
!       PsiL and PsiR are the quantum nubmers of the left/right
!            wavefunctions
!
!       Integrals are stored like:
!       Q1(i,k) -> <i+1|Qk|i>
!       Q2(i,k) -> <i|Qk^2|i>    , Q2(i+1,k) -> <i+2|Qk^2|i>
!       Q3(i,k) -> <i+1|Qk^3|i>  , Q3(i+1,k) -> <i+3|Qk^3|i>
!       Q4(i,k) -> <i|Qk^4|i>    , Q4(i+1,k) -> <i+2|Qk^4|i>
!                                , Q4(i+2,k) -> <i+4|Qk^4|i>
!       P2(i,k) -> <i|Pk^2|i>    , P2(i,k)   -> <i+2|Pk^2|i>
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
! nbas          : 1D int, basis functions per dimension
! nQ2           : int, number of quadratic force constants
! qQ2           : 1D int, QN of quadratic force constants
! Q2            : 1D real*8, quadratic force constants
! nQ3           : int, number of cubic force constants
! qQ3           : 1D int, QN of cubic force constants
! Q3            : 1D real*8, cubic force constants
! nQ4           : int, number of quartic force constants
! qQ4           : 1D int, QN of quartic force constants
! Q4            : 1D real*8, quartic force constants
! error         : int, exit code

SUBROUTINE H_HO_build_poly_incore(ndim,nbas,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                                  nQ4,qQ4,Q4,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Q2,Q3,Q4 
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, DIMENSION(0:), INTENT(IN) :: qQ2,qQ3,qQ4
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nQ2,nQ3,nQ4
  INTEGER, INTENT(IN) :: ndim

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Q1int,Q2int,Q3int,Q4int,&
                                               P2int
  INTEGER, DIMENSION(0:ndim-1) :: PsiL,PsiR,key
  INTEGER :: N,mbas
  INTEGER :: i,j,k

  error = 0
  N = PRODUCT(nbas)
  mbas = MAXVAL(nbas)
  Hij = 0.0D0

  ALLOCATE(Q1int(0:mbas-1,0:ndim-1))
  ALLOCATE(Q2int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(Q3int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(Q4int(0:3*mbas-1,0:ndim-1))
  ALLOCATE(P2int(0:2*mbas-1,0:ndim-1))

  !precalculate the needed integrals
  !we could do this on the fly, but this
  !is probably faster
  WRITE(*,*) "Calculating Integrals..."
  CALL ints_HO_Q1calc(ndim,nbas,Q1int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q2calc(ndim,nbas,Q2int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q3calc(ndim,nbas,Q3int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q4calc(ndim,nbas,Q4int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_P2calc(ndim,nbas,P2int,error)
  IF (error .NE. 0) RETURN 

  CALL ints_HO_key(ndim,nbas,key,error)

  DO j=0,N-1
    CALL ints_HO_qnum(ndim,j,nbas,key,PsiR,error)
    DO i=j,N-1
      CALL ints_HO_qnum(ndim,j,nbas,key,PsiL,error)
      CALL ints_HO_polyput(ndim,PsiL,PsiR,nQ2,qQ2,Q2,&
                           nQ3,qQ3,Q3,nQ4,qQ4,Q4,Q1int,Q2int,&
                           Q3int,Q4int,P2int,Hij(i,j),error)       
    END DO 
  END DO 

  DEALLOCATE(Q1int)
  DEALLOCATE(Q2int)
  DEALLOCATE(Q3int)
  DEALLOCATE(Q4int)

END SUBROUTINE H_HO_build_poly_incore
!------------------------------------------------------------
! H_HO_diag
!       - diagonalizes the hamiltonian in the HO basis
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! enum          : int, number of eigenvalues to calc
! mem           : int*8, memory in MB
! Hij           : 2D real*8, hamiltonian
! eval          : 1D real*8, array of eigenvalues
! Cij           : 2D real*8, coefficients 
! error         : int, error 

SUBROUTINE H_HO_diag(ndim,nbas,enum,mem,Hij,eval,Cij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Cij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eval
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  INTEGER, DIMENSION(:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,enum

  INTEGER :: memstat,lwork,N
  INTEGER :: i
  
  error = 0
  N = PRODUCT(nbas)

  CALL memory_Hdiag(ndim,nbas,enum,mem,memstat,lwork,error)
  IF (error .NE. 0) RETURN
  IF (memstat .NE. 2) THEN
    WRITE(*,*) "H_diag  : ERROR"
    WRITE(*,*) "Sorry, only incore diagonalization has been coded"
    error = 1
    RETURN
  END IF

  ALLOCATE(eval(0:enum-1))
  ALLOCATE(Cij(0:N-1,0:enum-1))
 
  WRITE(*,*) "Diagonalizing the Hamiltonian"
  CALL linal_dsyevx(N,enum,lwork,Hij,eval,Cij,error)
  CALL evec_print(ndim,nbas,N,enum,eval,Cij,error)
  
END SUBROUTINE H_HO_diag

!------------------------------------------------------------
END MODULE H_HO
!------------------------------------------------------------
