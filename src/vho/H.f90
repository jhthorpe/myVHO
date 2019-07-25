!------------------------------------------------------------
! H
!       - module containing subroutines dealing with the 
!         Hamiltonian
!------------------------------------------------------------
MODULE H
  USE V
  USE fcon
  USE basis
  USE ints
  USE memory
  USE evec

CONTAINS
!------------------------------------------------------------
! H_build
!       - builds the Hamiltonian
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

SUBROUTINE H_build(job,ndim,nbas,nabs,mem,q,W,Hij,error)
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
  WRITE(*,*) "Beginning Construction of Hamiltonian"
  WRITE(*,*) "Dimension of Hamiltonian",N
  WRITE(*,*)

  !analyze the memory situation 
  qmem = mem*1000000/8 !memory in qwords, 1 real*8 per qword
  CALL memory_Hbuild(job,mem,N,ndim,nbas,nabs,memstat,error) 
  IF (error .NE. 0) RETURN
  IF (memstat .NE. 2) THEN 
    WRITE(*,*) "Sorry, only incore memory coded"
    error = 1
    RETURN
  END IF

  !Read in basis set information
  CALL basis_get(ndim,basK,error)
  IF (error .NE. 0) RETURN

  !Get potential energies at abscissa 
  IF (job .EQ. 0 .OR. job .EQ. 1) THEN
    ALLOCATE(Vij(0:nabs-1,0:ndim-1))
    CALL V_get(job,ndim,nabs,q,Vij,error)
    IF (error .NE. 0) RETURN
  END IF

  !Read in force constant data 
  CALL fcon_get(ndim,nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                   nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,nQP,qQP,QP,&
                   nPQ,qPQ,PQ,error)
  IF (error .NE. 0) RETURN

  !Generate Hermite Polynomials
  IF (memstat .EQ. 2) THEN
    ALLOCATE(Herm(0:MAXVAL(nbas)-1,0:nabs-1))
    CALL H_Hermite_incore(nabs,nbas,q,Herm,error)
  END IF
  IF (error .NE. 0) RETURN

  !Evaluate integrals
  IF (job .NE. 1 .AND. job .NE. 0 ) THEN 
    WRITE(*,*) "H_build  : ERROR"
    WRITE(*,*) "Sorry, only jobtype 0,1 are supported" 
    error = 1
    RETURN
  END IF
  IF (memstat .EQ. 2) THEN
    ALLOCATE(Hij(0:N-1,0:N-1))
    CALL H_build_incore(ndim,nbas,nabs,q,W,Vij,basK,&
                    nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                    nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,&
                    nQP,qQP,QP,nPQ,qPQ,PQ,Herm,Hij,error)
  END IF

  DEALLOCATE(Vij)
  DEALLOCATE(Herm)

  CALL CPU_TIME(tf)
  WRITE(*,'(1x,A29,I2,A4,F6.1,1x,A1)') "H_build  : finished with status ", error," in ", tf-ti, "s"
  WRITE(*,*)

END SUBROUTINE H_build

!------------------------------------------------------------
! H_Hermite_incore
!       - preconstruct all the Hermite polynomials needed
!------------------------------------------------------------
! nabs          : int, nubmer of abscissa
! nbas          : 1D int, nubmer of basis functions
! q             : 1D real*8, abscissa
! Herm          : 2D real*8, hermite poly [order,abscissa]
! error         : int, exit code

SUBROUTINE H_Hermite_incore(nabs,nbas,q,Herm,error)
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

END SUBROUTINE H_Hermite_incore

!------------------------------------------------------------
! H_build_incore
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
SUBROUTINE H_build_incore(ndim,nbas,nabs,q,W,Vij,basK,&
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
  CALL ints_key(ndim,nbas,key,error)
  IF (error .NE. 0) RETURN

  !Generate normalization constants
  WRITE(*,*) "Calculating Normalization Constants..."
  ALLOCATE(norm(0:MAXVAL(nbas)-1))
  CALL ints_normcalc(nabs,nbas,W,Herm,norm,error)

  !Generate V integrals
  WRITE(*,*) "Calculating Potential Integrals..."
  ALLOCATE(VTint(0:MAXVAL(nbas)-1,0:MAXVAL(nbas)-1,0:ndim-1))
  CALL ints_VTcalc(ndim,nabs,nbas,q,W,basK,norm,Herm,Vij,VTint,error)

  !Perhaps this can be reversed? Loop over dimensions, and then
  ! over i,j ?
  !Go through vector by vector
  WRITE(*,*) "Filling in Hamiltonian..."
  DO j=0,N-1

    !generate quantum number of R
    CALL ints_qnum(ndim,j,nbas,key,PsiR,error)
    IF (error .NE. 0) RETURN

    !lower triangular
    DO i=j,N-1
      CALL ints_qnum(ndim,i,nbas,key,PsiL,error)
      IF (error .NE. 0) RETURN
     
      !Construct each part of the integral
      !Potential + Kinetic precalculated
      CALL ints_VTput(ndim,PsiL,PsiR,VTint,Hij(i,j),error)
   
      !force constants
      !CALL ints_Q1put(nabs,PsiL,PsiR,Q1int,Hij(i,j),error)

      !check the value
      CALL val_check(Hij(i,j),error) 
      IF (error .NE. 0) THEN
        WRITE(*,*) "H_build_incore  : ERROR"
        WRITE(*,*) "This integral had a bad value",N
        WRITE(*,*) "PsiL",PsiL
        WRITE(*,*) "PsiR",PsiR
      END IF
  
    END DO
  
  END DO

  WRITE(*,*)
  DEALLOCATE(norm)
  DEALLOCATE(VTint)

END SUBROUTINE H_build_incore

!------------------------------------------------------------
! H_diag
!       - diagonalizes the hamiltonian
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! enum          : int, number of eigenvalues to calc
! mem           : int*8, memory in MB
! Hij           : 2D real*8, hamiltonian
! eval          : 1D real*8, array of eigenvalues
! Cij           : 2D real*8, coefficients 
! error         : int, error 

SUBROUTINE H_diag(ndim,nbas,enum,mem,Hij,eval,Cij,error)
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
  
END SUBROUTINE H_diag

!------------------------------------------------------------
END MODULE H
!------------------------------------------------------------
