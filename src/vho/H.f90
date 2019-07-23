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

CONTAINS
!------------------------------------------------------------
! H_build
!       - builds the Hamiltonian
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! nabs          : 1D int, number of abscissa
! q             : 1D real*8, list of abscissa
! W             : 1D real*8, list of weights
! Hij             : 2D real*8, hamiltonian
! Herm          : 2D real*8, Hermite polynomials [quantum number,abscissa]
! error         : int, error code

SUBROUTINE H_build(job,ndim,nbas,nabs,mem,q,W,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: q,W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(IN) :: job,ndim
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Vij
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

  !analyze the memory situation 
  qmem = mem*1000000/8 !memory in qwords, 1 real*8 per qword
  CALL H_mem_build(job,qmem,N,ndim,nbas,nabs,memstat,error) 
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
    ALLOCATE(Vij(0:MAXVAL(nabs)-1,0:ndim-1))
    CALL V_get(job,ndim,nabs,q,Vij,error)
    IF (error .NE. 0) RETURN
  END IF

  !Read in force constant data 
  CALL fcon_get(ndim,nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                   nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,nQP,qQP,QP,&
                   nPQ,qPQ,PQ,error)
  IF (error .NE. 0) RETURN

  !Evaluate integrals
  IF (job .NE. 1 ) THEN 
    WRITE(*,*) "Sorry, only jobtype 1 is supported"
    error = 1
    RETURN
  ELSE IF (job .EQ. 1) THEN
    IF (memstat .EQ. 2) THEN
      CALL H_build_incore(ndim,nbas,nabs,q,W,Vij,basK,&
                      nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                      nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,&
                      nQP,qQP,QP,nPQ,qPQ,PQ,Hij,error)
    END IF
  END IF

  CALL CPU_TIME(tf)
  WRITE(*,'(1x,A29,I2,A4,F6.1,1x,A1)') "H_build  : finished with status ", error," in ", tf-ti, "s"

END SUBROUTINE H_build

!------------------------------------------------------------
! H_mem_build
!       - analyses the memory situation for building the 
!         hamiltonian
!------------------------------------------------------------
! job           : int, jobtype
! qmem          : int, memory in qwords (1 qword = 1 DP)
! N             : int, size of hamiltonian
! ndim          : int, number of dimensions
! nbas          : int, number of basis in each dimension
! nabs          : 1D int, number of abscissa
! memstat       : int, memory status
! error         : int, exit code

SUBROUTINE H_mem_build(job,qmem,N,ndim,nbas,nabs,memstat,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  REAL(KIND=8), INTENT(IN) :: qmem
  INTEGER, INTENT(INOUT) :: memstat,error
  INTEGER, INTENT(IN) :: job,N,ndim
  REAL(KIND=8) :: minmem,incoremem,basemem
  
  error = 0
  WRITE(*,*) "Dimension of Hamiltonian", N
  WRITE(*,*) 
  WRITE(*,*) "Starting memory analysis..." 

  basemem = 1000
  IF (job .EQ. 0 .OR. job .EQ. 1) THEN
    ! hamiltonian : N^2
    ! abscissa    : nabs
    ! weights     : nabs
    ! V           : nabs * ndim
    ! hermite     : either nabs*max(nbas) or max(nbas)
    minmem = N**2.0D0 + 2*MAXVAL(nabs) + MAXVAL(nbas) + MAXVAL(nabs)*ndim + basemem
    incoremem = N**2.0D0 + 2*MAXVAL(nabs) + MAXVAL(nbas)*MAXVAL(nabs) + MAXVAL(nabs)*ndim + basemem
  END IF

  WRITE(*,'(A20,F12.2)') "Available memory   ", qmem*8.0D0/1000000.0D0
  WRITE(*,'(A20,F12.2)') "Minimum memory     ", minmem*8.0D0/1000000.0D0
  WRITE(*,'(A20,F12.2)') "Incore memory      ", incoremem*8.0D0/1000000.0D0
  WRITE(*,*) "Values are in MB"

  WRITE(*,*)
  IF (qmem .LT. minmem) THEN
    WRITE(*,*) "H_mem_build  : ERROR"
    WRITE(*,*) "Not enough memory has been given"
    error = 2
    memstat = -1
    RETURN
  ELSE IF (qmem .LT. incoremem .AND. qmem .GE. minmem) THEN
    WRITE(*,*) "Building Hamiltonian with minimal memory"
    memstat = 0
    
  ELSE IF (qmem .GE. incoremem) THEN
    WRITE(*,*) "Building Hamiltonian with everything in core"
    memstat = 2
  ELSE
    WRITE(*,*) "H_mem_build  : ERROR" 
    WRITE(*,*) "Yours truely has somehow coded this badly..."
    error = 3
    RETURN
  END IF 
  WRITE(*,*) 

END SUBROUTINE H_mem_build

!------------------------------------------------------------
! H_Herm_incore
!       - preconstruct all the Hermite polynomials needed
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : 1D int, nubmer of dimensions
! nabs          : 1D int, nubmer of abscissa
! q             : 1D real*8, abscissa
! error         : int, exit code

SUBROUTINE H_Herm_incore(job,ndim,nbas,nabs,q,Herm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Herm
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: q
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim
  
  REAL(KIND=8) :: qval
  INTEGER :: i,j,imax

  imax = MAXVAL(nbas)-1
  
  Herm = 0.0D0
  !DO j=0,nabs-1
  DO j=0,-1  ! THIS IS VERY WRONG
    !qval = q(j)
    qval = 0 !!! THIS IS VERY WRONG
    IF (imax .EQ. 1) THEN
      Herm(0,j) = 1.0D0
      CYCLE
    ELSE IF (imax .EQ. 2) THEN
      Herm(0,j) = 1.0D0
      Herm(1,j) = 2.0*qval
    ELSE
      Herm(0,j) = 1.0D0
      Herm(1,j) = 2.0*qval
      DO i=1,imax-1
        Herm(i+1,j) = 2.0D0*qval*Herm(i,j) - 2.0D0*i*Herm(i-1,j)
      END DO
    END IF
  END DO

  WRITE(*,*) "WARNING WARNING WARNING"
  WRITE(*,*) "H_Herm_incore has not been adapted to the new setup"
  STOP

END SUBROUTINE H_Herm_incore

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
! nabs          : 1D int, number of abscissa
! q             : 2D real*8, abscissa [abscissa,dimension]
! W             : 2D real*8, weights  [weight,dimension]
! Vij           : 2D real*8, potential [potential,dimension]
! basK          : 1D real*8, basis function force constant
! nXY           : int, force constant X of order Y 
! qXY           : 1D int,  QN of FC X of order Y
! XY            : 1D real*8, FC X of order Y
! Hij           : 2D real*8, product basis hamiltonian
! error         : int, exit code
!------------------------------------------------------------
SUBROUTINE H_build_incore(ndim,nbas,nabs,q,W,Vij,basK,&
                     nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                     nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,&
                     nQP,qQP,QP,nPQ,qPQ,PQ,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: q,W,Vij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Q1,Q2,Q3,Q4,P1,&
                                             P2,QP,PQ,basK 
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs,qQ1,qQ2,qQ3,&
                                        qQ4,qP1,qP2,qQP,qPQ
  INTEGER, INTENT(IN) :: ndim,nQ1,nQ2,nQ3,nQ4,nP1,nP2,nQP,nPQ
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(0:MAXVAL(nabs)-1) :: HL,HR
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Norm,qmax,Wmax
  INTEGER, DIMENSION(0:ndim-1) :: PsiL,PsiR,key
  INTEGER :: i,j,N,nmax
  
  error = 0
  N = PRODUCT(nbas)
  nmax = MAXVAL(nabs)  
  j = MAXLOC(nabs,1)-1
  ALLOCATE(Norm(0:N-1))
  ALLOCATE(qmax(0:nmax-1))
  ALLOCATE(Wmax(0:nmax-1))
  qmax = q(0:nmax-1,j)
  Wmax = W(0:nmax-1,j)

  CALL ints_key(ndim,nbas,key,error)
  IF (error .NE. 0) RETURN

  !Go through vector by vector
  DO j=0,N-1

    !generate Hermite polynomials, and quantum numbers of R
    CALL ints_qnum(ndim,j,nbas,key,PsiR,error)
    IF (error .NE. 0) RETURN

    !CALL H_Hermite()
    IF (error .NE. 0) RETURN
    

  END DO

  DEALLOCATE(Norm)
  DEALLOCATE(qmax)
  DEALLOCATE(Wmax)

END SUBROUTINE H_build_incore

!------------------------------------------------------------
END MODULE H
!------------------------------------------------------------
