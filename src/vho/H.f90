!------------------------------------------------------------
! H
!       - module containing subroutines dealing with the 
!         Hamiltonian
!------------------------------------------------------------
MODULE H
  USE V

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
! Hij             : 2D real*8, hamiltonian
! Herm          : 2D real*8, Hermite polynomials [dimension,abscissa]
! error         : int, error code

SUBROUTINE H_build(job,ndim,nbas,nabs,mem,q,W,Hij,Herm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Hij,Herm
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q,W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(IN) :: job,ndim,nabs
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Vij
  REAL(KIND=8) :: qmem
  REAL(KIND=8) :: ti,tf
  INTEGER :: memstat
  INTEGER :: i,j,N

  CALL CPU_TIME(ti)
  error = 0
  N = PRODUCT(nbas)
  WRITE(*,*) "H_build  : called"

  !analyze the memory situation 
  qmem = mem*1000000/8 !memory in qwords, 1 real*8 per qword
  CALL H_mem_build(job,qmem,N,ndim,nbas,nabs,memstat,error) 
  IF (error .NE. 0) RETURN

  IF (memstat .NE. 2) THEN 
    WRITE(*,*) "Sorry, only incore memory coded"
    error = 1
    RETURN
  END IF

  !Get potential energies at abscissa 
  IF (job .EQ. 0 .OR. job .EQ. 1) THEN
    ALLOCATE(Vij(0:ndim-1,0:nabs-1))
    CALL V_get(job,ndim,nabs,q,Vij,error)
    IF (error .NE. 0) RETURN
  END IF

  !Generate Hermite polynomials
  IF (memstat .EQ. 2) THEN
    ALLOCATE(Herm(0:MAXVAL(nbas)-1,0:nabs-1))
    CALL H_Herm_incore(job,ndim,nbas,nabs,q,Herm,error)
    IF (error .NE. 0) RETURN
  END IF

  !Generate coupling information


  !Evaluate integrals
  IF (job .NE. 0) THEN 
    WRITE(*,*) "Sorry, only jobtype 1 is supported"
    error = 1
    RETURN
  ELSE IF (job .EQ. 1) THEN

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
! nabs          : int, number of abscissa
! memstat       : int, memory status
! error         : int, exit code

SUBROUTINE H_mem_build(job,qmem,N,ndim,nbas,nabs,memstat,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  REAL(KIND=8), INTENT(IN) :: qmem
  INTEGER, INTENT(INOUT) :: memstat,error
  INTEGER, INTENT(IN) :: job,N,nabs,ndim
  REAL(KIND=8) :: minmem,incoremem,basemem
  
  error = 0
  WRITE(*,*) "H_mem_build  : called"
  WRITE(*,*) "Dimension of Hamiltonian", N
  WRITE(*,*) "Starting memory analysis" 
  WRITE(*,*) 

  basemem = 100
  IF (job .EQ. 0 .OR. job .EQ. 1) THEN
    ! hamiltonian : N^2
    ! abscissa    : nabs
    ! weights     : nabs
    ! V           : nabs * ndim
    ! hermite     : either nabs*max(nbas) or max(nbas)
    minmem = N**2.0D0 + 2*nabs + MAXVAL(nbas) + nabs*ndim + basemem
    incoremem = N**2.0D0 + 2*nabs + MAXVAL(nbas)*nabs + nabs*ndim + basemem
  END IF

  WRITE(*,*) "Values are in MB"
  WRITE(*,'(A20,F12.2)') "Available memory   ", qmem*8.0D0/1000000.0D0
  WRITE(*,'(A20,F12.2)') "Minimum memory     ", minmem*8.0D0/1000000.0D0
  WRITE(*,'(A20,F12.2)') "Incore memory      ", incoremem*8.0D0/1000000.0D0

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
! ndim          : int, nubmer of dimensions
! nabs          : int, nubmer of abscissa
! q             : 1D real*8, abscissa
! error         : int, exit code

SUBROUTINE H_Herm_incore(job,ndim,nbas,nabs,q,Herm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Herm
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,nabs
  
  REAL(KIND=8) :: qval
  INTEGER :: i,j,imax

  imax = MAXVAL(nbas)-1
  
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
      DO i=1,imax-1
        Herm(i+1,j) = 2.0D0*qval*Herm(i,j) - 2.0D0*i*Herm(i-1,j)
      END DO
    END IF
  END DO

END SUBROUTINE H_Herm_incore
!------------------------------------------------------------

END MODULE H
!------------------------------------------------------------
