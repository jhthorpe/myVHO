!------------------------------------------------------------
! memory
!       - module for managing and analysing memory
!------------------------------------------------------------
MODULE memory
  USE val
  USE linal

CONTAINS
!------------------------------------------------------------
! memory_Hbuild 
!       - analyses the memory situation for building the 
!         hamiltonian
!------------------------------------------------------------
! job           : int, jobtype
! mem           : int, memory in MB
! N             : int, size of hamiltonian
! ndim          : int, number of dimensions
! nbas          : int, number of basis in each dimension
! nabs          : 1D int, number of abscissa
! memstat       : int, memory status
! error         : int, exit code

SUBROUTINE memory_Hbuild(job,mem,N,ndim,nbas,nabs,memstat,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: memstat,error
  INTEGER, INTENT(IN) :: job,N,ndim,nabs
  REAL(KIND=8) :: minmem,incoremem,basemem,qw2mb,qmem

  error = 0
  qw2mb = 8.0D0/1000000.0D0
  qmem = mem/qw2mb
  WRITE(*,*) "Starting Hbuild memory analysis..."
  CALL val_check(qmem,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "memory_Hbuild  : ERROR"
    WRITE(*,*) "qmem was too large to be analysed"
    error = 1
    RETURN
  END IF

  basemem = 1000
  IF (job .EQ. 0 .OR. job .EQ. 1) THEN
    ! hamiltonian : N^2
    ! abscissa    : nabs
    ! weights     : nabs
    ! V           : nabs * ndim
    ! hermite     : either nabs*max(nbas) or max(nbas)
    minmem = N**2.0D0 + 2*nabs + MAXVAL(nbas) + nabs*ndim + basemem
    incoremem = N**2.0D0 + 2*nabs + MAXVAL(nbas)*nabs + &
                nabs*ndim + basemem
  END IF

  WRITE(*,*) "JAMES, CHECK THIS IS STILL ALL GOOD"
  WRITE(*,'(A20,F12.2)') "Available memory   ", qmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Minimum memory     ", minmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Incore memory      ", incoremem*qw2mb
  WRITE(*,*) "Values are in MB"

  WRITE(*,*)
  IF (qmem .LT. minmem) THEN
    WRITE(*,*) "memory_Hbuild  : ERROR"
    WRITE(*,*) "Not enough memory has been given"
    error = 2
    memstat = -1
    RETURN
  ELSE IF (qmem .LT. incoremem .AND. qmem .GE. minmem) THEN
    WRITE(*,*) "Hamiltonian will be built with minimal memory"
    memstat = 0

  ELSE IF (qmem .GE. incoremem) THEN
    WRITE(*,*) "Hamiltonian will be built with incore memory"
    memstat = 2
  ELSE
    WRITE(*,*) "memory_Hbuild  : ERROR"
    WRITE(*,*) "Yours truely has somehow coded this badly..."
    error = 3
    RETURN
  END IF
  WRITE(*,*)

END SUBROUTINE memory_Hbuild

!------------------------------------------------------------
! memory_Hdiag
!       - figures out the memory situation for the 
!         diagonalization subroutine
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! enum          : int, number of eigenvalues to calc
! mem           : int*8, memory in MB
! memstat       : int, memory status, output
! lwork         : int, length of working vector allowed
! error         : int, error 

SUBROUTINE memory_Hdiag(ndim,nbas,enum,mem,memstat,lwork,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: memstat,error,lwork
  INTEGER, INTENT(IN) :: ndim,enum
  REAL(KIND=8) :: minmem,incoremem,basemem,qw2mb,qmem
  INTEGER :: N

  error = 0
  N = PRODUCT(nbas)
  qw2mb = 8.0D0/1000000.0D0

  qmem = mem/qw2mb
  CALL val_check(qmem,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "memory_Hdiag  : ERROR"
    WRITE(*,*) "qmem was too large to be analysed"
    error = 1
    RETURN
  END IF
 
  !Memory required
  ! in background...
  ! q           : nabs
  ! W           : nabs
  ! in diag routine...
  ! Hij         : N^2
  ! Cij         : N or N*enum 
  ! eval        : enum
  ! WORK        : 3*N, lwork 
  ! IWORK       : 5*N
  
  
  CALL linal_dsyevx_lwork(N,enum,lwork,error)
  IF (error .NE. 0) RETURN
  
  basemem = 1000
  minmem = basemem + N**2 + N + enum + 3*N + 5.0D0*N/2.0D0
  incoremem = basemem + N**2 +  N*enum + enum + lwork + 5.0D0*N/2.0D0
  
  WRITE(*,*) "Diagonalization Memory Analysis"
  WRITE(*,'(A20,F12.2)') "Available memory   ", qmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Minimum memory     ", minmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Incore memory      ", incoremem*qw2mb
  WRITE(*,*) "Values are in MB"

  WRITE(*,*)
  IF (qmem .LT. minmem) THEN
    WRITE(*,*) "memory_Hdiag  : ERROR"
    WRITE(*,*) "Not enough memory has been given"
    error = 2
    memstat = -1
    RETURN
  ELSE IF (qmem .LT. incoremem .AND. qmem .GE. minmem) THEN
    WRITE(*,*) "Hamiltonian will be built with minimal memory"
    memstat = 0

  ELSE IF (qmem .GE. incoremem) THEN
    WRITE(*,*) "Hamiltonian will be built with incore memory"
    memstat = 2
  ELSE
    WRITE(*,*) "memory_Hdiag  : ERROR"
    WRITE(*,*) "Yours truely has somehow coded this badly..."
    error = 3
    RETURN
  END IF
  WRITE(*,*)

END SUBROUTINE memory_Hdiag

!------------------------------------------------------------


END MODULE memory
