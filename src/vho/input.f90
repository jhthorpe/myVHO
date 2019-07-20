!------------------------------------------------------------
! input
!       - module for parsing input for xvho
!------------------------------------------------------------

MODULE input

CONTAINS

!------------------------------------------------------------
! input_jobinfo
!       - gets the job info
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! error         : int, error code

SUBROUTINE input_jobinfo(job,ndim,nbas,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: nbas
  INTEGER(KIND=8), INTENT(INOUT) :: mem
  INTEGER, INTENT(INOUT) :: job,ndim,error

  error = 0
  WRITE(*,*) "input_jobinfo : called"
  WRITE(*,*) "Reading input from vho.in"

  CALL input_read(job,ndim,nbas,mem,error)
  IF (error .NE. 0) RETURN
  CALL input_check(job,ndim,nbas,mem,error)
  IF (error .NE. 0) RETURN
  CALL input_write(job,ndim,nbas,mem,error)

  WRITE(*,*) "input_jobinfo : completed with status ", error
  WRITE(*,*) "======================================================="

END SUBROUTINE input_jobinfo

!------------------------------------------------------------
! input_read
!    - reads input from vho.in file
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB

SUBROUTINE input_read(job,ndim,nbas,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: nbas
  INTEGER(KIND=8), INTENT(INOUT) :: mem
  INTEGER, INTENT(INOUT) :: job,ndim, error
  LOGICAL :: ex
  error = 0
  INQUIRE(file='vho.in',EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You must create the input file, vho.in"
    error = 1
    RETURN
  END IF
  OPEN(file='vho.in',unit=100,status='old')
  READ(100,*) job
  READ(100,*) ndim
  ALLOCATE(nbas(0:ndim-1))
  READ(100,*) nbas 
  READ(100,*) mem
  CLOSE(unit=100)

END SUBROUTINE input_read

!------------------------------------------------------------
! input_check
!       - checks input of vho.in
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB

SUBROUTINE input_check(job,ndim,nbas,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:ndim-1), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim
  error = 0
  IF (job .NE. 0) THEN
    WRITE(*,*) "vho.in line #1"
    WRITE(*,*) "Jobtype", job," is not supported. Options are..."
    WRITE(*,*) "0 : generate abscissa and weights needed"
    error = 1
  END IF
  IF (ndim .LT. 1) THEN
    WRITE(*,*) "vho.in line #2"
    WRITE(*,*) "Must have at least one dimension"
    error = 2
  END IF
  IF (MINVAL(nbas) .LT. 1) THEN
    WRITE(*,*) "vho.in line #3"
    WRITE(*,*) "All dimensions need at least one basis function"
    error = 3
  END IF
  IF (mem .LT. 1) THEN
    WRITE(*,*) "vho line #4"
    WRITE(*,*) "Must have at least 1 MB of memory"
    error = 4
  END IF
END SUBROUTINE input_check

!------------------------------------------------------------
! input_write
!       - writes the input of vho.in
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! error         : int, error code

SUBROUTINE input_write(job,ndim,nbas,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:ndim-1), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim
  error = 0
  CALL EXECUTE_COMMAND_LINE('cat vho.in')
  IF (job .EQ. 0) THEN
    WRITE(*,*) "job     : 0 - abscissa and weights will be determined"
  END IF
  WRITE(*,*) "ndim    :",ndim
  WRITE(*,*) "nbas    :",nbas
  WRITE(*,*) "mem     :",mem 
END SUBROUTINE input_write
!------------------------------------------------------------

END MODULE input