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
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim
  error = 0
  IF (job .LT. -1 .OR. job .GT. 1) THEN
    WRITE(*,*) "vho.in line #1"
    WRITE(*,*) "Jobtype", job," is not supported. Options are..."
    WRITE(*,*) "-1 : print abscissa and weights needed"
    WRITE(*,*) " 0 : use precalculated abscissa and weights"
    WRITE(*,*) " 1 : cublic spline interpolation of potentials" 
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
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim
  error = 0
  CALL EXECUTE_COMMAND_LINE('cat vho.in')
  IF (job .EQ. -1) THEN
    WRITE(*,*) "job     : -1 - abscissa and weights will be printed"
  ELSE IF (job .EQ. 0) THEN
    WRITE(*,*) "job     : 0 - precalculated abscissa will be used"
  ELSE IF (job .EQ. 1) THEN 
    WRITE(*,*) "job     : 1 - potential interpolated by cubic spline" 
  END IF
  WRITE(*,*) "ndim    :",ndim
  WRITE(*,*) "nbas    :",nbas
  WRITE(*,*) "mem     :",mem 
END SUBROUTINE input_write

!------------------------------------------------------------
! input_fline
!       - checks if a file exists, and if so, how many lines are in it
!------------------------------------------------------------

SUBROUTINE input_fline(fline,fname,error)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(IN) :: fname
  INTEGER, INTENT(INOUT) :: fline,error
  !Internal
  INTEGER :: io
  LOGICAL :: ex
  error = 0
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You need to create the input file : ", TRIM(fname)
    error = 1
    fline = -1
    RETURN
  END IF
  fline = 0
  io = 0
  OPEN(unit=999,file=TRIM(fname),status='old',access='sequential')
  DO WHILE (io .EQ. 0)
    READ(999,*,iostat=io)
    IF (io .EQ. 0) fline = fline + 1
  END DO
  CLOSE(unit=999)

END SUBROUTINE input_fline
!------------------------------------------------------------

END MODULE input
