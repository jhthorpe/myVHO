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
! enum          : int, nubmer of eigenvalues to compute
! mem           : int*8, memory in MB
! error         : int, error code

SUBROUTINE input_jobinfo(job,ndim,nbas,enum,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: nbas
  INTEGER(KIND=8), INTENT(INOUT) :: mem
  INTEGER, INTENT(INOUT) :: job,ndim,enum,error

  error = 0
  WRITE(*,*) "input_jobinfo : called"
  WRITE(*,*) "Reading input from vho.in"

  CALL input_read(job,ndim,nbas,enum,mem,error)
  IF (error .NE. 0) RETURN
  CALL input_check(job,ndim,nbas,enum,mem,error)
  IF (error .NE. 0) RETURN
  CALL input_write(job,ndim,nbas,enum,mem,error)

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
! enum          : int, nubmer of eigenvalues to compute
! mem           : int*8, memory in MB

SUBROUTINE input_read(job,ndim,nbas,enum,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: nbas
  INTEGER(KIND=8), INTENT(INOUT) :: mem
  INTEGER, INTENT(INOUT) :: job,ndim,enum,error
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
  READ(100,*) enum
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
! enum          : int, nubmer of eigenvalues to compute
! mem           : int*8, memory in MB

SUBROUTINE input_check(job,ndim,nbas,enum,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,enum
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
  IF (enum .LT. 1 ) THEN
    WRITE(*,*) "vho.in line #4"
    WRITE(*,*) "Need at least one eigenvalue"
  END IF
  IF (mem .LT. 1) THEN
    WRITE(*,*) "vho line #5"
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
! enum          : int, nubmer of eigenvalues to compute
! mem           : int*8, memory in MB
! error         : int, error code

SUBROUTINE input_write(job,ndim,nbas,enum,mem,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,enum
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
  WRITE(*,*) "enum    :",enum
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
! input_nline
!       - checks how many lines are in a file 
!       - if file does not exist, returns 0
!------------------------------------------------------------

SUBROUTINE input_nline(nline,fname)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(IN) :: fname
  INTEGER, INTENT(INOUT) :: nline
  !Internal
  INTEGER :: io
  LOGICAL :: ex
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    nline = 0
    RETURN
  END IF
  nline = 0
  io = 0
  OPEN(unit=999,file=TRIM(fname),status='old',access='sequential')
  DO WHILE (io .EQ. 0)
    READ(999,*,iostat=io)
    IF (io .EQ. 0) nline = nline + 1
  END DO
  CLOSE(unit=999)

END SUBROUTINE input_nline

!------------------------------------------------------------
! input_QUAD_natoms
!       - gets number of atoms in input file
!------------------------------------------------------------
! natoms        : int, number of atoms
! error         : int, exit code

SUBROUTINE input_QUAD_natoms(natoms,error)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: natoms,error
  CHARACTER(LEN=1024) :: fname
  CHARACTER(LEN=4) :: dummy, line
  INTEGER :: fid
  error = 0
  fname = "QUADRATURE"
  fid = 99
  line = "    " 
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO WHILE( line .NE. "back")
    READ(fid,*) dummy, line
    IF (line .EQ. "freq") READ(fid,*)
  END DO
  DO WHILE (line .NE. "freq") 
    natoms = natoms + 1
    READ(fid,*) dummy, line
  END DO
  CLOSE(unit=fid)

  WRITE(*,*) "Number of atoms: ", natoms

END SUBROUTINE input_QUAD_natoms

!------------------------------------------------------------
! input_QUAD_q0
!       - reads undisplaced geometry from QUADRATURE file
!------------------------------------------------------------
! natoms        : int, number of atoms
! q0            : 2D real*8, undisplaced geometry   [atom,xyz]
! error         : int, error code

SUBROUTINE input_QUAD_q0(natoms,q0,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: q0
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: natoms
  REAL(KIND=8), DIMENSION(0:2) :: xyz
  CHARACTER(LEN=1024) :: fname
  CHARACTER(LEN=9) :: dummy,line
  INTEGER :: fid,i
  error = 0
  fname = "QUADRATURE"
  fid = 99
  line = "        "
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO WHILE (line .NE. "Reference")
    READ(fid,*) dummy, line
  END DO
  DO i=0,natoms-1
    READ(fid,*) xyz 
    q0(i,0:2) = xyz
  END DO
  CLOSE(unit=fid)
  
  WRITE(*,*) "Undisplaced coordinates" 
  DO i=0,natoms-1
    WRITE(*,'(1x,3(f14.10,2x))') q0(i,0:2)
  END DO
  WRITE(*,*)

END SUBROUTINE input_QUAD_q0

!------------------------------------------------------------
! input_QUAD_bas
!       - gets the basis function frequencies 
!------------------------------------------------------------
! ndim          : int, number of dimensions 
! bask          : 1D real*8, basis function force constants
! error         : int, exit code

SUBROUTINE input_QUAD_bask(ndim,bask,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: bask
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  CHARACTER(LEN=4) :: dummy, line
  INTEGER :: fid,i
  error = 0
  fname = "QUADRATURE"
  fid = 99
  line = "    " 
  bask = 0
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO i=0,ndim-1
    DO WHILE( line .NE. "freq")
      READ(fid,*) dummy, line
    END DO
    READ(fid,*) bask(i)
    line = "    "
  END DO
  CLOSE(unit=fid)
  WRITE(*,*) "Writing basis force constants to basis.in"
  OPEN(file='basis.in',unit=100,status='replace')
  DO i=0,ndim-1
    WRITE(100,*) i+1,bask(i)
  END DO
  CLOSE(unit=100)
END SUBROUTINE input_QUAD_bask

!------------------------------------------------------------
! input_QUAD_qmat
!       - read normal coordinate matricies from QUADRATURE
!------------------------------------------------------------
! ndim          : int, number of dimensions 
! natoms        : int, number of atoms
! qmat          : 3D real*8,q matricies, [atom,xyz,dim]
! error         : int, exit code

SUBROUTINE input_QUAD_qmat(ndim,natoms,qmat,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: qmat
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,natoms
  CHARACTER(LEN=1024) :: fname
  CHARACTER(LEN=4) :: dummy, line
  INTEGER :: fid,i,j
  error = 0
  fname = "QUADRATURE"
  fid = 99
  line = "    " 
  qmat = 0
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO i=0,ndim-1
    DO WHILE( line .NE. "back")
      READ(fid,*) dummy, line
      IF (line .EQ. "freq") READ(fid,*)
    END DO
    DO j=0,natoms-1
      READ(fid,*) qmat(j,0:2,i)
    END DO
    line = "    "
  END DO
  CLOSE(unit=fid)
END SUBROUTINE input_QUAD_qmat
!------------------------------------------------------------

END MODULE input
