!//////////////////////////////////////////////////////////////////
!///
!///		Module Containing subroutines for reading
!///		input for myVSCF
!///
!///		James H. Thorpe
!///
!//////////////////////////////////////////////////////////////////

MODULE input
  IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------
!	read_input
!		-reads input
!---------------------------------------------------------------------
! Variables
! nbasis	: int, number of harmonic oscillator basis functions
! Vq		: 1D real*8, 1D potential energy surface
! qmin		: real*8, minimum r
! qmax		: real*8, max r
! error		: bool, true if error

SUBROUTINE read_input(nbasis,Vq,q,qmin,qmax,npoints,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: Vq,q
  REAL(KIND=8), INTENT(INOUT) :: qmin,qmax
  INTEGER, INTENT(INOUT) :: nbasis, npoints
  LOGICAL, INTENT(INOUT)  :: error
  
  CHARACTER(LEN=1024) :: fname 
  INTEGER :: dummy,i
  
  error = .FALSE.
  fname = 'Vq'
  CALL EXECUTE_COMMAND_LINE('cat input')

  OPEN(file='input',unit=100,status='old')
  READ(100,*) nbasis
  CLOSE(unit=100)

  CALL getfline(npoints,fname,error)

  ALLOCATE(Vq(0:npoints-1)) 
  ALLOCATE(q(0:npoints-1))

  OPEN(file="Vq",unit=101,status='old')
  !read first line
  READ(101,*) dummy, q(0), Vq(0)
  qmin = q(0) 

  DO i=1,npoints-2
    READ(101,*) dummy, q(i), Vq(i) 
  END DO

  !read last line
  READ(101,*) dummy, q(npoints-1), Vq(npoints-1)
  qmax = q(npoints-1)
   
  CLOSE(unit=101)

  DO i=0,npoints-1
  WRITE(*,*) q(i), Vq(i)
  END DO

END SUBROUTINE read_input 

!---------------------------------------------------------------------
! getfline
!       James H. Thorpe
!       - checks if a file exists, and if so, how many lines are in it
!---------------------------------------------------------------------

SUBROUTINE getfline(fline,fname,flag)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(IN) :: fname
  INTEGER, INTENT(INOUT) :: fline
  LOGICAL, INTENT(INOUT) :: flag
  !Internal
  INTEGER :: io
  LOGICAL :: ex
  flag = .FALSE.
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You need to create the input file : ", TRIM(fname)
    flag = .TRUE.
    fline = -1
    RETURN
  END IF
  fline = 0
  io = 0
  OPEN(unit=999,file=TRIM(fname),status='old',access='sequential')
  DO WHILE (io .EQ. 0)
    READ(1,*,iostat=io)
    IF (io .EQ. 0) fline = fline + 1
  END DO
  CLOSE(unit=999)

END SUBROUTINE getfline

!---------------------------------------------------------------------


END MODULE input
