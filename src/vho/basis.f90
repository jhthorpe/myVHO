!------------------------------------------------------------
! basis 
!       - module dealing with basis information
!------------------------------------------------------------
MODULE basis
  USE input

CONTAINS
!------------------------------------------------------------
! basis_get
!       - gets basis information
!------------------------------------------------------------
! ndim          : int, number of dimensions
! basK          : 1D real*8, force constants of basis set
! error         : int, exit code

SUBROUTINE basis_get(ndim,basK,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: basK
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=2014) :: fname
  REAL(KIND=8) :: val
  INTEGER :: fline,i,j
  LOGICAL :: ex
  error = 0
  fname = "basis.in"
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "basis_get  : ERROR"
    WRITE(*,*) "You need the input file basis.in"
    error = 1
    RETURN
  ELSE
    ALLOCATE(basK(0:ndim-1))
    basK = -1.0D0
    CALL input_fline(fline,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "basis_get  : ERROR"
      WRITE(*,*) "You need the input file basis.in"
      error = 1
      RETURN
    END IF
    IF (fline .NE. ndim) THEN
      WRITE(*,*) "basis_get  : ERROR"
      WRITE(*,*) "All basis k must be provided"
      error = 1
      RETURN
    END IF
    OPEN(file=TRIM(fname),unit=100,status='old')
    DO i=0,ndim-1
      READ(100,*) j,val
      IF (j .LT. 1 .OR. j .GT. ndim .OR. val .LE. 0.0D0) THEN
        WRITE(*,*) "basis_get  : ERROR"
        WRITE(*,*) "Input ", i,", of basis.in has a bad value"
        error = 1
      END IF
      basK(j-1) = val 
    END DO
    CLOSE(unit=100)
  END IF

  IF (ANY(basK .LE. 0.0D0)) THEN
    WRITE(*,*) "basis_get  : ERROR"
    WRITE(*,*) "All basis k must be provided, and > 0"
    error = 1
  END IF

  WRITE(*,*) "basis is..."
  WRITE(*,*) basK

END SUBROUTINE basis_get
!------------------------------------------------------------

END MODULE basis
!------------------------------------------------------------
