!------------------------------------------------------------
!  quad
!       - module containing subroutines for dealing with 
!         quadratic force constants             
!------------------------------------------------------------
MODULE quad
  USE input

CONTAINS

!------------------------------------------------------------
! quad_get
!       -reads in quadratic force constants and labeling
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational numbering offset
! phi2          : 1D real*8, 2nd order force constants
! h2l           : 1D int, converts harmonic levels -> labels
! l2h           : 1D int, converts labels -> harmonic levels
! error         : int, exit code

SUBROUTINE quad_get(nvib,voff,phi2,h2l,l2h,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: phi2
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: h2l,l2h
  INTEGER, INTENT(INOUT) :: nvib,voff,error

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Rtemp
  INTEGER, DIMENSION(:), ALLOCATABLE :: Itemp
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,fid,fline,vmax
  LOGICAL :: ex,found

  error = 0
  fname = 'quadratic'
  fid = 100  

  !Check file exits
  INQUIRE(FILE=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "quad_get  : ERROR"
    WRITE(*,*) "Could not find the 'quadratic' file"
    error = error + 1
    RETURN
  END IF
 
  !Get number of lines
  CALL input_fline(fline,fname,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "quad_get  : ERROR"
    WRITE(*,*) "There is some problem in",TRIM(fname)
    RETURN
  END IF
  ALLOCATE(Itemp(0:fline-1))
  ALLOCATE(Rtemp(0:fline-1))
  Itemp = 0
  Rtemp = 0.0D0
  nvib = 0
  vmax = 0

  !Read Force Constants
  OPEN(FILE=TRIM(fname),UNIT=fid,STATUS='old')
  DO j=0,fline-1
    READ(fid,*) i,val
    Itemp(j) = i
    Rtemp(j) = val
    !update voff
    IF (j .EQ. 0) voff = i
    IF (i .LT. voff) voff = i
    IF (i .GT. vmax) vmax = i 
    !update number of unique vibrational modes
    IF (ALL(i .NE. Itemp(0:j-1))) nvib = nvib + 1
  END DO
  CLOSE(UNIT=fid)

  !Put unique values into place and check
  ALLOCATE(phi2(0:nvib-1))
  ALLOCATE(h2l(0:nvib-1))
  DO j=0,fline-1
    IF (Itemp(j)-voff .LT. 0 .OR. Itemp(j)-voff .GT. nvib-1) THEN
      WRITE(*,*) "quad_get  : ERROR"
      WRITE(*,*) "There is an issue with the mode numbering in 'quadratic'"
      WRITE(*,*) Itemp(j),Rtemp(j)
      WRITE(*,*) "^^^ something is out of bounds"
      error = error + 1
      RETURN
    END IF
    phi2(Itemp(j)-voff) = Rtemp(j) 
  END DO

  !read h2l
  fname = 'modenumber'
  INQUIRE(FILE=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    OPEN(FILE=TRIM(fname),UNIT=fid,STATUS='old')
    DO j=0,nvib-1
      READ(fid,*) i,k 
      h2l(i-voff) = k
    END DO
    CLOSE(UNIT=fid)
  ELSE
    h2l = (/ (i,i=1,nvib) /) 
  END IF
  !check the values
  DO j=0,nvib-1
    IF (h2l(j) .LT. 1 .OR. h2l(j) .GT. nvib) THEN
      WRITE(*,*) "quad_get  : ERROR"
      WRITE(*,*) "There is a problem with the following label"
      WRITE(*,*) j+1,h2l(j)
      error = error + 1
    END IF
  END DO
  IF (error .NE. 0) RETURN

  !generate l2h
  ALLOCATE(l2h(0:nvib-1))
  DO i=0,nvib-1
    l2h(h2l(i)-1) = i+1 
  END DO

!  WRITE(*,*) "TESTING TESTING TESTING"
!  DO WHILE(.TRUE.) 
!    READ(*,*) i
!    WRITE(*,*) l2h(i-1)
!  END DO

  !Print Quadratic Force Constants
  WRITE(*,*) "Harmonic Labeling"
  DO i=0,nvib-1
    WRITE(*,'(2x,I2,2x,I2,2x,F18.13)') i+1,h2l(i),phi2(i)
  END DO
  WRITE(*,*)

  DEALLOCATE(Itemp)
  DEALLOCATE(Rtemp)

END SUBROUTINE quad_get
!------------------------------------------------------------


END MODULE quad
!------------------------------------------------------------
