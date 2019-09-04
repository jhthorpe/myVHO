!------------------------------------------------------------
! quar
!       - module containing subroutines for dealing with 
!         quartic force constants             
!------------------------------------------------------------
MODULE quar
  USE input

CONTAINS

!------------------------------------------------------------
! quar_read
!       -reads in quartic force constants and labeling
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational numbering offset
! phi4          : 4D real*8, fourth order force constants
! error         : int, exit code

SUBROUTINE quar_get(nvib,voff,phi4,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: phi4
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nvib,voff

  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,l,n,fid,fline
  LOGICAL :: ex

  error = 0
  fname = 'quartic'
  fid = 101  

  !Check file exits
  INQUIRE(FILE=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "quar_get  : ERROR"
    WRITE(*,*) "Could not find the 'quar' file"
    error = error + 1
    RETURN
  END IF
 
  !Get number of lines
  CALL input_fline(fline,fname,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "quar_get  : ERROR"
    WRITE(*,*) "There is some problem in",TRIM(fname)
    RETURN
  END IF
  ALLOCATE(phi4(0:nvib-1,0:nvib-1,0:nvib-1,0:nvib-1))
  phi4 = 0.0D0

  !Read Force Constants
  OPEN(FILE=TRIM(fname),UNIT=fid,STATUS='old')
  DO n=0,fline-1
    READ(fid,*) i,j,k,l,val
    i = i - voff
    j = j - voff
    k = k - voff
    l = l - voff
    IF (ANY([i,j,k,l] .LT. 0) .OR. ANY([i,j,k,l] .GT. nvib-1)) THEN
      WRITE(*,*) "quar_get  : ERROR"
      WRITE(*,*) "There is an issue with the mode numbering in 'quartic'"
      WRITE(*,*) i+voff,j+voff,k+voff,l+voff,val 
      WRITE(*,*) "^^^ something is out of bounds"
      error = error + 1
    END IF
    phi4(i, j, k, l) = val
    phi4(i, j, l, k) = val
    phi4(i, k, j, l) = val
    phi4(i, k, l, j) = val
    phi4(i, l, j, k) = val
    phi4(i, l, k, j) = val
    phi4(j, i, k, l) = val
    phi4(j, i, l, k) = val
    phi4(j, k, i, l) = val
    phi4(j, k, l, i) = val
    phi4(j, l, i, k) = val
    phi4(j, l, k, i) = val
    phi4(k, i, j, l) = val
    phi4(k, i, l, j) = val
    phi4(k, j, i, l) = val
    phi4(k, j, l, i) = val
    phi4(k, l, i, j) = val
    phi4(k, l, j, i) = val
    phi4(l, i, j, k) = val
    phi4(l, i, k, j) = val
    phi4(l, j, i, k) = val
    phi4(l, j, k, i) = val
    phi4(l, k, i, j) = val
    phi4(l, k, j, i) = val
  END DO
  CLOSE(UNIT=fid)
  IF (error .NE. 0) RETURN

END SUBROUTINE quar_get
!------------------------------------------------------------


END MODULE quar
!------------------------------------------------------------
