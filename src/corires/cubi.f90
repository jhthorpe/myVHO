!------------------------------------------------------------
! cubi
!       - module containing subroutines for dealing with 
!         cubic force constants             
!------------------------------------------------------------
MODULE cubi
  USE input

CONTAINS

!------------------------------------------------------------
! cubi_read
!       -reads in cubic force constants and labeling
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational numbering offset
! phi3          : 1D real*8, 2nd order force constants
! error         : int, exit code

SUBROUTINE cubi_get(nvib,voff,phi3,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: phi3
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nvib,voff

  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,l,fid,fline
  LOGICAL :: ex

  error = 0
  fname = 'cubic'
  fid = 101  

  !Check file exits
  INQUIRE(FILE=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "cubi_get  : ERROR"
    WRITE(*,*) "Could not find the 'cubic' file"
    error = error + 1
    RETURN
  END IF
 
  !Get number of lines
  CALL input_fline(fline,fname,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "cubi_get  : ERROR"
    WRITE(*,*) "There is some problem in",TRIM(fname)
    RETURN
  END IF
  ALLOCATE(phi3(0:nvib-1,0:nvib-1,0:nvib-1))
  phi3 = 0.0D0

  !Read Force Constants
  OPEN(FILE=TRIM(fname),UNIT=fid,STATUS='old')
  DO l=0,fline-1
    READ(fid,*) i,j,k,val
    i = i - voff
    j = j - voff
    k = k - voff
    IF (ANY([i,j,k] .LT. 0) .OR. ANY([i,j,k] .GT. nvib-1)) THEN
      WRITE(*,*) "cubi_get  : ERROR"
      WRITE(*,*) "There is an issue with the mode numbering in 'cubic'"
      WRITE(*,*) i+voff,j+voff,k+voff,val 
      WRITE(*,*) "^^^ something is out of bounds"
      error = error + 1
    END IF
    phi3(i,j,k) = val
    phi3(i,k,j) = val
    phi3(j,i,k) = val
    phi3(j,k,i) = val
    phi3(k,i,j) = val
    phi3(k,j,i) = val
  END DO
  CLOSE(UNIT=fid)
  IF (error .NE. 0) RETURN

  !Print Cubic Force Constants
!  WRITE(*,*) "Cubic Force Constants"
!  DO i=0,nvib-1
!     DO j=0,nvib-1
!       DO k=0,nvib-1
!        WRITE(*,'(2x,3(I2,1x),2x,F24.13)') i+voff,j+voff,k+voff,phi3(i,j,k)
!       END DO
!     END DO
!  END DO
!  WRITE(*,*)

END SUBROUTINE cubi_get
!------------------------------------------------------------


END MODULE cubi
!------------------------------------------------------------
