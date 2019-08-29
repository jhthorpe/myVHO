!------------------------------------------------------------
! cori
!       - module containing subroutines involving the 
!         Coriolis zetas
!------------------------------------------------------------
MODULE cori
  USE input

CONTAINS

!------------------------------------------------------------
! cori_get
!       - reads in, sorts, and stores the coriolis coupling
!         constants as a 3 index matrix
!       - the qantum numbers of the constants are stored
!         [i,j,alpha], all indexed from zero
!------------------------------------------------------------
! nvib          : int, number of dimensions
! voff          : int, vibrational numbering offset
! cori          : 2D real*8, coriolis zetas [zeta,rotation]
! error         : int, exit code

SUBROUTINE cori_get(nvib,voff,cori,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: cori
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nvib,voff
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  LOGICAL :: ex
  INTEGER :: fline,fid
  INTEGER :: i,j,a,k

  error = 0
  fname = 'coriolis'
  fid = 300

  ALLOCATE(cori(0:nvib-1,0:nvib-1,0:2))
  cori = 0.0D0

  !If only one dimension
  INQUIRE(file=TRIM(fname),exist=ex)
  IF (nvib .EQ. 1 .OR. .NOT. ex) THEN
    RETURN
  END IF

  !Initialize
  CALL input_nline(fline,fname)

  !read in data
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO k=0,fline-1
    READ(fid,*) a,i,j,val
    i = i - voff
    j = j - voff
    a = a - 1
    IF (a .LT. 0 .OR. a .GT. 2 .OR. i .LT. 0 .OR. i .GT. nvib-1 .OR. &
      j .LT. 0 .OR. j .GT. nvib-1) THEN
      WRITE(*,*) "cori_get  : ERROR"
      WRITE(*,*) "In coriolis, line ",k,", has a bad value"
      error = error + 1
    END IF
    cori(i,j,a) = val
  END DO
  CLOSE(unit=fid)
  IF (error .NE. 0) RETURN

!  DO a=0,2
!    DO i=0,nvib-1
!      DO j=i+1,nvib-1
!        IF (ABS(cori(i,j,a)) .GT. 1.0D-15) THEN
!         WRITE(*,'(2x,I1,2x,2(I4,2x),F24.15)') a+1,i+voff,j+voff,cori(i,j,a)
!         WRITE(*,'(2x,I1,2x,2(I4,2x),F24.15)')a+1,j+voff,i+voff,cori(j,i,a)
!        END IF
!      END DO
!    END DO
!  END DO
!  WRITE(*,*)

END SUBROUTINE cori_get

!------------------------------------------------------------
END MODULE cori
!------------------------------------------------------------
