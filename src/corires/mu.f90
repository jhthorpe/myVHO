!------------------------------------------------------------
! mu
!       - module containing subroutines involving the 
!         mu constants 
!------------------------------------------------------------
MODULE mu
  USE input
  
CONTAINS 

!------------------------------------------------------------
! mu_get
!       - generates mu0,mu1, and mu2
!       mu1 -> μ_{α,β}^r
!       mu2 -> μ_{α,β}^{r,s}
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational offset
! Be            : real*8, rotational constants in cm-1
! mu1           : 3D real*8, μ_{α,β}^r (vib,rot,rot)
! mu2           : 4D real*8, μ_{α,β}^{r,s} (vib,vib,rot,rot)
! error         : int, exit code

SUBROUTINE mu_get(nvib,voff,Be,mu1,mu2,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: mu2
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: mu1
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Be
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nvib,voff
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: didq
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val,cc,pi
  INTEGER :: i,j,a,b,g,fid,fline 
  LOGICAL :: ex

  fname = 'didq'
  fid = 200
  cc = 137.035999084
  pi = 3.1415926535897932

  !Check didq file exists
  INQUIRE(FILE=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "mu_get  : ERROR"
    WRITE(*,*) "Could not find the 'didq' file"
    error = error + 1
    RETURN
  END IF

  !Get number of lines
  CALL input_fline(fline,fname,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "mu_get  : ERROR"
    WRITE(*,*) "There is some problem in",TRIM(fname)
    RETURN
  END IF

  !Read didq
  ALLOCATE(didq(0:nvib-1,0:2,0:2))
  OPEN(FILE=TRIM(fname),UNIT=fid,STATUS='old')
  DO j=0,fline-1 
    READ(fid,*) a,b,i,val
    a = a - 1
    b = b -1
    i = i - voff
    IF (ANY([a,b] .LT. 0) .OR. ANY([a,b] .GT. 2) .OR. &
        i .LT. 0 .OR. i .GT. nvib-1) THEN
      WRITE(*,*) "mu_get  : ERROR"
      WRITE(*,*) "There is a bad value at line",j+1,"in the 'didq' file"
      error = error + 1
    END IF
    didq(i,a,b) = val
  END DO
  CLOSE(UNIT=fid)
  IF (error .NE. 0) RETURN

  !Account for factor of hbar^2/(2*h*c) in cfour's didq
  didq = didq/(4.0D0*pi*cc)

  !Generate orders of mu
  ALLOCATE(mu1(0:nvib-1,0:nvib-1,0:2))
  ALLOCATE(mu2(0:nvib-1,0:nvib-1,0:2,0:2))
  mu1 = 0.0D0
  mu2 = 0.0D0

  !Order 1
  DO a=0,2
    DO b=0,2
      DO i=0,nvib-1
        mu1(i,a,b) = -0.25D0*didq(i,a,b)*Be(a)*Be(b)
      END DO
    END DO
  END DO

  !Order 2
  DO a=0,2
    DO b=0,2
      DO i=0,nvib-1
        DO j=0,nvib-1
          DO g=0,2
            mu2(i,j,a,b) = mu2(i,j,a,b) + 0.125D0*Be(a)*Be(b)*Be(g)*&
                           (didq(i,a,g)*didq(j,g,b) + didq(j,a,g)*didq(i,g,b))
          END DO
        END DO
      END DO
    END DO
  END DO
  mu2 = 0.75D0*mu2

  DEALLOCATE(didq)

END SUBROUTINE mu_get
!------------------------------------------------------------

END MODULE mu
!------------------------------------------------------------
