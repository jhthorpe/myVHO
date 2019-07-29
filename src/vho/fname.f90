!---------------------------------------------------------------------
! fname
!       - module for dealing with file naming schemes
!---------------------------------------------------------------------
MODULE fname

CONTAINS

!---------------------------------------------------------------------
! fname_Vin
!       - constructs fname for potential files, Vx.in
!---------------------------------------------------------------------
! id            : int, id
! fname         : int, filename
! error         : int, exit code 

SUBROUTINE fname_Vin(id,fname,error)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(INOUT) :: fname
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: id
  !Internal
  CHARACTER(LEN=1024) :: str_fmt

  error = 0 

  IF (id .LT. 10) THEN
    str_fmt = "(A1,I1,A3)"
  ELSE IF (id .GE. 10 .AND. id .LT. 100) THEN
    str_fmt = "(A1,I2,A3)"
  ELSE IF (id .GE. 100 .AND. id .LT. 1000) THEN
    str_fmt = "(A1,I3,A3)"
  ELSE IF (id .GE. 1000 .AND. id .LT. 10000) THEN
    str_fmt = "(A1,I4,A3)"
  ELSE
    WRITE(*,*) "ERROR"
    WRITE(*,*) "fname_Vin  : need more cases"
    error = 0 
    RETURN
  END IF

  WRITE(fname,str_fmt) "V",id,".in"

END SUBROUTINE fname_Vin

!---------------------------------------------------------------------
! fname_splinedat
!       - constructs fname for potential files, splinex.dat
!---------------------------------------------------------------------
! id            : int, id
! fname         : int, filename
! error         : int, exit code 

SUBROUTINE fname_splinedat(id,fname,error)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(INOUT) :: fname
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: id
  !Internal
  CHARACTER(LEN=1024) :: str_fmt

  error = 0 

  IF (id .LT. 10) THEN
    str_fmt = "(A6,I1,A4)"
  ELSE IF (id .GE. 10 .AND. id .LT. 100) THEN
    str_fmt = "(A6,I2,A4)"
  ELSE IF (id .GE. 100 .AND. id .LT. 1000) THEN
    str_fmt = "(A6,I3,A4)"
  ELSE IF (id .GE. 1000 .AND. id .LT. 10000) THEN
    str_fmt = "(A6,I4,A4)"
  ELSE
    WRITE(*,*) "ERROR"
    WRITE(*,*) "fname_splinedat  : need more cases"
    error = 0 
    RETURN
  END IF

  WRITE(fname,str_fmt) "spline",id,".dat"

END SUBROUTINE fname_splinedat

!---------------------------------------------------------------------
! fname_pointstxt
!       - constructs fname for potential files, splinex.dat
!---------------------------------------------------------------------
! id            : int, id
! fname         : int, filename
! error         : int, exit code 

SUBROUTINE fname_pointstxt(id,fname,error)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(INOUT) :: fname
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: id
  !Internal
  CHARACTER(LEN=1024) :: str_fmt
  error = 0 
  IF (id .LT. 10) THEN
    str_fmt = "(A6,I1,A4)"
  ELSE IF (id .GE. 10 .AND. id .LT. 100) THEN
    str_fmt = "(A6,I2,A4)"
  ELSE IF (id .GE. 100 .AND. id .LT. 1000) THEN
    str_fmt = "(A6,I3,A4)"
  ELSE IF (id .GE. 1000 .AND. id .LT. 10000) THEN
    str_fmt = "(A6,I4,A4)"
  ELSE
    WRITE(*,*) "ERROR"
    WRITE(*,*) "fname_splinedat  : need more cases"
    error = 0 
    RETURN
  END IF
  WRITE(fname,str_fmt) "points",id,".txt"
END SUBROUTINE fname_pointstxt

!---------------------------------------------------------------------
END MODULE fname
!---------------------------------------------------------------------
