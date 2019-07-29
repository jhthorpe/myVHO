!------------------------------------------------------------
! val
!       - modules for checking values 
!------------------------------------------------------------
MODULE val

CONTAINS

!------------------------------------------------------------
! val_check
!       - check for NaN or inf
!------------------------------------------------------------
SUBROUTINE val_check(val,error)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: val
  INTEGER, INTENT(INOUT) :: error
  REAL(KIND=8) :: s,infty
  infty = HUGE(s)
  error = 0
  IF (val .GT. infty) THEN
    WRITE(*,*) 
    WRITE(*,*) "ERROR"
    WRITE(*,*) "val_check  : val is infinty" 
    error = 42 
  ELSE IF (val .NE. val) THEN
    WRITE(*,*) 
    WRITE(*,*) "ERROR"
    WRITE(*,*) "val_check  : val is NaN" 
    error = 43 
  END IF
END SUBROUTINE val_check
!------------------------------------------------------------

END MODULE val
