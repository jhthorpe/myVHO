MODULE val

CONTAINS

SUBROUTINE checkval(val,error)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: val
  LOGICAL, INTENT(INOUT) :: error
  REAL(KIND=8) :: s,infty
  infty = HUGE(s)
  error = .FALSE.
  IF (val .GT. infty) THEN
    WRITE(*,*) 
    WRITE(*,*) "val:checkval -- val is infinty" 
    error = .TRUE.
  ELSE IF (val .NE. val) THEN
    WRITE(*,*)
    WRITE(*,*) "val:checkval -- val is NaN" 
    error = .TRUE.
  END IF
END SUBROUTINE

END MODULE val
