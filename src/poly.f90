!//////////////////////////////////////////////////////////////////
!///
!///            Module that contains subroutines for 
!///            fitting and evaluating polynomial functions
!///            to an arbitrary order
!///
!//////////////////////////////////////////////////////////////////

MODULE poly
  USE linal
  USE val

CONTAINS
!---------------------------------------------------------------------
! poly_fit
!       - subrotuine that fits a function to an increasing number of
!         polynomials, until a threshold RMS is met
!---------------------------------------------------------------------
! Variables
! Fx            : 1D real*8, function to fit
! x             : 1D real*8, list of x values
! N             : int, number of points
! tol           : real*8, RMS target tollerance
! ord           : int, on exit, order of the polynomial used 
! coef          : 1D real*8, list of coefficients used
! error         : bool, true on exit if there was a problem

SUBROUTINE poly_fit(Fx,x,N,tol,ord,coef,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: coef
  REAL(KIND=8), DIMENSION(0:N-1), INTENT(IN) :: Fx,x
  REAL(KIND=8), INTENT(IN) :: tol
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(INOUT) :: ord
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: b
  REAL(KIND=8) :: rms
  INTEGER :: i,j

  error = .FALSE.

  WRITE(*,*) 
  WRITE(*,*) "Generating polynomial fit of surface"

  ALLOCATE(A(0:N-1,0:N-1))
  ALLOCATE(b(0:N-1))
  A = 0
  b = 0

  !loop until we get below tolerance
  DO i=0,N-1 
    ord = i
    A(ord,0:N-1) = (/ (x(j)**ord,j=0,N-1) /)  !add on each row

    CALL lsqr(A(0:ord,0:N-1),Fx(0:N-1),N,ord+1,b(0:ord),error) 

    IF (error) THEN
      WRITE(*,*) "poly:poly_fit -- lsqr exited with error" 
      error = .TRUE. 
      GOTO 11
    END IF 

    CALL poly_rms(Fx,x,N,coef,ord,rms,error)
    WRITE(*,*) "At iteration ", ord," rms is", rms

    IF (rms .LT. tol) EXIT
  END DO

  IF (rms .GT. tol) THEN
    error = .TRUE.
    WRITE(*,*) "poly:poly_fit -- fail to get below RMS target"
    DEALLOCATE(coef)
    RETURN
  ELSE
    WRITE(*,'(2x,A17,I4,A11)') "RMS converged in ", ord, " iterations"
    OPEN(file='fit.dat',unit=105,status='replace')
    WRITE(105,*) ord
    DO i=0,ord
      WRITE(105,*) i, coef(i) 
    END DO
    CLOSE(unit=105)
    WRITE(*,*) "Coefficients written to fit.dat"
  END IF 

  ALLOCATE(coef(0:ord))
  coef(0:ord) = b(0:ord)

11 DEALLOCATE(A)
   DEALLOCATE(b)

END SUBROUTINE poly_fit

!---------------------------------------------------------------------
! poly_rms
!       - returns rms of poynomial fit
!---------------------------------------------------------------------
SUBROUTINE poly_rms(Fx,x,N,coef,ord,rms,error)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Fx,x,coef
  REAL(KIND=8), INTENT(INOUT) :: rms
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,ord

  REAL(KIND=8) :: val
  INTEGER :: i

  error = .FALSE.

  rms = 0.0D0
  DO i=0,N-1
    CALL poly_eval(ord,coef,x(i),val,error) 
    IF (error) RETURN
    rms = rms + (Fx(i) - val)**2.0
  END DO

END SUBROUTINE poly_rms

!---------------------------------------------------------------------
! poly_eval
!       - given order and coefficients, returns the value of 
!         a fitted polynomial at x
!---------------------------------------------------------------------
SUBROUTINE poly_eval(ord,coef,x,val,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: coef
  REAL(KIND=8), INTENT(INOUT) :: val
  REAL(KIND=8), INTENT(IN) :: x
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ord

  REAL(KIND=8) :: t
  INTEGER :: i
  
  error = .FALSE.

  val = 0
  DO i=0,ord !remember that ord = 1 is a 1st order polynomial, a*x + b
    val = val + coef(i)*x**i
  END DO

  !check for inf/nan
  CALL checkval(val,error)
  
  IF (error) THEN
    WRITE(*,*) "poly:poly_eval -- bad value in evaulation" 
  END IF

END SUBROUTINE poly_eval

!---------------------------------------------------------------------

END MODULE poly
